#!/usr/bin/env python3
# coding: utf-8

"""
Comprehensive RNA-seq Differential Expression Analysis Pipeline

This script performs complete differential expression analysis using PyDESeq2,
with beautiful visualizations and robust statistical testing. The pipeline 
transforms raw sequencing counts into biological insights through:

1. Quality control assessment
2. Statistical differential expression analysis  
3. Rich data visualizations
4. Results export and summarization

Author: Your Name
Usage: python code.py
Dependencies: pandas, numpy, seaborn, matplotlib, scanpy, pydeseq2, gseapy, adjustText
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import random  # Used for sophisticated volcano plot aesthetics

# Statistical analysis libraries
from pydeseq2.dds import DeseqDataSet  # Core DESeq2 functionality
from pydeseq2.ds import DeseqStats     # Statistical results extraction

# Optional: Gene Set Enrichment Analysis
import gseapy as gp
from gseapy.plot import gseaplot

# Optional: Automatic text label positioning in plots
from adjustText import adjust_text

# Set random seed for reproducible visualizations
random.seed(42)
np.random.seed(42)

#######################################################
# 1. EXPERIMENTAL DESIGN: Sample Information Setup
#######################################################
"""
Define your experimental setup here. This is the foundation of your analysis
- each sample gets assigned to an experimental condition.

CUSTOMIZE THIS SECTION for your experiment:
- Update sample names to match your count matrix columns
- Update conditions to reflect your experimental design
"""

sample_info_dict = {
    "Sample": ["Kelley_17", "Kelley_18", "Kelley_19", "Kelley_20",
               "Kelley_21", "Kelley_22", "Kelley_23", "Kelley_24"],
    "Condition": ["C", "C", "C", "C",                    # Control samples
                  "mut", "mut", "mut", "mut"]            # Mutant/treatment samples
}

# Convert to pandas DataFrame for easier manipulation
sample_info = pd.DataFrame(sample_info_dict)
print("EXPERIMENTAL DESIGN:")
print("=" * 50)
print(sample_info)
print("=" * 50, "\n")

# Automatically identify sample groups for downstream analysis
ctrl_list = sample_info.loc[sample_info['Condition'] == 'C', 'Sample'].tolist()
variant_list = sample_info.loc[sample_info['Condition'] == 'mut', 'Sample'].tolist()

print(f"📋 Control group samples ({len(ctrl_list)}): {ctrl_list}")
print(f"📋 Mutant group samples ({len(variant_list)}): {variant_list}")
print(f"📊 Total samples in analysis: {len(ctrl_list) + len(variant_list)}\n")

#######################################################
# 2. DATA LOADING: Import and Clean Count Matrix
#######################################################
"""
Load the raw RNA-seq count matrix and perform initial data cleaning.
The count matrix should have genes as rows and samples as columns.
"""

print("🔄 Loading RNA-seq count data...")

# Load count matrix from CSV file
read_df = pd.read_csv('read_counts_table.csv')
print(f"📁 Initial count matrix shape: {read_df.shape}")
print(f"📝 Available columns: {list(read_df.columns)}")

# Set gene names as row indices (genes = rows, samples = columns)
read_df = read_df.set_index('gname')
print(f"🧬 Using 'gname' column as gene identifiers")

# Remove genes with zero counts across ALL samples
# These genes provide no information for differential analysis
genes_before = len(read_df)
read_df = read_df[read_df.sum(axis=1) > 0]
genes_removed = genes_before - len(read_df)
print(f"🗑️ Removed {genes_removed} genes with zero counts across all samples")

#######################################################
# 3. DATA FILTERING: Focus on Experimental Samples
#######################################################
"""
Filter the count matrix to include only samples relevant to our analysis
and prepare metadata for statistical testing.
"""

print("\n🎯 Filtering data for experimental samples...")

# Keep only the samples defined in our experimental design
all_samples = ctrl_list + variant_list
read_df = read_df[all_samples]
print(f"📊 Filtered count matrix dimensions: {read_df.shape}")
print(f"✅ Final sample list: {list(read_df.columns)}")

# Prepare sample metadata for DESeq2 requirements
sample_info.set_index('Sample', inplace=True)  # Use sample names as index
sample_info["Condition"] = sample_info["Condition"].astype("category")  # Convert to categorical
# Order categories with control first (important for fold-change interpretation)
sample_info["Condition"] = sample_info["Condition"].cat.reorder_categories(["C", "mut"], ordered=True)

print("✅ Sample metadata prepared for statistical analysis\n")

#######################################################
# QUALITY CONTROL: Sample Health Assessment  
#######################################################
"""
Before diving into differential expression, let's check if our samples
are healthy and comparable. This early QC can save hours of debugging later!
"""

print("🔍 Performing quality control assessment...")

# Calculate key quality metrics for each sample
qc_metrics = pd.DataFrame({
    'Total_Reads': read_df.sum(axis=0),        # Total sequencing depth per sample
    'Detected_Genes': (read_df > 0).sum(axis=0),  # Number of genes with >0 reads
    'Condition': sample_info['Condition']       # Experimental condition for comparison
})

print("📊 Quality Control Metrics:")
print("-" * 40)
print(qc_metrics)
print("-" * 40)

# Create informative QC visualizations
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Plot 1: Total reads per sample (sequencing depth)
sns.barplot(data=qc_metrics.reset_index(), x='index', y='Total_Reads', 
            hue='Condition', ax=axes[0])
axes[0].set_title('Total Reads per Sample\n(Sequencing Depth)')
axes[0].set_xlabel('Sample')
axes[0].set_ylabel('Total Read Counts')
axes[0].tick_params(axis='x', rotation=45)

# Plot 2: Detected genes per sample (library complexity)
sns.barplot(data=qc_metrics.reset_index(), x='index', y='Detected_Genes', 
            hue='Condition', ax=axes[1])
axes[1].set_title('Detected Genes per Sample\n(Library Complexity)')
axes[1].set_xlabel('Sample')
axes[1].set_ylabel('Number of Detected Genes')
axes[1].tick_params(axis='x', rotation=45)

# Save and display QC plots
plt.tight_layout()
plt.savefig("qc_metrics.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.show()

print("💾 QC metrics plot saved as 'qc_metrics.png'")
print("✅ Quality control assessment complete!\n")

#######################################################
# 4. STATISTICAL ANALYSIS: DESeq2 Differential Expression
#######################################################
"""
The core statistical analysis using DESeq2 algorithm. This sophisticated
method accounts for count data characteristics and provides robust
differential expression testing.
"""

print("🧮 Setting up DESeq2 statistical analysis...")

# DESeq2 expects samples as rows and genes as columns (opposite of typical format)
transposed_df = read_df.T
print(f"🔄 Transposed data shape: {transposed_df.shape} (samples × genes)")

# Create DESeqDataSet object with count data, sample info, and experimental design
my_dds = DeseqDataSet(
    counts=transposed_df,           # Count matrix (samples × genes)
    metadata=sample_info,           # Sample information with conditions
    design_factors=["Condition"]    # Which factor(s) to test for differences
)

print("📋 DESeq2 dataset created:")
print(my_dds)

print("\n🚀 Running DESeq2 analysis (this may take a few moments)...")

# Run the complete DESeq2 pipeline:
# 1. Estimate size factors (normalize for library size differences)
# 2. Estimate dispersions (model count variability)  
# 3. Fit negative binomial GLM and test for differences
my_dds.deseq2()

print("✅ DESeq2 analysis completed!")
print(my_dds)

#######################################################
# 5. RESULTS EXTRACTION: Getting the Statistical Evidence
#######################################################
"""
Extract the statistical results from DESeq2 analysis.
We're asking: "Which genes are significantly different between mutant and control?"
"""

print("\n📊 Extracting differential expression results...")

# Define the statistical contrast: mutant vs control
# This tells DESeq2 exactly which comparison we want
my_stats = DeseqStats(my_dds, contrast=["Condition", "mut", "C"])
stats_summary = my_stats.summary()  # Get summary statistics
results_df = my_stats.results_df    # Get detailed results for each gene

print("📈 Statistical Summary:")
print(stats_summary)
print("\n🔍 First few results:")
print(results_df.head())

#######################################################
# 6. GENE ANNOTATION: Adding Biological Context
#######################################################
"""
Transform cryptic gene IDs into meaningful gene names and descriptions.
This makes results much more interpretable and biologically relevant.
"""

print("\n📚 Adding gene annotations for biological context...")

try:
    # Load gene annotation file
    gene_info = pd.read_csv('Arabidopsis_gene_annotation.tsv', sep='\t')
    gene_info.set_index('Nomenclature ID', inplace=True)
    
    # Filter annotation file to include only genes with symbols
    # Note: This filtering step ensures we only keep well-annotated genes
    gene_info['check_symbol'] = gene_info['Symbol'].isnull()  # Fixed from original
    gene_info = gene_info.loc[gene_info['check_symbol'] == False]
    
    # Merge statistical results with gene annotations
    complete_df = pd.merge(
        left=results_df,        # Statistical results
        right=gene_info,        # Gene annotations
        how="inner",           # Only keep genes present in both datasets
        left_index=True,       # Match on gene IDs
        right_index=True
    )
    
    # Add gene symbols to main results dataframe
    results_df['Symbol'] = complete_df['Symbol']
    print(f"✅ Added gene symbols for {len(complete_df)} genes")
    
except FileNotFoundError:
    print("⚠️ Gene annotation file not found - using gene IDs as symbols")
    results_df['Symbol'] = results_df.index
except Exception as e:
    print(f"⚠️ Error loading annotations: {e}")
    results_df['Symbol'] = results_df.index

#######################################################
# 7. SIGNIFICANCE FILTERING: Finding Meaningful Changes
#######################################################
"""
Apply filters to identify genes with both statistical significance
and biological relevance. We use two criteria:
1. Statistical: adjusted p-value < 0.1 (10% false discovery rate)
2. Biological: |log2FoldChange| > 0.5 (at least 1.4-fold change)
"""

print("\n🎯 Identifying significantly changed genes...")

# Filter out lowly expressed genes (these are unreliable for statistical testing)
genes_before_filter = len(results_df)
results_df = results_df[results_df.baseMean >= 10]
genes_after_filter = len(results_df)
print(f"🔧 Filtered out {genes_before_filter - genes_after_filter} lowly expressed genes (baseMean < 10)")

# Apply significance criteria to identify truly interesting genes
significant_hits = results_df[
    (results_df.padj < 0.1) &                    # Statistically significant (FDR < 10%)
    (abs(results_df.log2FoldChange) > 0.5)       # Biologically meaningful (>1.4-fold change)
]

print(f"\n🎉 DISCOVERY SUMMARY:")
print(f"📊 Total genes analyzed: {len(results_df)}")
print(f"⭐ Significantly changed genes: {len(significant_hits)}")

if len(significant_hits) > 0:
    upregulated = (significant_hits.log2FoldChange > 0).sum()
    downregulated = (significant_hits.log2FoldChange < 0).sum()
    print(f"📈 Upregulated in mutant: {upregulated}")
    print(f"📉 Downregulated in mutant: {downregulated}")
    print("\n🔥 Most dramatic changes:")
    print(significant_hits.nlargest(5, 'log2FoldChange')[['Symbol', 'log2FoldChange', 'padj']])
else:
    print("❌ No genes met significance criteria - consider adjusting thresholds")

#######################################################
# 8. DATA EXPLORATION: Visualizing the Big Picture
#######################################################
"""
Create comprehensive visualizations to explore overall data patterns.
These plots help you understand sample relationships and identify
potential issues before focusing on individual genes.
"""

print("\n🎨 Creating exploratory data visualizations...")

try:
    # Access normalized count data from DESeq2 results
    # Different PyDESeq2 versions may store this differently
    norm_data = my_dds.layers["normed_counts"]
    print("✅ Successfully accessed normalized counts")
    
except AttributeError:
    print("⚠️ Could not access 'my_dds.layers[\"normed_counts\"]'")
    print("   This might be a PyDESeq2 version issue - check documentation")
    norm_data = None

# Only proceed with visualizations if we have normalized data
if norm_data is not None:
    
    # BUILD SCANPY OBJECT FOR ADVANCED ANALYSIS
    # AnnData is the standard format for single-cell analysis tools
    adata = anndata.AnnData(
        X=norm_data,                                    # Normalized expression matrix
        obs=my_dds.obs,                                # Sample metadata  
        var=pd.DataFrame(index=my_dds.var.index)       # Gene metadata
    )
    
    # Add additional data layers for different transformations
    adata.layers["counts"] = norm_data              # Store original normalized counts
    adata.layers["log1p"] = np.log1p(norm_data)     # Log-transformed counts
    
    print("🔬 Created AnnData object for advanced analysis")

    # PRINCIPAL COMPONENT ANALYSIS (PCA)
    # PCA reveals the major sources of variation in your data
    print("🎯 Performing Principal Component Analysis...")
    
    sc.pp.scale(adata, max_value=10)    # Scale data to unit variance (prevents high-expression genes from dominating)
    sc.tl.pca(adata)                    # Compute principal components
    
    # Create PCA plot colored by experimental condition
    sc.pl.pca(adata, color='Condition', size=150, show=False)
    plt.title("PCA: Sample Relationships Across Transcriptome\n(Samples should cluster by condition)")
    plt.savefig("pca_plot.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    print("💾 PCA plot saved as 'pca_plot.png'")

    # SAMPLE CORRELATION HEATMAP
    # Shows how similar each sample is to every other sample
    print("🔥 Creating sample correlation heatmap...")
    
    corr_mat = np.corrcoef(norm_data)  # Calculate pairwise correlations between samples
    
    # Create clustered heatmap with sample names as labels
    sns.clustermap(
        corr_mat,
        row_cluster=True, col_cluster=True,      # Cluster both rows and columns
        cmap="vlag",                             # Blue-white-red color scheme
        xticklabels=adata.obs.index,            # Sample names on x-axis
        yticklabels=adata.obs.index,            # Sample names on y-axis
        figsize=(8, 8)
    )
    plt.suptitle("Sample-to-Sample Correlation Matrix\n(Similar samples = warmer colors)", y=1.02)
    plt.savefig("correlation_heatmap.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    print("💾 Correlation heatmap saved as 'correlation_heatmap.png'")

    # MA PLOT (Mean vs. Average)
    # Shows relationship between gene expression level and fold-change magnitude
    print("📈 Creating MA plot...")
    
    plt.figure(figsize=(8, 6))
    
    # Create scatter plot with fold-change colored by magnitude
    scatter = plt.scatter(
        x=results_df['baseMean'],           # X-axis: average expression level
        y=results_df['log2FoldChange'],     # Y-axis: fold-change (mutant vs control)
        c=results_df['log2FoldChange'],     # Color by fold-change magnitude
        alpha=0.6,                          # Semi-transparent points
        cmap='RdBu_r',                      # Red-blue color scheme
        s=20                                # Point size
    )
    
    # Add reference lines
    plt.xscale('log')                                    # Log scale for expression level
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.7)     # No change line
    plt.axhline(y=0.5, color='gray', linestyle=':', alpha=0.5)     # Significance thresholds
    plt.axhline(y=-0.5, color='gray', linestyle=':', alpha=0.5)
    
    # Styling and labels
    plt.xlabel("Mean Normalized Expression (log scale)")
    plt.ylabel("Log2 Fold Change (mutant vs control)")
    plt.title("MA Plot: Expression Level vs Fold Change")
    plt.colorbar(scatter, label='Log2 Fold Change')
    plt.grid(True, alpha=0.3)
    
    plt.savefig("ma_plot.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    print("💾 MA plot saved as 'ma_plot.png'")

else:
    print("⚠️ Skipping exploratory visualizations due to normalized data access issue")

#######################################################
# 9. VOLCANO PLOT: The Crown Jewel Visualization
#######################################################
"""
Create a sophisticated volcano plot that combines statistical significance
with biological effect size. This is often the most impactful figure
in differential expression papers!

The plot uses weighted random sampling to create visually appealing
gene categories while maintaining statistical rigor.
"""

print("\n🌋 Creating the masterpiece: Volcano Plot...")

# Prepare data for volcano plot
results_df.rename(columns={"Symbol": "symbol"}, inplace=True)  # Standardize column name
results_df = results_df.dropna(subset=["padj", "symbol", "log2FoldChange", "baseMean"])  # Remove incomplete data
results_df["nlog10"] = -np.log10(results_df["padj"])  # Transform p-values for plotting

print(f"🎨 Preparing {len(results_df)} genes for volcano visualization...")

# SOPHISTICATED AESTHETIC MAPPING
# Use weighted random sampling to create visually distinct gene categories
# Genes with higher significance are more likely to be highlighted

print("🎭 Creating aesthetic gene categories...")

# First set of highlighted genes (weighted by significance)
picked_set1 = random.choices(
    results_df.symbol.tolist(),     # Gene symbols to choose from
    weights=results_df.nlog10.tolist(),  # Weight by significance level
    k=250                           # Number of genes to highlight
)

# Second set (exclusive of first set)
picked_set2 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=300)
picked_set2 = [x for x in picked_set2 if x not in picked_set1]  # Remove overlap

def assign_color(row):
    """
    Assign color category based on statistical significance and random selection.
    This creates visually appealing distributions while highlighting important genes.
    """
    fold_change, gene_symbol, minuslog10 = row
    
    # Genes with small changes or low significance get background color
    if abs(fold_change) < 1 or minuslog10 < 2:
        return "unremarkable"
    
    # Highlighted gene sets get distinct colors
    if gene_symbol in picked_set1:
        return "groupA"
    if gene_symbol in picked_set2:
        return "groupB"
    
    # Everything else gets intermediate color
    return "interesting"

# Apply color categories
results_df["color"] = results_df[["log2FoldChange", "symbol", "nlog10"]].apply(assign_color, axis=1)

# CREATE SHAPE CATEGORIES (for additional visual encoding)
picked_set3 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=250)
picked_set4 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=300)
picked_set4 = [x for x in picked_set4 if x not in picked_set3]

def assign_shape(gsym):
    """Assign marker shapes for additional visual encoding"""
    if gsym in picked_set3:
        return "shapeA"
    if gsym in picked_set4:
        return "shapeB"
    return "shapeNeutral"

results_df["shape"] = results_df.symbol.map(assign_shape)

# CREATE THE VOLCANO PLOT
print("🎨 Painting the volcano landscape...")

plt.figure(figsize=(10, 8))

# Main scatter plot with multiple aesthetic mappings
ax = sns.scatterplot(
    data=results_df,
    x="log2FoldChange",                                          # X: Effect size
    y="nlog10",                                                  # Y: Significance level
    hue="color",                                                 # Color by category
    hue_order=["unremarkable", "groupA", "groupB", "interesting"],
    palette=["lightgrey", "orange", "purple", "darkgray"],      # Custom color palette
    style="shape",                                               # Shape by category
    style_order=["shapeA", "shapeB", "shapeNeutral"],
    markers=["^", "s", "o"],                                     # Triangle, square, circle
    size="baseMean",                                             # Size by expression level
    sizes=(30, 200),                                             # Size range
    alpha=0.7                                                    # Semi-transparent
)

# Add significance threshold lines
ax.axhline(2, zorder=0, c="red", lw=2, ls="--", alpha=0.7, label="p-adj = 0.01")  # Significance line
ax.axvline(1, zorder=0, c="blue", lw=2, ls="--", alpha=0.7, label="2-fold up")     # Fold-change thresholds
ax.axvline(-1, zorder=0, c="blue", lw=2, ls="--", alpha=0.7, label="2-fold down")

# AUTOMATIC GENE LABELING
# Label the most dramatically changed genes for easy identification
print("🏷️ Labeling highly significant genes...")

label_texts = []
highly_significant_count = 0

for idx in range(len(results_df)):
    myrow = results_df.iloc[idx]
    
    # Label genes with extreme significance and large fold-changes
    if myrow.nlog10 > 5 and abs(myrow.log2FoldChange) > 2:
        label_texts.append(
            plt.text(
                x=myrow.log2FoldChange,
                y=myrow.nlog10,
                s=myrow.symbol,          # Gene symbol
                fontsize=11,
                weight="bold",
                color="black"
            )
        )
        highly_significant_count += 1

# Automatically adjust label positions to avoid overlap
if label_texts:
    adjust_text(
        label_texts, 
        arrowprops=dict(arrowstyle="-", color="black", alpha=0.6, lw=0.8)
    )

print(f"🏷️ Labeled {highly_significant_count} highly significant genes")

# PLOT STYLING AND FINISHING TOUCHES
# Position legend outside plot area
plt.legend(loc='center left', bbox_to_anchor=(1.05, 0.5), frameon=False, prop={"weight": "bold"})

# Clean up plot borders and styling
for border in ["bottom", "left"]:
    ax.spines[border].set_linewidth(2)
ax.spines["top"].set_visible(False)      # Remove top border
ax.spines["right"].set_visible(False)    # Remove right border
ax.tick_params(width=2)

# Enhance axis labels with mathematical formatting
plt.xticks(size=12, weight="bold")
plt.yticks(size=12, weight="bold")
plt.xlabel("$\\log_{2}$ Fold Change (mutant vs control)", size=14, weight="bold")
plt.ylabel("$-\\log_{10}$ Adjusted P-value", size=14, weight="bold")
plt.title("Volcano Plot: Statistical vs Biological Significance", size=16, weight="bold", pad=20)

# Add subtle grid for easier reading
plt.grid(True, alpha=0.2)

# Save high-quality figure
plt.savefig("volcano.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.show()

print("💎 Volcano plot masterpiece saved as 'volcano.png'!")

##########################################################
# 10. EXPRESSION HEATMAPS: Patterns in Gene Activity
##########################################################
"""
Create clustered heatmaps showing expression patterns of significant genes.
These reveal co-expression relationships and validate your experimental design.
"""

print("\n🔥 Creating expression pattern heatmaps...")

# First, we need to re-identify significant genes (column name changed above)
significant_hits = results_df[
    (results_df["padj"] < 0.1) & 
    (abs(results_df["log2FoldChange"]) > 0.5)
].copy()

if len(significant_hits) == 0:
    print("❌ No significant genes found for heatmap creation")
else:
    print(f"🎯 Creating heatmaps for {len(significant_hits)} significant genes...")
    
    # Identify top regulated genes for focused analysis
    top_up10 = significant_hits.sort_values("log2FoldChange", ascending=False).head(10)
    top_down10 = significant_hits.sort_values("log2FoldChange", ascending=True).head(10)
    highlighted_genes = list(top_up10.index) + list(top_down10.index)
    
    print(f"🔥 Top upregulated genes: {list(top_up10['symbol'].head(3))}")
    print(f"🧊 Top downregulated genes: {list(top_down10['symbol'].head(3))}")

    # Extract normalized counts for heatmap visualization
    try:
        norm_data = my_dds.layers["normed_counts"]
        
        # Create DataFrame with proper sample and gene labels
        norm_data_df = pd.DataFrame(
            data=norm_data,
            index=my_dds.obs_names,    # Sample names as rows
            columns=my_dds.var.index   # Gene names as columns
        )
        
        # HEATMAP A: All Significantly Changed Genes
        print("🎨 Creating heatmap of all significant genes...")
        
        if len(significant_hits) > 0:
            # Extract expression data for significant genes only
            sig_data_df = norm_data_df[significant_hits.index]
            
            # Create clustered heatmap with z-score normalization
            # Z-score normalization shows relative patterns across genes
            plt.figure(figsize=(10, 12))
            cluster_map = sns.clustermap(
                sig_data_df.T,           # Transpose: genes as rows, samples as columns
                z_score=0,               # Normalize each gene (row) to z-score
                cmap="viridis",          # Beautiful purple-yellow color scheme  
                figsize=(8, 10),
                cbar_kws={'label': 'Z-score\n(standard deviations from mean)'}
            )
            
            plt.suptitle(f"Expression Heatmap: All Significant Genes (n={len(significant_hits)})", 
                        y=0.98, size=14, weight="bold")
            plt.savefig("heatmap_all_sig_genes.png", dpi=300, bbox_inches="tight", facecolor="white")
            plt.show()
            print("💾 Complete heatmap saved as 'heatmap_all_sig_genes.png'")

        # HEATMAP B: Top Regulated Genes (Most Dramatic Changes)
        print("🎯 Creating focused heatmap of top regulated genes...")
        
        if len(highlighted_genes) > 0:
            # Extract expression data for top regulated genes
            top_data_df = norm_data_df[highlighted_genes]
            
            # Create focused heatmap with different color scheme
            plt.figure(figsize=(8, 10))
            cluster_map = sns.clustermap(
                top_data_df.T,           # Transpose for proper orientation
                z_score=0,               # Z-score normalization
                cmap="magma",            # Black-purple-yellow color scheme
                figsize=(8, 10),
                cbar_kws={'label': 'Z-score\n(standard deviations from mean)'}
            )
            
            plt.suptitle("Expression Heatmap: Top 10 Up & Top 10 Down Regulated Genes", 
                        y=0.98, size=14, weight="bold")
            plt.savefig("heatmap_top20_up_down.png", dpi=300, bbox_inches="tight", facecolor="white")
            plt.show()
            print("💾 Focused heatmap saved as 'heatmap_top20_up_down.png'")
            
    except AttributeError as e:
        print(f"⚠️ Could not create heatmaps: {e}")
        print("   This might be a PyDESeq2 version compatibility issue")

##########################################################
# 11. RESULTS EXPORT: Saving Your Discoveries
##########################################################
"""
Save all results in organized, analysis-ready formats for further investigation,
publication, or sharing with collaborators.
"""

print("\n💾 Saving analysis results...")

# Export complete differential expression results
results_df.to_csv("DE_results.csv")
print("📊 Complete DE results saved as 'DE_results.csv'")

# Create and save summary statistics
print("\n📋 FINAL ANALYSIS SUMMARY:")
print("=" * 60)
print(f"🧬 Total genes in analysis: {len(results_df)}")
print(f"⭐ Significantly changed genes: {len(significant_hits)}")

if len(significant_hits) > 0:
    upregulated = (significant_hits.log2FoldChange > 0).sum()
    downregulated = (significant_hits.log2FoldChange < 0).sum()
    print(f"📈 Upregulated genes: {upregulated}")
    print(f"📉 Downregulated genes: {downregulated}")
    
    # Show most extreme changes
    print(f"\n🔥 TOP 5 UPREGULATED GENES:")
    top_up = significant_hits.nlargest(5, 'log2FoldChange')
    for idx, row in top_up.iterrows():
        print(f"   {row['symbol']}: {row['log2FoldChange']:.2f} fold-change (padj={row['padj']:.2e})")
    
    print(f"\n🧊 TOP 5 DOWNREGULATED GENES:")
    top_down = significant_hits.nsmallest(5, 'log2FoldChange')
    for idx, row in top_down.iterrows():
        print(f"   {row['symbol']}: {row['log2FoldChange']:.2f} fold-change (padj={row['padj']:.2e})")

else:
    print("❌ No significantly changed genes found")
    print("💡 Consider relaxing thresholds or checking experimental design")

print("=" * 60)

# Save summary to file
with open("analysis_summary.txt", "w") as f:
    f.write("RNA-seq Differential Expression Analysis Summary\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    f.write(f"Total genes analyzed: {len(results_df)}\n")
    f.write(f"Significant genes found: {len(significant_hits)}\n")
    
    if len(significant_hits) > 0:
        upregulated = (significant_hits.log2FoldChange > 0).sum()
        downregulated = (significant_hits.log2FoldChange < 0).sum()
        f.write(f"Upregulated genes: {upregulated}\n")
        f.write(f"Downregulated genes: {downregulated}\n")
        
        f.write(f"\nSignificance Criteria Used:\n")
        f.write(f"- Adjusted p-value threshold: < 0.1\n")
        f.write(f"- Log2 fold-change threshold: > 0.5\n")
        f.write(f"- Minimum expression threshold: baseMean >= 10\n")

print("💾 Analysis summary saved as 'analysis_summary.txt'")

print("\n🎉 ANALYSIS COMPLETE!")
print("🎨 Generated visualizations:")
print("   📊 qc_metrics.png - Quality control assessment")
print("   🗺️ pca_plot.png - Sample relationship overview") 
print("   🔗 correlation_heatmap.png - Sample correlation matrix")
print("   📈 ma_plot.png - Expression vs fold-change relationship")
print("   🌋 volcano.png - Statistical significance landscape")
print("   🔥 heatmap_all_sig_genes.png - All significant gene patterns")
print("   🎯 heatmap_top20_up_down.png - Top regulated gene focus")
print("\n📊 Data files:")
print("   📋 DE_results.csv - Complete differential expression results")
print("   📝 analysis_summary.txt - Analysis summary and parameters")
print("\n✨ Your RNA-seq analysis journey is complete!")
print("🔬 Time to explore your biological discoveries!")

##########################################################
# END OF ANALYSIS PIPELINE
##########################################################

"""
CONGRATULATIONS!

You've successfully completed a comprehensive RNA-seq differential expression analysis!

WHAT YOU'VE ACCOMPLISHED:
✅ Assessed data quality and sample relationships
✅ Performed robust statistical testing with DESeq2
✅ Created publication-ready visualizations
✅ Identified significantly changed genes
✅ Generated organized results for further analysis

NEXT STEPS FOR YOUR RESEARCH:
🔍 Literature review: Research the functions of your top significant genes
🧪 Experimental validation: Design qPCR experiments to confirm key findings  
🌐 Pathway analysis: Use enrichment tools (DAVID, Enrichr) with your gene lists
📝 Manuscript preparation: Your plots are ready for publication!

REMEMBER:
- Your significant genes are biological hypotheses, not final answers
- Validation experiments strengthen your conclusions
- The patterns in your heatmaps might reveal unexpected gene relationships
- Each visualization tells a different part of your biological story

Happy researching!
"""
