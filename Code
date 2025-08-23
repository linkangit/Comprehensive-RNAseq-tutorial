#!/usr/bin/env python3
# coding: utf-8

"""
RNA-seq analysis script using PyDESeq2, Scanpy, and other libraries.

Usage:
    python rna_analysis_improved.py [--config config.yaml]
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scanpy as sc
import anndata
import argparse
import logging
import os
import sys
from pathlib import Path

# PyDESeq2
try:
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
except ImportError:
    print("Error: PyDESeq2 not installed. Install with: pip install pydeseq2")
    sys.exit(1)

# GSEA (optional)
try:
    import gseapy as gp
    from gseapy.plot import gseaplot
    GSEA_AVAILABLE = True
except ImportError:
    print("Warning: GSEA not available. Install with: pip install gseapy")
    GSEA_AVAILABLE = False

# For label adjustment
try:
    from adjustText import adjust_text
    ADJUST_TEXT_AVAILABLE = True
except ImportError:
    print("Warning: adjustText not available. Install with: pip install adjustText")
    ADJUST_TEXT_AVAILABLE = False

# Set random seed for reproducibility
np.random.seed(42)

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

#######################################################
# Configuration Parameters
#######################################################

class Config:
    """Configuration parameters for RNA-seq analysis"""
    
    # File paths
    READ_COUNTS_FILE = 'read_counts_table.csv'
    GENE_ANNOTATION_FILE = 'Arabidopsis_gene_annotation.tsv'
    OUTPUT_DIR = 'results'
    
    # Analysis parameters
    MIN_BASE_MEAN = 10
    PADJ_THRESHOLD = 0.1
    LOG2FC_THRESHOLD = 0.5
    VOLCANO_HIGHLIGHT_PADJ = 5  # -log10(padj) threshold for labeling
    VOLCANO_HIGHLIGHT_LOG2FC = 2  # log2FC threshold for labeling
    
    # Sample information
    SAMPLE_INFO = {
        "Sample": ["Kelley_17", "Kelley_18", "Kelley_19", "Kelley_20",
                   "Kelley_21", "Kelley_22", "Kelley_23", "Kelley_24"],
        "Condition": ["C", "C", "C", "C",
                      "mut", "mut", "mut", "mut"]
    }
    
    # Plot settings
    FIGURE_DPI = 300
    HEATMAP_FIGSIZE = (8, 10)
    VOLCANO_FIGSIZE = (8, 8)

def check_file_exists(filepath):
    """Check if file exists and raise error if not"""
    if not os.path.exists(filepath):
        raise FileNotFoundError(f"Required file not found: {filepath}")

def create_output_dir(output_dir):
    """Create output directory if it doesn't exist"""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

#######################################################
# 1. Prepare sample information (metadata)
#######################################################

def prepare_sample_info(config):
    """Prepare and validate sample information"""
    logger.info("Preparing sample information...")
    
    sample_info = pd.DataFrame(config.SAMPLE_INFO)
    logger.info(f"Sample info:\n{sample_info}")
    
    # Validate sample info
    if sample_info['Sample'].duplicated().any():
        raise ValueError("Duplicate sample names found in sample information")
    
    # Identify control vs. mutant groups
    ctrl_list = sample_info.loc[sample_info['Condition'] == 'C', 'Sample'].tolist()
    variant_list = sample_info.loc[sample_info['Condition'] == 'mut', 'Sample'].tolist()
    
    logger.info(f"Control samples: {ctrl_list}")
    logger.info(f"Mutant samples: {variant_list}")
    
    if len(ctrl_list) == 0 or len(variant_list) == 0:
        raise ValueError("Both control and mutant samples are required")
    
    return sample_info, ctrl_list, variant_list

#######################################################
# 2. Load and clean up read-count data
#######################################################

def load_and_clean_counts(config, all_samples):
    """Load and clean read count data"""
    logger.info("Loading read count data...")
    
    check_file_exists(config.READ_COUNTS_FILE)
    read_df = pd.read_csv(config.READ_COUNTS_FILE)
    
    logger.info(f"Initial read count table shape: {read_df.shape}")
    logger.info(f"Columns: {list(read_df.columns)}")
    
    # Check if 'gname' column exists
    if 'gname' not in read_df.columns:
        raise ValueError("'gname' column not found in read count data")
    
    # Make 'gname' column the index
    read_df = read_df.set_index('gname')
    
    # Check if all required samples are in the data
    missing_samples = [s for s in all_samples if s not in read_df.columns]
    if missing_samples:
        raise ValueError(f"Missing samples in read count data: {missing_samples}")
    
    # Filter columns based on sample list
    read_df = read_df[all_samples]
    
    # Exclude genes with zero total counts
    initial_genes = len(read_df)
    read_df = read_df[read_df.sum(axis=1) > 0]
    logger.info(f"Removed {initial_genes - len(read_df)} genes with zero counts")
    
    logger.info(f"Final read count table shape: {read_df.shape}")
    
    return read_df

#######################################################
# 3. Quality control metrics
#######################################################

def calculate_qc_metrics(read_df, sample_info):
    """Calculate and display quality control metrics"""
    logger.info("Calculating quality control metrics...")
    
    qc_metrics = pd.DataFrame({
        'Total_Reads': read_df.sum(axis=0),
        'Detected_Genes': (read_df > 0).sum(axis=0),
        'Condition': sample_info.set_index('Sample')['Condition']
    })
    
    logger.info("Quality control metrics:")
    logger.info(f"\n{qc_metrics}")
    
    # Plot QC metrics
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Total reads per sample
    sns.barplot(data=qc_metrics.reset_index(), x='index', y='Total_Reads', 
                hue='Condition', ax=axes[0])
    axes[0].set_title('Total Reads per Sample')
    axes[0].tick_params(axis='x', rotation=45)
    
    # Detected genes per sample
    sns.barplot(data=qc_metrics.reset_index(), x='index', y='Detected_Genes', 
                hue='Condition', ax=axes[1])
    axes[1].set_title('Detected Genes per Sample')
    axes[1].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{Config.OUTPUT_DIR}/qc_metrics.png", dpi=Config.FIGURE_DPI, bbox_inches="tight")
    plt.show()
    
    return qc_metrics

#######################################################
# 4. Prepare data for DESeq2 and run analysis
#######################################################

def run_deseq2_analysis(read_df, sample_info):
    """Run DESeq2 differential expression analysis"""
    logger.info("Running DESeq2 analysis...")
    
    # Prepare sample info for DESeq2
    sample_info_deseq = sample_info.set_index('Sample')
    sample_info_deseq["Condition"] = sample_info_deseq["Condition"].astype("category")
    sample_info_deseq["Condition"] = sample_info_deseq["Condition"].cat.reorder_categories(["C", "mut"], ordered=True)
    
    # PyDESeq2 expects samples in rows and genes in columns
    transposed_df = read_df.T
    
    try:
        my_dds = DeseqDataSet(
            counts=transposed_df,
            metadata=sample_info_deseq,
            design_factors=["Condition"]
        )
        
        logger.info(f"DESeq2 dataset created: {my_dds}")
        
        my_dds.deseq2()
        logger.info("DESeq2 analysis completed")
        
        return my_dds
        
    except Exception as e:
        logger.error(f"Error in DESeq2 analysis: {e}")
        raise

#######################################################
# 5. Extract and process DE results
#######################################################

def extract_de_results(my_dds, config):
    """Extract and process differential expression results"""
    logger.info("Extracting DE results...")
    
    try:
        my_stats = DeseqStats(my_dds, contrast=["Condition", "mut", "C"])
        stats_summary = my_stats.summary()
        results_df = my_stats.results_df.copy()
        
        logger.info(f"DE analysis summary:\n{stats_summary}")
        
        return results_df
        
    except Exception as e:
        logger.error(f"Error extracting DE results: {e}")
        raise

#######################################################
# 6. Add gene annotations
#######################################################

def add_gene_annotations(results_df, config):
    """Add gene annotations to results"""
    logger.info("Adding gene annotations...")
    
    if not os.path.exists(config.GENE_ANNOTATION_FILE):
        logger.warning(f"Gene annotation file not found: {config.GENE_ANNOTATION_FILE}")
        results_df['Symbol'] = results_df.index
        return results_df
    
    try:
        gene_info = pd.read_csv(config.GENE_ANNOTATION_FILE, sep='\t')
        
        if 'Nomenclature ID' not in gene_info.columns:
            logger.warning("'Nomenclature ID' column not found in annotation file")
            results_df['Symbol'] = results_df.index
            return results_df
            
        gene_info.set_index('Nomenclature ID', inplace=True)
        
        # Filter out rows without symbols (fixed from original)
        if 'Symbol' in gene_info.columns:
            gene_info = gene_info.dropna(subset=['Symbol'])
        
        # Merge annotations
        complete_df = pd.merge(
            left=results_df,
            right=gene_info[['Symbol']] if 'Symbol' in gene_info.columns else pd.DataFrame(index=gene_info.index),
            how="left",
            left_index=True,
            right_index=True
        )
        
        # Use gene ID as symbol if no symbol available
        if 'Symbol' in complete_df.columns:
            complete_df['Symbol'] = complete_df['Symbol'].fillna(complete_df.index)
        else:
            complete_df['Symbol'] = complete_df.index
            
        results_df = complete_df
        logger.info(f"Added annotations for {(~results_df['Symbol'].isnull()).sum()} genes")
        
        return results_df
        
    except Exception as e:
        logger.warning(f"Error loading gene annotations: {e}")
        results_df['Symbol'] = results_df.index
        return results_df

#######################################################
# 7. Filter results and identify significant genes
#######################################################

def filter_and_identify_significant(results_df, config):
    """Filter results and identify significant genes"""
    logger.info("Filtering results and identifying significant genes...")
    
    # Filter out low-expression genes
    initial_count = len(results_df)
    results_df = results_df[results_df.baseMean >= config.MIN_BASE_MEAN]
    logger.info(f"Filtered out {initial_count - len(results_df)} low-expression genes")
    
    # Remove genes with missing values
    results_df = results_df.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    
    # Find significant genes
    significant_mask = (results_df.padj < config.PADJ_THRESHOLD) & \
                      (abs(results_df.log2FoldChange) > config.LOG2FC_THRESHOLD)
    significant_hits = results_df[significant_mask].copy()
    
    logger.info(f"Identified {len(significant_hits)} significant genes")
    logger.info(f"Upregulated: {(significant_hits.log2FoldChange > 0).sum()}")
    logger.info(f"Downregulated: {(significant_hits.log2FoldChange < 0).sum()}")
    
    return results_df, significant_hits

#######################################################
# 8. Visualization functions
#######################################################

def create_pca_plot(my_dds, config):
    """Create PCA plot"""
    logger.info("Creating PCA plot...")
    
    try:
        # Try different ways to access normalized counts
        norm_data = None
        if hasattr(my_dds, 'layers') and 'normed_counts' in my_dds.layers:
            norm_data = my_dds.layers["normed_counts"]
        elif hasattr(my_dds, 'X'):
            norm_data = my_dds.X
        else:
            logger.warning("Could not find normalized counts for PCA")
            return None
        
        if norm_data is not None:
            # Build AnnData object
            adata = anndata.AnnData(
                X=norm_data, 
                obs=my_dds.obs, 
                var=pd.DataFrame(index=my_dds.var.index)
            )
            
            # Log transform and scale
            adata.layers["log1p"] = np.log1p(norm_data)
            sc.pp.scale(adata, max_value=10)
            sc.tl.pca(adata)
            
            # Plot PCA
            sc.pl.pca(adata, color='Condition', size=150, show=False)
            plt.title("PCA of Samples")
            plt.savefig(f"{config.OUTPUT_DIR}/pca_plot.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
            plt.show()
            
            return adata
            
    except Exception as e:
        logger.warning(f"Error creating PCA plot: {e}")
        return None

def create_correlation_heatmap(my_dds, config):
    """Create sample correlation heatmap"""
    logger.info("Creating correlation heatmap...")
    
    try:
        # Get normalized counts
        norm_data = None
        if hasattr(my_dds, 'layers') and 'normed_counts' in my_dds.layers:
            norm_data = my_dds.layers["normed_counts"]
        elif hasattr(my_dds, 'X'):
            norm_data = my_dds.X
            
        if norm_data is not None:
            # Calculate correlation matrix (samples x samples)
            corr_mat = np.corrcoef(norm_data)
            
            # Create heatmap
            plt.figure(figsize=(8, 6))
            sns.heatmap(
                corr_mat,
                annot=True,
                cmap="viridis",
                xticklabels=my_dds.obs.index,
                yticklabels=my_dds.obs.index,
                square=True
            )
            plt.title("Sample Correlation Heatmap")
            plt.tight_layout()
            plt.savefig(f"{config.OUTPUT_DIR}/correlation_heatmap.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
            plt.show()
            
    except Exception as e:
        logger.warning(f"Error creating correlation heatmap: {e}")

def create_ma_plot(results_df, config):
    """Create MA plot"""
    logger.info("Creating MA plot...")
    
    plt.figure(figsize=(8, 6))
    
    # Color points by significance
    significant_mask = (results_df.padj < config.PADJ_THRESHOLD) & \
                      (abs(results_df.log2FoldChange) > config.LOG2FC_THRESHOLD)
    
    # Plot non-significant points
    plt.scatter(
        x=results_df.loc[~significant_mask, 'baseMean'],
        y=results_df.loc[~significant_mask, 'log2FoldChange'],
        alpha=0.3,
        c='lightgray',
        s=20,
        label='Non-significant'
    )
    
    # Plot significant points
    plt.scatter(
        x=results_df.loc[significant_mask, 'baseMean'],
        y=results_df.loc[significant_mask, 'log2FoldChange'],
        alpha=0.7,
        c='red',
        s=30,
        label='Significant'
    )
    
    plt.xscale('log')
    plt.axhline(y=0, color='black', linestyle='--', alpha=0.7)
    plt.axhline(y=config.LOG2FC_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
    plt.axhline(y=-config.LOG2FC_THRESHOLD, color='blue', linestyle='--', alpha=0.5)
    
    plt.xlabel("Mean Normalized Counts (log scale)")
    plt.ylabel("Log2 Fold Change")
    plt.title("MA Plot")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(f"{config.OUTPUT_DIR}/ma_plot.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
    plt.show()

def create_volcano_plot(results_df, config):
    """Create volcano plot"""
    logger.info("Creating volcano plot...")
    
    # Prepare data
    plot_data = results_df.copy()
    plot_data["nlog10_padj"] = -np.log10(plot_data["padj"])
    
    # Define significance categories
    def get_significance_category(row):
        padj_thresh = -np.log10(config.PADJ_THRESHOLD)
        if (abs(row['log2FoldChange']) >= config.LOG2FC_THRESHOLD and 
            row['nlog10_padj'] >= padj_thresh):
            if row['log2FoldChange'] > 0:
                return 'Upregulated'
            else:
                return 'Downregulated'
        return 'Not significant'
    
    plot_data['category'] = plot_data.apply(get_significance_category, axis=1)
    
    # Create plot
    plt.figure(figsize=config.VOLCANO_FIGSIZE)
    
    # Define colors
    colors = {'Not significant': 'lightgray', 'Upregulated': 'red', 'Downregulated': 'blue'}
    
    for category in colors.keys():
        subset = plot_data[plot_data['category'] == category]
        plt.scatter(
            subset['log2FoldChange'],
            subset['nlog10_padj'],
            c=colors[category],
            alpha=0.6,
            s=20,
            label=f"{category} ({len(subset)})"
        )
    
    # Add threshold lines
    plt.axhline(-np.log10(config.PADJ_THRESHOLD), color='black', linestyle='--', alpha=0.7)
    plt.axvline(config.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)
    plt.axvline(-config.LOG2FC_THRESHOLD, color='black', linestyle='--', alpha=0.7)
    
    # Add labels for highly significant genes
    if ADJUST_TEXT_AVAILABLE:
        label_texts = []
        highly_significant = plot_data[
            (plot_data['nlog10_padj'] > config.VOLCANO_HIGHLIGHT_PADJ) &
            (abs(plot_data['log2FoldChange']) > config.VOLCANO_HIGHLIGHT_LOG2FC)
        ]
        
        for idx, row in highly_significant.iterrows():
            label_texts.append(
                plt.text(
                    row['log2FoldChange'],
                    row['nlog10_padj'],
                    row['Symbol'],
                    fontsize=10,
                    weight="bold"
                )
            )
        
        if label_texts:
            adjust_text(label_texts, arrowprops=dict(arrowstyle="-", color="black", alpha=0.7))
    
    plt.xlabel("Log2 Fold Change")
    plt.ylabel("-Log10 Adjusted P-value")
    plt.title("Volcano Plot")
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.savefig(f"{config.OUTPUT_DIR}/volcano_plot.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
    plt.show()

def create_expression_heatmaps(my_dds, significant_hits, config):
    """Create expression heatmaps"""
    logger.info("Creating expression heatmaps...")
    
    try:
        # Get normalized counts
        norm_data = None
        if hasattr(my_dds, 'layers') and 'normed_counts' in my_dds.layers:
            norm_data = my_dds.layers["normed_counts"]
        elif hasattr(my_dds, 'X'):
            norm_data = my_dds.X
            
        if norm_data is None:
            logger.warning("Could not find normalized counts for heatmaps")
            return
            
        norm_data_df = pd.DataFrame(
            data=norm_data,
            index=my_dds.obs_names,
            columns=my_dds.var.index
        )
        
        # Heatmap 1: All significant genes
        if len(significant_hits) > 0:
            sig_data = norm_data_df[significant_hits.index]
            
            plt.figure(figsize=config.HEATMAP_FIGSIZE)
            sns.clustermap(
                sig_data.T, 
                z_score=0, 
                cmap="RdBu_r", 
                figsize=config.HEATMAP_FIGSIZE,
                cbar_kws={'label': 'Z-score'}
            )
            plt.suptitle(f"Heatmap: All Significant Genes (n={len(significant_hits)})", y=0.98)
            plt.savefig(f"{config.OUTPUT_DIR}/heatmap_all_significant.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
            plt.show()
            
            # Heatmap 2: Top regulated genes
            if len(significant_hits) >= 20:
                top_up = significant_hits.nlargest(10, 'log2FoldChange')
                top_down = significant_hits.nsmallest(10, 'log2FoldChange')
                top_genes = list(top_up.index) + list(top_down.index)
                
                top_data = norm_data_df[top_genes]
                
                plt.figure(figsize=config.HEATMAP_FIGSIZE)
                sns.clustermap(
                    top_data.T, 
                    z_score=0, 
                    cmap="RdBu_r", 
                    figsize=config.HEATMAP_FIGSIZE,
                    cbar_kws={'label': 'Z-score'}
                )
                plt.suptitle("Heatmap: Top 10 Up and Down Regulated Genes", y=0.98)
                plt.savefig(f"{config.OUTPUT_DIR}/heatmap_top20.png", dpi=config.FIGURE_DPI, bbox_inches="tight")
                plt.show()
        else:
            logger.info("No significant genes found for heatmap")
            
    except Exception as e:
        logger.warning(f"Error creating heatmaps: {e}")

#######################################################
# 9. GSEA Analysis (optional)
#######################################################

def run_gsea_analysis(significant_hits, config):
    """Run GSEA analysis if available"""
    if not GSEA_AVAILABLE:
        logger.info("GSEA not available, skipping pathway analysis")
        return
    
    logger.info("Running GSEA analysis...")
    
    try:
        # Prepare gene list for GSEA
        gene_list = significant_hits[['Symbol', 'log2FoldChange']].copy()
        gene_list = gene_list.sort_values('log2FoldChange', ascending=False)
        
        # Run enrichment analysis
        enr = gp.enrichr(
            gene_list=gene_list['Symbol'].tolist(),
            gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
            organism='plant',
            outdir=f"{config.OUTPUT_DIR}/gsea"
        )
        
        logger.info("GSEA analysis completed")
        
        # Save results
        if hasattr(enr, 'results'):
            enr.results.to_csv(f"{config.OUTPUT_DIR}/gsea_results.csv", index=False)
        
    except Exception as e:
        logger.warning(f"Error in GSEA analysis: {e}")

#######################################################
# 10. Save results
#######################################################

def save_results(results_df, significant_hits, config):
    """Save analysis results"""
    logger.info("Saving results...")
    
    # Save all results
    results_df.to_csv(f"{config.OUTPUT_DIR}/DE_results_all.csv")
    
    # Save significant results
    if len(significant_hits) > 0:
        significant_hits.to_csv(f"{config.OUTPUT_DIR}/DE_results_significant.csv")
    
    # Create summary report
    summary = {
        'Total genes analyzed': len(results_df),
        'Significant genes': len(significant_hits),
        'Upregulated genes': (significant_hits.log2FoldChange > 0).sum() if len(significant_hits) > 0 else 0,
        'Downregulated genes': (significant_hits.log2FoldChange < 0).sum() if len(significant_hits) > 0 else 0,
        'Parameters used': {
            'Min base mean': config.MIN_BASE_MEAN,
            'Padj threshold': config.PADJ_THRESHOLD,
            'Log2FC threshold': config.LOG2FC_THRESHOLD
        }
    }
    
    # Save summary
    with open(f"{config.OUTPUT_DIR}/analysis_summary.txt", 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    
    logger.info("Results saved successfully")

#######################################################
# Main analysis function
#######################################################

def main():
    """Main analysis pipeline"""
    # Initialize configuration
    config = Config()
    
    # Create output directory
    create_output_dir(config.OUTPUT_DIR)
    
    try:
        # 1. Prepare sample information
        sample_info, ctrl_list, variant_list = prepare_sample_info(config)
        all_samples = ctrl_list + variant_list
        
        # 2. Load and clean count data
        read_df = load_and_clean_counts(config, all_samples)
        
        # 3. Calculate QC metrics
        qc_metrics = calculate_qc_metrics(read_df, sample_info)
        
        # 4. Run DESeq2 analysis
        my_dds = run_deseq2_analysis(read_df, sample_info)
        
        # 5. Extract DE results
        results_df = extract_de_results(my_dds, config)
        
        # 6. Add gene annotations
        results_df = add_gene_annotations(results_df, config)
        
        # 7. Filter and identify significant genes
        results_df, significant_hits = filter_and_identify_significant(results_df, config)
        
        # 8. Create visualizations
        adata = create_pca_plot(my_dds, config)
        create_correlation_heatmap(my_dds, config)
        create_ma_plot(results_df, config)
        create_volcano_plot(results_df, config)
        create_expression_heatmaps(my_dds, significant_hits, config)
        
        # 9. Run GSEA analysis
        if len(significant_hits) > 0:
            run_gsea_analysis(significant_hits, config)
        
        # 10. Save results
        save_results(results_df, significant_hits, config)
        
        logger.info("Analysis completed successfully!")
        logger.info(f"Results saved to: {config.OUTPUT_DIR}")
        
    except Exception as e:
        logger.error(f"Analysis failed: {e}")
        raise

if __name__ == "__main__":
    main()
