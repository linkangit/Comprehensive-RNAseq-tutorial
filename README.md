# RNA-seq Analysis Tutorial

This tutorial walks through a complete RNA-seq differential expression analysis using Python with PyDESeq2, Scanpy, and visualization libraries. The pipeline is designed for research workflows with a focus on simplicity and effectiveness.

## Table of Contents
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Input Data Requirements](#input-data-requirements)
4. [Installation](#installation)
5. [Code Walkthrough](#code-walkthrough)
6. [Understanding the Outputs](#understanding-the-outputs)
7. [Customizing for Your Data](#customizing-for-your-data)
8. [Troubleshooting](#troubleshooting)

## Overview

This analysis pipeline performs differential expression analysis comparing control vs. mutant samples using a streamlined approach. The workflow includes:

- **Sample metadata preparation** and validation
- **Count data loading** and filtering
- **DESeq2 differential expression** analysis
- **Comprehensive visualizations**: PCA, correlation heatmap, MA plot, and volcano plot
- **Expression heatmaps** for significant genes
- **Results export** in CSV format

### Key Features
- Clean, linear workflow easy to understand and modify
- Rich visualizations with multiple aesthetic mappings
- Minimal dependencies and straightforward execution
- Research-focused approach without over-engineering

## Prerequisites

### Required Knowledge
- Basic understanding of RNA-seq data analysis
- Familiarity with Python and pandas
- Understanding of differential expression concepts

### Required Software
- Python 3.7 or higher
- Standard data science libraries

## Input Data Requirements

You need **exactly 2 files** in your working directory:

### 1. `read_counts_table.csv` (REQUIRED)
A CSV file containing raw read counts:
- **First column:** `gname` - Gene identifiers (e.g., AT1G01010)
- **Subsequent columns:** Sample names (must match sample metadata in script)
- **Values:** Integer read counts

```csv
gname,Kelley_17,Kelley_18,Kelley_19,Kelley_20,Kelley_21,Kelley_22,Kelley_23,Kelley_24
AT1G01010,523,445,678,512,1205,1456,1123,1087
AT1G01020,89,92,76,85,45,38,52,61
AT1G01030,1234,1456,1123,1345,234,198,267,312
```

### 2. `Arabidopsis_gene_annotation.tsv` (REQUIRED)
A tab-separated file with gene annotations:
- **Nomenclature ID:** Gene identifiers matching the count matrix
- **Symbol:** Gene symbols/names
- **Additional columns:** Optional annotation information

```tsv
Nomenclature ID	Symbol	Description
AT1G01010	NAC001	NAC domain containing protein 1
AT1G01020	ARV1	ARV1 family protein
AT1G01030	NGA3	NGATHA3
```

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/linkangit/Comprehensive-RNAseq-tutorial.git
cd Comprehensive-RNAseq-tutorial
```

### 2. Install Required Packages
```bash
pip install pandas numpy seaborn matplotlib scanpy anndata pydeseq2 gseapy adjustText
```

### 3. Verify Installation
```python
python -c "import pydeseq2, scanpy, pandas; print('Ready to analyze!')"
```

## Code Walkthrough

Let's walk through each section of the analysis script:

### Section 1: Sample Information Setup

```python
sample_info_dict = {
    "Sample": ["Kelley_17", "Kelley_18", "Kelley_19", "Kelley_20",
               "Kelley_21", "Kelley_22", "Kelley_23", "Kelley_24"],
    "Condition": ["C", "C", "C", "C",
                  "mut", "mut", "mut", "mut"]
}

sample_info = pd.DataFrame(sample_info_dict)
```

**What happens here:**
- Creates a dictionary linking each sample to its experimental condition
- `"C"` represents control samples, `"mut"` represents mutant samples
- Converts to DataFrame for easier manipulation
- Identifies control and mutant sample groups

**Key Point:** You'll need to modify this section to match your sample names and conditions.

### Section 2: Count Data Loading

```python
read_df = pd.read_csv('read_counts_table.csv')
read_df = read_df.set_index('gname')
read_df = read_df[read_df.sum(axis=1) > 0]
```

**What happens here:**
- Loads the count matrix from CSV
- Sets gene names (`gname`) as row indices
- Removes genes with zero counts across all samples (these can't be analyzed statistically)

**Result:** Clean count matrix with genes as rows, samples as columns

### Section 3: Data Filtering and Preparation

```python
all_samples = ctrl_list + variant_list
read_df = read_df[all_samples]

sample_info.set_index('Sample', inplace=True)
sample_info["Condition"] = sample_info["Condition"].astype("category")
sample_info["Condition"] = sample_info["Condition"].cat.reorder_categories(["C", "mut"], ordered=True)
```

**What happens here:**
- Filters count matrix to include only the samples in your experiment
- Prepares sample metadata for DESeq2 by setting proper index
- Converts conditions to categorical data with proper ordering (control first, then mutant)

### Section 3.5: Quality Control Visualization

```python
# Calculate QC metrics
qc_metrics = pd.DataFrame({
    'Total_Reads': read_df.sum(axis=0),
    'Detected_Genes': (read_df > 0).sum(axis=0),
    'Condition': sample_info['Condition']
})

# Create QC plots
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
plt.savefig("qc_metrics.png", dpi=300, bbox_inches="tight", facecolor="white")
plt.show()
```

**What happens here:**
- Calculates total read counts per sample to check for library size differences
- Counts detected genes (reads > 0) per sample to assess library complexity
- Creates side-by-side bar plots comparing control vs. mutant samples
- Helps identify potential outlier samples or batch effects early
- Saves high-quality plot for inclusion in reports

**Why this matters:** These QC metrics help you spot issues before running expensive statistical analysis. Samples with very different total reads or detected genes might indicate technical problems.

### Section 4: DESeq2 Analysis

```python
transposed_df = read_df.T

my_dds = DeseqDataSet(
    counts=transposed_df,
    metadata=sample_info,
    design_factors=["Condition"]
)

my_dds.deseq2()
```

**What happens here:**
- Transposes the count matrix (DESeq2 expects samples as rows, genes as columns)
- Creates a DESeqDataSet object with your count data, sample metadata, and experimental design
- Runs the DESeq2 algorithm which:
  - Estimates size factors (normalizes for library size differences)
  - Estimates gene-wise dispersions (accounts for count variability)
  - Performs statistical testing using negative binomial distribution

### Section 5: Extract Results

```python
my_stats = DeseqStats(my_dds, contrast=["Condition", "mut", "C"])
results_df = my_stats.results_df
```

**What happens here:**
- Defines the statistical contrast: mutant vs. control
- Extracts differential expression results including:
  - `log2FoldChange`: Effect size (positive = upregulated in mutant)
  - `padj`: Adjusted p-values (corrected for multiple testing)
  - `baseMean`: Average normalized expression across samples

### Section 6: Gene Annotation Integration

```python
gene_info = pd.read_csv('Arabidopsis_gene_annotation.tsv', sep='\t')
gene_info.set_index('Nomenclature ID', inplace=True)

complete_df = pd.merge(
    left=results_df,
    right=gene_info,
    how="inner",
    left_index=True,
    right_index=True
)

results_df['Symbol'] = complete_df['Symbol']
```

**What happens here:**
- Loads gene annotation file
- Merges statistical results with gene symbols and descriptions
- Adds readable gene names to the results for easier interpretation

### Section 7: Significance Filtering

```python
results_df = results_df[results_df.baseMean >= 10]
significant_hits = results_df[(results_df.padj < 0.1) & (abs(results_df.log2FoldChange) > 0.5)]
```

**What happens here:**
- Filters out lowly expressed genes (baseMean < 10)
- Identifies statistically and biologically significant genes:
  - Adjusted p-value < 0.1 (10% false discovery rate)
  - Absolute log2 fold-change > 0.5 (at least 1.4-fold change)

### Section 8: Exploratory Visualizations

#### PCA Analysis
```python
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata)
sc.pl.pca(adata, color='Condition', size=150, show=False)
```

**Purpose:** Shows overall sample relationships and clustering by condition

#### Correlation Heatmap
```python
corr_mat = np.corrcoef(norm_data)
sns.clustermap(corr_mat, cmap="vlag")
```

**Purpose:** Visualizes sample-to-sample correlations and potential batch effects

#### MA Plot
```python
plt.scatter(x=results_df['baseMean'], y=results_df['log2FoldChange'])
plt.xscale('log')
```

**Purpose:** Shows relationship between gene expression level and fold-change magnitude

### Section 9: Advanced Volcano Plot

The volcano plot section is particularly sophisticated:

```python
# Create significance categories with weighted random sampling
picked_set1 = random.choices(results_df.symbol.tolist(), weights=results_df.nlog10.tolist(), k=250)

def assign_color(row):
    fold_change, gene_symbol, minuslog10 = row
    if abs(fold_change) < 1 or minuslog10 < 2:
        return "unremarkable"
    if gene_symbol in picked_set1:
        return "groupA"
    return "interesting"
```

**What's clever here:**
- Uses weighted random sampling based on significance levels
- Creates visually distinct categories for different gene groups
- Combines multiple aesthetic mappings (color, shape, size) for rich visualization
- Automatically labels highly significant genes

### Section 10: Expression Heatmaps

```python
# All significant genes
sns.clustermap(sign_data_df.T, z_score=0, cmap="viridis")

# Top regulated genes  
top_up10 = significant_hits.sort_values("log2FoldChange", ascending=False).head(10)
top_down10 = significant_hits.sort_values("log2FoldChange", ascending=True).head(10)
```

**What happens here:**
- Creates clustered heatmaps showing expression patterns
- Z-score normalization makes patterns visible across different expression levels
- Separate heatmap for most highly regulated genes (top 10 up + top 10 down)

## Understanding the Outputs

After running the script, you'll get:

### Generated Files
- **`qc_metrics.png`** - Quality control plots showing read counts and detected genes per sample
- **`volcano.png`** - Volcano plot showing statistical vs. biological significance
- **`heatmap_all_sig_genes.png`** - Clustered heatmap of all significant genes
- **`heatmap_top20_up_down.png`** - Heatmap of most highly regulated genes
- **`DE_results.csv`** - Complete differential expression results

### Interactive Plots (Displayed During Execution)
- **QC metrics** - Bar plots comparing read counts and gene detection between conditions
- **PCA plot** - Sample clustering and outlier detection
- **Correlation heatmap** - Sample-to-sample relationships
- **MA plot** - Expression level vs. fold-change relationship

### Key Result Columns

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized expression across all samples |
| `log2FoldChange` | Log2 fold-change (mutant vs. control) |
| `pvalue` | Raw p-value from statistical test |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `Symbol` | Gene symbol from annotation file |

### Interpreting Results

**Significant Genes:**
- `padj < 0.1` AND `|log2FoldChange| > 0.5`
- Positive log2FoldChange = upregulated in mutant
- Negative log2FoldChange = downregulated in mutant

**Volcano Plot:**
- X-axis: Log2 fold-change (effect size)
- Y-axis: -log10(adjusted p-value) (significance)
- Different colors/shapes represent gene categories
- Labeled genes: highly significant (padj < 1e-5, |log2FC| > 2)

## Customizing for Your Data

### 1. Update Sample Information
Edit the `sample_info_dict` to match your samples:

```python
sample_info_dict = {
    "Sample": ["Sample_1", "Sample_2", "Sample_3", "Sample_4",
               "Sample_5", "Sample_6", "Sample_7", "Sample_8"],
    "Condition": ["control", "control", "control", "control",
                  "treatment", "treatment", "treatment", "treatment"]
}
```

**Important:** 
- Sample names must exactly match column headers in your count matrix
- Use consistent condition labels (e.g., "control"/"treatment" or "WT"/"mutant")

### 2. Adjust File Names (if needed)
If your files have different names, update these lines:

```python
read_df = pd.read_csv('your_count_file.csv')  # Line ~45
gene_info = pd.read_csv('your_annotation_file.tsv', sep='\t')  # Line ~85
```

### 3. Modify Significance Thresholds (optional)
```python
# Current thresholds
results_df = results_df[results_df.baseMean >= 10]  # Minimum expression
significant_hits = results_df[(results_df.padj < 0.1) & (abs(results_df.log2FoldChange) > 0.5)]

# Example: More stringent
significant_hits = results_df[(results_df.padj < 0.05) & (abs(results_df.log2FoldChange) > 1.0)]
```

## Folder Structure

Your working directory should look like this:

```
your_project_folder/
├── code.py                            # The analysis script
├── read_counts_table.csv              # Your count matrix
├── Arabidopsis_gene_annotation.tsv    # Gene annotations
└── (results will be generated here)
```

## Running the Analysis

```bash
# Navigate to your project folder
cd your_project_folder

# Ensure your data files are present
ls  # Should show: code.py, read_counts_table.csv, Arabidopsis_gene_annotation.tsv

# Run the analysis
python code.py
```

## Detailed Code Explanation

### The Beauty of This Approach

This script follows a clean, research-oriented workflow that prioritizes:

1. **Immediate feedback** - Print statements show progress and intermediate results
2. **Visual exploration** - Multiple plot types reveal different aspects of the data
3. **Flexible filtering** - Easy to adjust thresholds and criteria
4. **Rich visualizations** - The volcano plot uses sophisticated aesthetic mappings

### Key Sections Explained

#### Sample Metadata (Lines 25-35)
```python
sample_info_dict = {...}
sample_info = pd.DataFrame(sample_info_dict)
ctrl_list = sample_info.loc[sample_info['Condition'] == 'C', 'Sample'].tolist()
variant_list = sample_info.loc[sample_info['Condition'] == 'mut', 'Sample'].tolist()
```
Creates structured sample information and automatically identifies control/mutant groups.

#### Count Data Processing (Lines 40-50)
```python
read_df = pd.read_csv('read_counts_table.csv')
read_df = read_df.set_index('gname')
read_df = read_df[read_df.sum(axis=1) > 0]
```
Loads, indexes, and cleans count data by removing unexpressed genes.

#### DESeq2 Setup (Lines 65-75)
```python
transposed_df = read_df.T
my_dds = DeseqDataSet(counts=transposed_df, metadata=sample_info, design_factors=["Condition"])
my_dds.deseq2()
```
The core statistical analysis using the DESeq2 algorithm for robust differential expression testing.

#### Volcano Plot Magic (Lines 120-180)
The volcano plot section is particularly elegant:
- **Weighted random sampling** creates visually distinct gene categories
- **Multiple aesthetics** (color, shape, size) encode different information
- **Automatic labeling** highlights the most significant genes
- **Professional styling** with custom axis formatting

## Expected Results

### Typical Output
- **Console output:** Progress messages and summary statistics
- **Interactive plots:** PCA, correlation heatmap, MA plot (displayed during execution)
- **Saved files:** Volcano plot and heatmaps as PNG files
- **Data export:** Complete results in CSV format

### What Good Results Look Like

1. **PCA Plot:** Samples should cluster by condition, not by batch
2. **Correlation Heatmap:** Within-condition samples should be more similar
3. **MA Plot:** Significant genes should be distributed across expression levels
4. **Volcano Plot:** Clear separation of significant genes from background
5. **Heatmaps:** Distinct expression patterns between conditions

### Troubleshooting Common Issues

#### No Significant Genes Found
```python
# Check your thresholds
print(f"Genes with padj < 0.2: {(results_df.padj < 0.2).sum()}")
print(f"Genes with |log2FC| > 0.25: {(abs(results_df.log2FoldChange) > 0.25).sum()}")

# Relax criteria if needed
significant_hits = results_df[(results_df.padj < 0.2) & (abs(results_df.log2FoldChange) > 0.25)]
```

#### File Loading Errors
- Verify file names match exactly
- Check that sample names in CSV match the script
- Ensure annotation file uses tab separation

#### Memory Issues with Large Datasets
```python
# Pre-filter low-count genes more aggressively
read_df = read_df[read_df.sum(axis=1) >= 50]  # Increase from default filtering
```

## Customization Examples

### Different Organism
For non-Arabidopsis data, update the annotation file and GSEA organism:

```python
# Change annotation file name
gene_info = pd.read_csv('your_organism_annotation.tsv', sep='\t')

# Update GSEA organism (if using)
# organism='human', 'mouse', 'rat', etc.
```

### Different Experimental Design
For time-course or multi-factor experiments:

```python
sample_info_dict = {
    "Sample": ["Sample_1", "Sample_2", "Sample_3", "Sample_4"],
    "Treatment": ["Control", "Control", "Drug", "Drug"],
    "Timepoint": ["6h", "24h", "6h", "24h"]
}

# Update DESeq2 design
my_dds = DeseqDataSet(
    counts=transposed_df,
    metadata=sample_info,
    design_factors=["Treatment", "Timepoint"]  # Multi-factor design
)
```

### Adjust Visualization Parameters
```python
# Volcano plot labeling threshold
if myrow.nlog10 > 3 and abs(myrow.log2FoldChange) > 1.5:  # More stringent

# Heatmap parameters
sns.clustermap(sign_data_df.T, z_score=0, cmap="RdBu_r", figsize=(12, 8))  # Different colors/size
```

## Advanced Tips

### 1. Batch Effect Detection
If samples don't cluster by condition in PCA:
- Check for batch effects in sample collection/processing
- Consider adding batch as a covariate in the design formula

### 2. Multiple Comparisons
For multiple treatment groups:
```python
# Example: comparing multiple mutants to control
my_stats_mut1 = DeseqStats(my_dds, contrast=["Condition", "mut1", "C"])
my_stats_mut2 = DeseqStats(my_dds, contrast=["Condition", "mut2", "C"])
```

### 3. Export for Further Analysis
```python
# Save normalized counts for other tools
norm_counts = my_dds.layers["normed_counts"]
norm_df = pd.DataFrame(norm_counts, index=my_dds.obs_names, columns=my_dds.var.index)
norm_df.to_csv("normalized_counts.csv")
```

## Why This Approach Works

This script exemplifies good research code because it:

- **Gets things done** without unnecessary complexity
- **Shows intermediate results** so you can verify each step
- **Creates publication-quality figures** with minimal code
- **Is easy to modify** for different experimental designs
- **Follows established workflows** from the bioinformatics community

The weighted random sampling in the volcano plot is particularly clever - it creates visually appealing plots while maintaining statistical rigor.

## Running Your Analysis

Execute the complete pipeline:

```bash
python code.py
```

The analysis typically takes 2-5 minutes depending on dataset size. You'll see progress messages and plots displayed during execution, with final files saved to your working directory.

## Citation

If you use this pipeline in your research, please cite:

- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
- **PyDESeq2**: Muzellec, B., Teyssier, M., Girard, E., et al. (2023). PyDESeq2: a Python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39(12), 2068-2069.

---

**Happy analyzing! 🧬📊**
