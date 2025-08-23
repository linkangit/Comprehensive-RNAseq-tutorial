# Comprehensive RNA-seq Differential Expression Analysis Tutorial

This tutorial provides a complete walkthrough of RNA-seq differential expression analysis using Python, featuring PyDESeq2, Scanpy, and other bioinformatics libraries. The pipeline includes quality control, statistical analysis, visualization, and pathway enrichment analysis.

## Table of Contents
1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Input Data Requirements](#input-data-requirements)
4. [Installation](#installation)
5. [Pipeline Walkthrough](#pipeline-walkthrough)
6. [Output Files](#output-files)
7. [Interpretation Guide](#interpretation-guide)
8. [Troubleshooting](#troubleshooting)

## Overview

This analysis pipeline performs comprehensive differential expression analysis comparing control vs. mutant samples. The workflow includes:

- **Quality Control**: Sample-level statistics and visualization
- **Differential Expression**: Statistical analysis using DESeq2 algorithm
- **Visualization**: PCA, MA plots, volcano plots, and heatmaps
- **Functional Analysis**: Gene set enrichment analysis (GSEA)
- **Results Export**: Comprehensive result tables and summary statistics

### Key Features
- Robust error handling and logging
- Reproducible results with fixed random seeds
- Modular design for easy customization
- Professional-quality visualizations
- Comprehensive documentation

## Prerequisites

### Required Knowledge
- Basic understanding of RNA-seq data analysis concepts
- Familiarity with Python programming
- Understanding of differential expression analysis principles

### Required Software
- Python 3.8 or higher
- Git (for cloning the repository)

## Input Data Requirements

Based on the RNA-seq analysis script, here are the **exact inputs** that need to be in your working folder:

## Required Input Files

### 1. `read_counts_table.csv` (REQUIRED)
**Format:** CSV file with raw read counts
- **First column:** Must be named `gname` containing gene identifiers
- **Subsequent columns:** Sample names (must match the sample names in the script configuration)
- **Data:** Integer read counts

**Example structure:**
```csv
gname,Kelley_17,Kelley_18,Kelley_19,Kelley_20,Kelley_21,Kelley_22,Kelley_23,Kelley_24
AT1G01010,523,445,678,512,1205,1456,1123,1087
AT1G01020,89,92,76,85,45,38,52,61
AT1G01030,1234,1456,1123,1345,234,198,267,312
AT1G01040,45,67,23,89,156,234,189,201
...
```

### 2. `Arabidopsis_gene_annotation.tsv` (OPTIONAL)
**Format:** Tab-separated file with gene annotations
- **Required column:** `Nomenclature ID` (must match gene IDs in the count matrix)
- **Required column:** `Symbol` (gene symbols/names)
- **Additional columns:** Any other annotation information

**Example structure:**
```tsv
Nomenclature ID	Symbol	Description	Gene_type
AT1G01010	NAC001	NAC domain containing protein 1	protein_coding
AT1G01020	ARV1	ARV1 family protein	protein_coding
AT1G01030	NGA3	NGATHA3	protein_coding
AT1G01040	ASU1	ABSCISIC ACID INSENSITIVE4	protein_coding
...
```

## Sample Information (Defined in Script)

The sample information is **hardcoded in the script** in the `Config` class. You need to **modify this section** to match your actual samples:

```python
SAMPLE_INFO = {
    "Sample": ["Kelley_17", "Kelley_18", "Kelley_19", "Kelley_20",
               "Kelley_21", "Kelley_22", "Kelley_23", "Kelley_24"],
    "Condition": ["C", "C", "C", "C",
                  "mut", "mut", "mut", "mut"]
}
```

**Important Notes:**
- Sample names in `SAMPLE_INFO["Sample"]` must **exactly match** the column headers in `read_counts_table.csv`
- Use `"C"` for control samples and `"mut"` for mutant/treatment samples
- You can have different numbers of replicates, but you need at least 2 samples per condition

## Folder Structure

Your working directory should look like this:

```
your_project_folder/
├── code.py        # The analysis script
├── read_counts_table.csv          # REQUIRED: Your count matrix
├── Arabidopsis_gene_annotation.tsv # OPTIONAL: Gene annotations
└── results/                       # Will be created automatically
```

## What You Need to Customize

### 1. Update Sample Information
Edit the `Config` class in the script to match your samples:

```python
class Config:
    # ... other parameters ...
    
    # UPDATE THIS SECTION FOR YOUR DATA
    SAMPLE_INFO = {
        "Sample": ["Your_Sample_1", "Your_Sample_2", "Your_Sample_3", 
                   "Your_Sample_4", "Your_Sample_5", "Your_Sample_6"],
        "Condition": ["control", "control", "control", 
                      "treatment", "treatment", "treatment"]
    }
```

### 2. Update File Paths (if different)
If your files have different names, update these in the `Config` class:

```python
class Config:
    # File paths - UPDATE IF YOUR FILES HAVE DIFFERENT NAMES
    READ_COUNTS_FILE = 'your_counts_file.csv'
    GENE_ANNOTATION_FILE = 'your_annotation_file.tsv'
    OUTPUT_DIR = 'results'
```

## Quick Checklist Before Running

✅ **File naming:**
- [ ] Count matrix file is named `read_counts_table.csv` (or update script)
- [ ] Annotation file is named `Arabidopsis_gene_annotation.tsv` (or update script)

✅ **File format:**
- [ ] Count matrix has `gname` as first column name
- [ ] Count matrix has sample names matching your `SAMPLE_INFO`
- [ ] Annotation file has `Nomenclature ID` and `Symbol` columns (if using)

✅ **Sample information:**
- [ ] Updated `SAMPLE_INFO` in the script to match your samples
- [ ] Sample names match exactly between count matrix and script
- [ ] At least 2 samples per condition group

✅ **File location:**
- [ ] All files are in the same directory as the script
- [ ] No extra spaces or special characters in file names

## Installation

### 1. Clone the Repository
```bash
git clone <your-repository-url>
cd rna-seq-analysis
```

### 2. Create Virtual Environment (Recommended)
```bash
python -m venv rna_env
source rna_env/bin/activate  # On Windows: rna_env\Scripts\activate
```

### 3. Install Required Packages
```bash
# Core dependencies
pip install pandas numpy matplotlib seaborn
pip install scanpy anndata
pip install pydeseq2

# Optional but recommended
pip install gseapy adjustText
```

### 4. Verify Installation
```python
python -c "import pydeseq2, scanpy, pandas; print('All packages installed successfully!')"
```

## Pipeline Walkthrough

### Step 1: Configuration and Setup

The analysis begins by setting up configuration parameters and creating output directories:

```python
class Config:
    # File paths
    READ_COUNTS_FILE = 'read_counts_table.csv'
    GENE_ANNOTATION_FILE = 'Arabidopsis_gene_annotation.tsv'
    OUTPUT_DIR = 'results'
    
    # Analysis parameters
    MIN_BASE_MEAN = 10          # Minimum mean expression for inclusion
    PADJ_THRESHOLD = 0.1        # Adjusted p-value threshold
    LOG2FC_THRESHOLD = 0.5      # Log2 fold-change threshold
```

**What happens here:**
- Sets file paths and analysis parameters
- Creates reproducible analysis with fixed random seed
- Initializes logging system for progress tracking

### Step 2: Sample Information Preparation

```python
def prepare_sample_info(config):
    sample_info = pd.DataFrame(config.SAMPLE_INFO)
    ctrl_list = sample_info.loc[sample_info['Condition'] == 'C', 'Sample'].tolist()
    variant_list = sample_info.loc[sample_info['Condition'] == 'mut', 'Sample'].tolist()
```

**What happens here:**
- Creates a structured DataFrame with sample metadata
- Identifies control and mutant sample groups
- Validates that both groups have samples
- Prepares data for downstream statistical analysis

### Step 3: Data Loading and Quality Control

```python
def load_and_clean_counts(config, all_samples):
    read_df = pd.read_csv(config.READ_COUNTS_FILE)
    read_df = read_df.set_index('gname')
    read_df = read_df[all_samples]  # Filter to relevant samples
    read_df = read_df[read_df.sum(axis=1) > 0]  # Remove zero-count genes
```

**What happens here:**
- Loads the count matrix from CSV file
- Sets gene names as row indices
- Filters to include only specified samples
- Removes genes with zero counts across all samples
- Validates data integrity and sample presence

### Step 4: Quality Control Metrics

```python
def calculate_qc_metrics(read_df, sample_info):
    qc_metrics = pd.DataFrame({
        'Total_Reads': read_df.sum(axis=0),
        'Detected_Genes': (read_df > 0).sum(axis=0),
        'Condition': sample_info.set_index('Sample')['Condition']
    })
```

**What happens here:**
- Calculates total read counts per sample
- Counts detected genes (with >0 reads) per sample
- Creates visualizations comparing control vs. mutant groups
- Helps identify potential batch effects or outlier samples

**Generated Output:**
- `qc_metrics.png` - Bar plots showing read counts and gene detection rates

### Step 5: DESeq2 Differential Expression Analysis

```python
def run_deseq2_analysis(read_df, sample_info):
    transposed_df = read_df.T  # DESeq2 expects samples as rows
    my_dds = DeseqDataSet(
        counts=transposed_df,
        metadata=sample_info_deseq,
        design_factors=["Condition"]
    )
    my_dds.deseq2()
```

**What happens here:**
- Transposes count matrix (DESeq2 expects samples in rows, genes in columns)
- Creates DESeqDataSet object with experimental design
- Runs DESeq2 algorithm for:
  - Size factor estimation (normalization)
  - Dispersion estimation
  - Statistical testing using negative binomial distribution
- Handles low-count filtering automatically

### Step 6: Results Extraction and Gene Annotation

```python
def extract_de_results(my_dds, config):
    my_stats = DeseqStats(my_dds, contrast=["Condition", "mut", "C"])
    results_df = my_stats.results_df
```

**What happens here:**
- Defines statistical contrast (mutant vs. control)
- Extracts results including:
  - `log2FoldChange`: Effect size
  - `pvalue`: Raw p-values
  - `padj`: Benjamini-Hochberg adjusted p-values
  - `baseMean`: Mean normalized expression
- Merges with gene annotations if available

### Step 7: Statistical Filtering

```python
def filter_and_identify_significant(results_df, config):
    results_df = results_df[results_df.baseMean >= config.MIN_BASE_MEAN]
    significant_mask = (results_df.padj < config.PADJ_THRESHOLD) & \
                      (abs(results_df.log2FoldChange) > config.LOG2FC_THRESHOLD)
    significant_hits = results_df[significant_mask]
```

**What happens here:**
- Filters out low-expression genes (baseMean < 10)
- Applies statistical significance criteria:
  - Adjusted p-value < 0.1 (FDR < 10%)
  - Absolute log2 fold-change > 0.5 (1.4-fold change)
- Counts upregulated and downregulated genes

### Step 8: Visualization Generation

#### Principal Component Analysis (PCA)
```python
def create_pca_plot(my_dds, config):
    adata = anndata.AnnData(X=norm_data, obs=my_dds.obs)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pl.pca(adata, color='Condition')
```

**Purpose:** Visualizes overall sample relationships and potential batch effects
**Output:** `pca_plot.png`

#### MA Plot
```python
def create_ma_plot(results_df, config):
    plt.scatter(x=results_df['baseMean'], y=results_df['log2FoldChange'])
    plt.xscale('log')
```

**Purpose:** Shows relationship between gene expression level and fold-change
**Output:** `ma_plot.png`

#### Volcano Plot
```python
def create_volcano_plot(results_df, config):
    plot_data["nlog10_padj"] = -np.log10(plot_data["padj"])
    plt.scatter(plot_data['log2FoldChange'], plot_data['nlog10_padj'])
```

**Purpose:** Visualizes statistical significance vs. biological significance
**Output:** `volcano_plot.png`

#### Expression Heatmaps
```python
def create_expression_heatmaps(my_dds, significant_hits, config):
    sns.clustermap(sig_data.T, z_score=0, cmap="RdBu_r")
```

**Purpose:** Shows expression patterns of significant genes across samples
**Outputs:** 
- `heatmap_all_significant.png` - All significant genes
- `heatmap_top20.png` - Top 10 up + top 10 down regulated genes

### Step 9: Pathway Enrichment Analysis (Optional)

```python
def run_gsea_analysis(significant_hits, config):
    enr = gp.enrichr(
        gene_list=gene_list['Symbol'].tolist(),
        gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
        organism='plant'
    )
```

**What happens here:**
- Uses significantly differentially expressed genes
- Queries Gene Ontology and KEGG pathway databases
- Identifies overrepresented biological processes and pathways
- Saves results for functional interpretation

### Step 10: Results Export

```python
def save_results(results_df, significant_hits, config):
    results_df.to_csv(f"{config.OUTPUT_DIR}/DE_results_all.csv")
    significant_hits.to_csv(f"{config.OUTPUT_DIR}/DE_results_significant.csv")
```

**What happens here:**
- Saves complete results table with all genes
- Creates separate file with only significant genes
- Generates analysis summary with key statistics
- Organizes all outputs in structured directory

## Output Files

After successful completion, the `results/` directory will contain:

### Visualization Files
- **`qc_metrics.png`** - Quality control plots showing read counts and detected genes per sample
- **`pca_plot.png`** - Principal component analysis showing sample clustering
- **`correlation_heatmap.png`** - Sample-to-sample correlation matrix
- **`ma_plot.png`** - MA plot showing expression level vs. fold-change relationship
- **`volcano_plot.png`** - Volcano plot highlighting significant genes
- **`heatmap_all_significant.png`** - Clustered heatmap of all significant genes
- **`heatmap_top20.png`** - Heatmap of top 10 up- and down-regulated genes

### Data Files
- **`DE_results_all.csv`** - Complete differential expression results for all genes
- **`DE_results_significant.csv`** - Results filtered to significant genes only
- **`analysis_summary.txt`** - Summary statistics and parameters used
- **`gsea/`** - Directory containing pathway enrichment results (if available)

### Key Columns in Results Files

| Column | Description |
|--|-|
| `baseMean` | Mean normalized expression across all samples |
| `log2FoldChange` | Log2 fold-change (mutant vs. control) |
| `lfcSE` | Standard error of log2 fold-change |
| `stat` | Test statistic |
| `pvalue` | Raw p-value |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `Symbol` | Gene symbol (from annotation file) |

## Interpretation Guide

### Understanding the Results

#### 1. Quality Control Assessment
- **Read count distribution**: Should be similar across samples
- **Detected genes**: Should be consistent within condition groups
- **PCA plot**: Samples should cluster by condition, not by batch

#### 2. Statistical Results
- **Significant genes**: padj < 0.1 AND |log2FoldChange| > 0.5
- **Upregulated**: Positive log2FoldChange values (higher in mutant)
- **Downregulated**: Negative log2FoldChange values (lower in mutant)

#### 3. Visualization Interpretation

**MA Plot:**
- X-axis: Average expression level
- Y-axis: Log2 fold-change
- Red points: Significant genes
- Horizontal lines: Fold-change thresholds

**Volcano Plot:**
- X-axis: Log2 fold-change (effect size)
- Y-axis: -log10(adjusted p-value) (significance)
- Upper right: Significantly upregulated
- Upper left: Significantly downregulated

**Heatmaps:**
- Rows: Genes
- Columns: Samples
- Colors: Z-score normalized expression
- Clustering shows gene co-expression patterns

### Typical Results Interpretation

1. **No significant genes found:**
   - Check if biological difference is expected
   - Consider less stringent thresholds
   - Verify sample labeling

2. **Many significant genes (>5000):**
   - Strong biological effect
   - Consider more stringent thresholds
   - Check for batch effects

3. **Good results (100-2000 significant genes):**
   - Reasonable biological response
   - Proceed with pathway analysis
   - Validate key genes

## Troubleshooting

### Common Issues and Solutions

#### 1. Import Errors
**Problem:** `ModuleNotFoundError: No module named 'pydeseq2'`
**Solution:**
```bash
pip install pydeseq2
# or if using conda:
conda install -c bioconda pydeseq2
```

#### 2. File Not Found Errors
**Problem:** `FileNotFoundError: read_counts_table.csv`
**Solution:**
- Verify file exists in working directory
- Check file name spelling
- Use absolute path if necessary

#### 3. Sample Mismatch
**Problem:** `Missing samples in read count data`
**Solution:**
- Verify sample names match exactly between count matrix and sample info
- Check for typos or extra spaces
- Update sample names in Config class

#### 4. Memory Issues
**Problem:** `MemoryError` with large datasets
**Solution:**
- Increase available RAM
- Filter to protein-coding genes only
- Process in smaller batches

#### 5. No Significant Genes
**Problem:** Empty significant results
**Solution:**
- Check if biological effect is expected
- Reduce stringency: increase padj threshold to 0.2
- Reduce fold-change threshold to 0.25
- Examine MA and volcano plots for overall patterns

#### 6. PyDESeq2 Version Issues
**Problem:** `AttributeError: 'DeseqDataSet' object has no attribute 'layers'`
**Solution:**
```bash
pip install --upgrade pydeseq2
# Check version compatibility in documentation
```

### Performance Optimization

For large datasets (>50,000 genes), consider:

1. **Pre-filtering:**
```python
# Filter genes with low counts across all samples
keep = read_df.sum(axis=1) >= 10
read_df = read_df[keep]
```

2. **Subset analysis:**
```python
# Focus on protein-coding genes
protein_coding_genes = gene_info[gene_info['Gene_type'] == 'protein_coding'].index
read_df = read_df.loc[read_df.index.isin(protein_coding_genes)]
```

### Getting Help

If you encounter issues not covered here:

1. Check the [PyDESeq2 documentation](https://pydeseq2.readthedocs.io/)
2. Review [Scanpy tutorials](https://scanpy.readthedocs.io/)
3. Open an issue on this repository with:
   - Error message
   - Python version
   - Package versions (`pip list`)
   - Minimal reproducible example

## Running the Analysis

To execute the complete pipeline:

```bash
# Ensure you're in the directory with your data files
ls  # Should show: read_counts_table.csv, Arabidopsis_gene_annotation.tsv

# Run the analysis
python code.py

# Check results
ls results/  # View generated outputs
```

The analysis typically takes 2-10 minutes depending on dataset size and computer specifications.



## Citation

If you use this pipeline in your research, please cite the underlying methods:

- **DESeq2**: Love, M.I., Huber, W., Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15, 550.
- **PyDESeq2**: Muzellec, B., Teyssier, M., Girard, E., et al. (2023). PyDESeq2: a Python package for bulk RNA-seq differential expression analysis. Bioinformatics, 39(12), 2068-2069.



**Happy analyzing! 🧬📊**
