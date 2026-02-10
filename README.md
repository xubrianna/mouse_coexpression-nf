# Mouse Coexpression Pipeline

A Nextflow pipeline for computing gene-gene coexpression matrices from pre-filtered single-cell RNA-seq data.

## Pipeline Overview

The pipeline starts with **already filtered h5ad files** (split by cell type) and consists of 3 main steps:

1. **GENERATE_GENE_LISTS**: Creates gene intersection lists across studies
   - Stringent: genes found in ALL studies per cell type
   - Less stringent: genes found in ≥20% of studies per cell type
   - Computes pairwise intersections between cell types

2. **SELECT_GENES_AND_CPM**: Filters to selected genes and normalizes
   - Selects genes based on strategy (intersection/union/cell-type-specific)
   - Adds missing genes with zero expression
   - Applies CPM normalization

3. **COMPUTE_COEXPRESSION**: Computes gene-gene correlation matrices (runs on SLURM)
   - Per-sample Pearson correlation
   - Rank normalization
   - Aggregation across samples
   - Outputs correlation matrices and edge lists

## Input Requirements

The pipeline expects **pre-filtered h5ad files** following the naming convention:
```
CellType_StudyName.h5ad
```

Examples:
- `Ast_Tasic2018.h5ad`
- `Mic_ROSMAP.h5ad`
- `Exc_Zeng2023.h5ad`

These files should already be:
- Filtered for genes (e.g., expressed in >2% of cells)
- Filtered for cells (e.g., bottom 2% removed)
- Split by cell type
- Gene names following ENSMUS format

## Quick Start

```bash
# Run pipeline (COMPUTE_COEXPRESSION automatically uses SLURM)
nextflow run main.nf --input_dir /path/to/filtered/h5ad

# Override parameters
nextflow run main.nf \
    --input_dir /path/to/raw/h5ad \
    --outdir my_results \
    --gene_selection_strategy intersection
```

## Parameters

### Input/Output
- `input_dir`: Directory with **filtered** h5ad files (default: data/0.raw_data/flt/)
- `outdir`: Output directory (default: results/)

### Step 1: Gene Lists
- `gene_frequency_threshold`: Minimum fraction of studies for "less stringent" genes (default: 0.2)

### Step 2: Gene Selection
- `gene_selection_strategy`: Strategy for selecting genes
  - `intersection`: Genes common to all cell types (most conservative)
  - `union`: All genes from any cell type
  - Cell type name (e.g., `Mic`): Genes from specific cell type

### Step 3: Coexpression
- `correlation_type`: Correlation method (default: pearson)
- `replace_nans`: Replace NaN values in correlation matrix (default: true)
- `min_cells`: Minimum number of cells (optional)

### Environment
- `conda_env`: Path to conda environment (default: /home/bxu/miniconda3/envs/python)

## Output Structure

```
results/
├── 01_gene_lists/        # Gene intersection lists
│   └── genes-lists-0.2.json
├── 02_cpm_normalized/    # CPM normalized files
│   ├── Ast_Study1.h5ad
│   └── ...
├── 03_coexpression/      # Correlation matrices
│   ├── CellType_Study/
│   │   ├── agg-mtx/
│   │   └── edge-list/
│   └── ...
└── pipeline_info/        # Execution reports
    ├── execution_timeline.html
    ├── execution_report.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```

## Requirements

- Nextflow ≥23.04.0
- SLURM cluster access (for COMPUTE_COEXPRESSION step)
- Conda environment with:
  - Python 3.10
  - scanpy
  - anndata
  - pandas
  - numpy
  - scipy
  - bottleneck
  - psutil

## Execution

The pipeline automatically runs:
- **GENERATE_GENE_LISTS** and **SELECT_GENES_AND_CPM**: locally
- **COMPUTE_COEXPRESSION**: on SLURM cluster

```bash
nextflow run main.nf --input_dir /path/to/filtered/h5ad
```

### SLURM Configuration

The COMPUTE_COEXPRESSION process is configured to use SLURM by default. To customize SLURM settings, update `nextflow.config`:

```groovy
process {
    withName: 'COMPUTE_COEXPRESSION' {
        executor = 'slurm'
        // Add partition if needed
        // queue = 'your_partition'
        // clusterOptions = '--account=your_account --time=24:00:00'
    }
}
```

## Examples

### Process all files with intersection strategy
```bash
nextflow run main.nf --gene_selection_strategy intersection
```

### Use union of all genes
```bash
nextflow run main.nf --gene_selection_strategy union
```

### Use microglia-specific genes
```bash
nextflow run main.nf --gene_selection_strategy Mic
```

### Custom filtering thresholds
```bashgene frequency threshold
```bash
nextflow run main.nf --gene_frequency_threshold 0.3
```

### Process specific input directory
```bash
nextflow run main.nf --input_dir /path/to/filtered/h5ad/files

## Troubleshooting

### Pipeline fails at GENERATE_GENE_LISTS
- Ensure all filtered files follow naming convention: `CellType_StudyName.h5ad`
- Check that multiple studies are present for each cell type

### Memory errors
- Increase memory in `nextflow.config` for the failing process
- Check available system resources

### Conda environment issues
- Ensure `conda_env` parameter points to correct environment
- Verify all required packages are installed:
  ```bash
  conda activate python
  python -c "import scanpy, anndata, scipy, bottleneck, psutil"
  ```

## Author

Brianna Xu

## Version

1.0.0
