# R Scripts for Single-Cell Analysis

This directory contains R scripts that interface with the scPipeline package for single-cell RNA-seq analysis.

## Setup

### Install R Dependencies

Run the following in R to install all required packages:

```r
# Install CRAN packages
install.packages(c("optparse", "jsonlite", "Seurat", "dplyr", "magrittr", "rlang"))
install.packages("scPipeline")  # From CRAN: https://cran.r-project.org/web/packages/scPipeline/index.html

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "batchelor",
    "ReactomeGSA",
    "celldex",
    "SingleR",
    "SummarizedExperiment",
    "biomaRt"
))
```

### Verify Installation

Test that R and required packages are available:

```bash
# Check R is installed
Rscript --version

# Test package loading
Rscript -e "library(scPipeline); library(Seurat); cat('All packages loaded successfully\n')"
```

## Scripts

### Main Scripts

- **scpipeline_analysis.R**: Main workflow script that performs end-to-end analysis
- **find_markers.R**: Find marker genes for cell clusters
- **differential_expression.R**: Perform differential expression analysis
- **annotate_cell_types.R**: Annotate cell types using reference datasets
- **batch_correction.R**: Correct for batch effects
- **pathway_analysis.R**: Perform pathway enrichment analysis

## Usage

These scripts are called automatically by the Python `single_cell_analysis.py` tool. They accept command-line arguments and output JSON results.

### Example Manual Usage

```bash
# Run complete analysis
Rscript scpipeline_analysis.R \
  --data-file /path/to/data \
  --data-format 10x \
  --output-dir /path/to/output \
  --steps all

# Find markers
Rscript find_markers.R \
  --seurat-object seurat_object.rds \
  --ident-column seurat_clusters \
  --min-pct 0.25 \
  --logfc-threshold 0.25
```

## Important Notes

### scPipeline API Compatibility

The R scripts reference scPipeline functions that may need to be adapted based on the actual scPipeline API. The scPipeline package is a wrapper around Seurat, and function names may differ from what's written here.

**If you encounter function errors**, you may need to:

1. Check the actual scPipeline function names in the package documentation
2. Adapt the function calls to match the actual API
3. Or use Seurat functions directly if scPipeline doesn't provide the exact wrapper

### Example Adaptation

If `scPipeline::preprocess_data()` doesn't exist, you might need to use Seurat directly:

```r
# Instead of:
seurat_obj <- scPipeline::preprocess_data(seurat_obj, ...)

# Use Seurat directly:
seurat_obj <- Seurat::NormalizeData(seurat_obj)
seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
# etc.
```

## Troubleshooting

### R Script Execution Errors

1. **Package not found**: Ensure all packages are installed (see Setup above)
2. **Function not found**: Check scPipeline documentation and adapt function calls
3. **Permission denied**: Make scripts executable: `chmod +x *.R`
4. **JSON parsing errors**: Check that scripts output valid JSON

### Testing Scripts

Test individual scripts with sample data:

```bash
# Create a test Seurat object first (in R)
Rscript -e "
library(Seurat)
data <- Seurat::CreateSeuratObject(counts = matrix(rpois(1000, 1), nrow=100))
saveRDS(data, 'test_seurat.rds')
"

# Test marker finding
Rscript find_markers.R \
  --seurat-object test_seurat.rds \
  --ident-column seurat_clusters
```

## References

- [scPipeline CRAN](https://cran.r-project.org/web/packages/scPipeline/index.html)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [SingleR Documentation](https://bioconductor.org/packages/release/bioc/html/SingleR.html)




