#!/usr/bin/env Rscript
# Main scPipeline analysis script
# This script performs end-to-end single-cell analysis using scPipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(scPipeline)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--data-file"), type="character", default=NULL,
              help="Path to input data file"),
  make_option(c("--data-format"), type="character", default="10x",
              help="Data format: 10x, h5, csv, or seurat"),
  make_option(c("--output-dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--steps"), type="character", default="all",
              help="Comma-separated list of steps: preprocessing,markers,differential,pathways,annotation,batch_correction,all"),
  make_option(c("--min-cells"), type="integer", default=3,
              help="Minimum number of cells expressing a gene"),
  make_option(c("--min-features"), type="integer", default=200,
              help="Minimum number of features per cell"),
  make_option(c("--max-mito"), type="double", default=0.2,
              help="Maximum mitochondrial percentage"),
  make_option(c("--nfeatures"), type="integer", default=2000,
              help="Number of variable features"),
  make_option(c("--dims"), type="character", default="1:30",
              help="Dimensions for PCA/UMAP"),
  make_option(c("--resolution"), type="double", default=0.5,
              help="Clustering resolution")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Parse steps
steps <- strsplit(opt$steps, ",")[[1]]
if ("all" %in% steps) {
  steps <- c("preprocessing", "markers", "differential", "pathways", "annotation")
}

# Create output directory
dir.create(opt$`output-dir`, showWarnings = FALSE, recursive = TRUE)

# Initialize results list
results <- list(
  status = "success",
  steps_completed = c(),
  output_files = c(),
  errors = c()
)

tryCatch({
  # Load data
  if (is.null(opt$`data-file`)) {
    stop("Data file is required")
  }
  
  seurat_obj <- NULL
  
  if (opt$`data-format` == "10x") {
    # Load 10x data
    seurat_obj <- scPipeline::load_10x_data(opt$`data-file`)
  } else if (opt$`data-format` == "h5") {
    # Load H5 file
    seurat_obj <- scPipeline::load_h5_data(opt$`data-file`)
  } else if (opt$`data-format` == "seurat") {
    # Load Seurat object
    seurat_obj <- readRDS(opt$`data-file`)
  } else {
    stop(paste("Unsupported data format:", opt$`data-format`))
  }
  
  # Save Seurat object
  seurat_path <- file.path(opt$`output-dir`, "seurat_object.rds")
  saveRDS(seurat_obj, seurat_path)
  results$output_files <- c(results$output_files, seurat_path)
  
  # Preprocessing step
  if ("preprocessing" %in% steps || "all" %in% steps) {
    cat("Performing preprocessing...\n")
    seurat_obj <- scPipeline::preprocess_data(
      seurat_obj,
      min.cells = opt$`min-cells`,
      min.features = opt$`min-features`,
      max.mito = opt$`max-mito`,
      nfeatures = opt$`nfeatures`
    )
    results$steps_completed <- c(results$steps_completed, "preprocessing")
  }
  
  # Find markers
  if ("markers" %in% steps || "all" %in% steps) {
    cat("Finding marker genes...\n")
    markers <- scPipeline::find_markers(seurat_obj)
    markers_path <- file.path(opt$`output-dir`, "markers.csv")
    write.csv(markers, markers_path, row.names = FALSE)
    results$output_files <- c(results$output_files, markers_path)
    results$markers_summary <- list(
      total_markers = nrow(markers),
      top_markers = head(markers, 10)
    )
    results$steps_completed <- c(results$steps_completed, "markers")
  }
  
  # Differential expression
  if ("differential" %in% steps || "all" %in% steps) {
    cat("Performing differential expression analysis...\n")
    de_results <- scPipeline::differential_expression(seurat_obj)
    de_path <- file.path(opt$`output-dir`, "differential_expression.csv")
    write.csv(de_results, de_path, row.names = FALSE)
    results$output_files <- c(results$output_files, de_path)
    results$de_summary <- list(
      total_genes = nrow(de_results),
      significant_genes = sum(de_results$p_val_adj < 0.05, na.rm = TRUE)
    )
    results$steps_completed <- c(results$steps_completed, "differential")
  }
  
  # Pathway analysis
  if ("pathways" %in% steps || "all" %in% steps) {
    cat("Performing pathway enrichment...\n")
    pathways <- scPipeline::pathway_analysis(seurat_obj)
    pathways_path <- file.path(opt$`output-dir`, "pathways.csv")
    write.csv(pathways, pathways_path, row.names = FALSE)
    results$output_files <- c(results$output_files, pathways_path)
    results$pathways_summary <- list(
      total_pathways = nrow(pathways),
      top_pathways = head(pathways, 10)
    )
    results$steps_completed <- c(results$steps_completed, "pathways")
  }
  
  # Cell type annotation
  if ("annotation" %in% steps || "all" %in% steps) {
    cat("Annotating cell types...\n")
    annotations <- scPipeline::annotate_cell_types(seurat_obj)
    annotation_path <- file.path(opt$`output-dir`, "cell_annotations.csv")
    write.csv(annotations, annotation_path, row.names = FALSE)
    results$output_files <- c(results$output_files, annotation_path)
    results$annotation_summary <- list(
      cell_types = unique(annotations$cell_type),
      total_cells = nrow(annotations)
    )
    results$steps_completed <- c(results$steps_completed, "annotation")
  }
  
  # Save final Seurat object
  final_seurat_path <- file.path(opt$`output-dir`, "seurat_object_final.rds")
  saveRDS(seurat_obj, final_seurat_path)
  results$output_files <- c(results$output_files, final_seurat_path)
  
  # Convert results to JSON
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})




