#!/usr/bin/env Rscript
# Perform batch correction using scPipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(scPipeline)
})

option_list <- list(
  make_option(c("--seurat-object"), type="character", default=NULL,
              help="Path to Seurat object RDS file"),
  make_option(c("--batch-column"), type="character", default="batch",
              help="Column name for batch information"),
  make_option(c("--method"), type="character", default="harmony",
              help="Batch correction method: harmony, seurat, or batchelor")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

results <- list(status = "success")

tryCatch({
  # Load Seurat object
  seurat_obj <- readRDS(opt$`seurat-object`)
  
  # Perform batch correction
  corrected_obj <- scPipeline::batch_correction(
    seurat_obj,
    batch.column = opt$`batch-column`,
    method = opt$method
  )
  
  # Save corrected object
  output_path <- sub("\\.rds$", "_corrected.rds", opt$`seurat-object`)
  saveRDS(corrected_obj, output_path)
  
  results$corrected_object_path <- output_path
  results$summary <- list(
    method = opt$method,
    batches = unique(seurat_obj@meta.data[[opt$`batch-column`]]),
    total_cells = ncol(corrected_obj)
  )
  
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})




