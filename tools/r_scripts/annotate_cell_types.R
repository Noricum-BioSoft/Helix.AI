#!/usr/bin/env Rscript
# Annotate cell types using scPipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(scPipeline)
})

option_list <- list(
  make_option(c("--seurat-object"), type="character", default=NULL,
              help="Path to Seurat object RDS file"),
  make_option(c("--reference"), type="character", default="HumanPrimaryCellAtlasData",
              help="Reference dataset name from celldex"),
  make_option(c("--method"), type="character", default="SingleR",
              help="Annotation method: SingleR or scType")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

results <- list(status = "success")

tryCatch({
  # Load Seurat object
  seurat_obj <- readRDS(opt$`seurat-object`)
  
  # Annotate cell types
  annotations <- scPipeline::annotate_cell_types(
    seurat_obj,
    reference = opt$reference,
    method = opt$method
  )
  
  results$annotations <- annotations
  results$summary <- list(
    total_cells = nrow(annotations),
    cell_types = unique(annotations$cell_type),
    cell_type_counts = as.list(table(annotations$cell_type))
  )
  
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})




