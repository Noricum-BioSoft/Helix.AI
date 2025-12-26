#!/usr/bin/env Rscript
# Find marker genes using scPipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(scPipeline)
})

option_list <- list(
  make_option(c("--seurat-object"), type="character", default=NULL,
              help="Path to Seurat object RDS file"),
  make_option(c("--ident-column"), type="character", default="seurat_clusters",
              help="Column name for cluster identity"),
  make_option(c("--min-pct"), type="double", default=0.25,
              help="Minimum percentage of cells expressing the gene"),
  make_option(c("--logfc-threshold"), type="double", default=0.25,
              help="Minimum log fold change threshold"),
  make_option(c("--only-pos"), action="store_true", default=FALSE,
              help="Only return positive markers")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

results <- list(status = "success")

tryCatch({
  # Load Seurat object
  seurat_obj <- readRDS(opt$`seurat-object`)
  
  # Find markers
  markers <- scPipeline::find_markers(
    seurat_obj,
    ident.1 = NULL,  # Find markers for all clusters
    ident.column = opt$`ident-column`,
    min.pct = opt$`min-pct`,
    logfc.threshold = opt$`logfc-threshold`,
    only.pos = opt$`only-pos`
  )
  
  results$markers <- markers
  results$summary <- list(
    total_markers = nrow(markers),
    clusters = unique(markers$cluster),
    top_markers = head(markers, 20)
  )
  
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})




