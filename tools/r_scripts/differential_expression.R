#!/usr/bin/env Rscript
# Perform differential expression analysis using scPipeline

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(scPipeline)
})

option_list <- list(
  make_option(c("--seurat-object"), type="character", default=NULL,
              help="Path to Seurat object RDS file"),
  make_option(c("--group1"), type="character", default=NULL,
              help="First group identifier"),
  make_option(c("--group2"), type="character", default=NULL,
              help="Second group identifier"),
  make_option(c("--ident-column"), type="character", default="seurat_clusters",
              help="Column name for group identity"),
  make_option(c("--min-pct"), type="double", default=0.1,
              help="Minimum percentage of cells expressing the gene"),
  make_option(c("--logfc-threshold"), type="double", default=0.25,
              help="Minimum log fold change threshold")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

results <- list(status = "success")

tryCatch({
  # Load Seurat object
  seurat_obj <- readRDS(opt$`seurat-object`)
  
  # Perform differential expression
  de_results <- scPipeline::differential_expression(
    seurat_obj,
    group1 = opt$group1,
    group2 = opt$group2,
    ident.column = opt$`ident-column`,
    min.pct = opt$`min-pct`,
    logfc.threshold = opt$`logfc-threshold`
  )
  
  results$differential_expression <- de_results
  results$summary <- list(
    total_genes = nrow(de_results),
    significant_genes = sum(de_results$p_val_adj < 0.05, na.rm = TRUE),
    upregulated = sum(de_results$avg_log2FC > 0 & de_results$p_val_adj < 0.05, na.rm = TRUE),
    downregulated = sum(de_results$avg_log2FC < 0 & de_results$p_val_adj < 0.05, na.rm = TRUE),
    top_genes = head(de_results[order(de_results$p_val_adj), ], 20)
  )
  
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})

