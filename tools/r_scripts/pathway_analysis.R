#!/usr/bin/env Rscript
# Perform pathway enrichment analysis using scPipeline

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
  make_option(c("--database"), type="character", default="Reactome",
              help="Pathway database: Reactome, GO, or KEGG")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

results <- list(status = "success")

tryCatch({
  # Load Seurat object
  seurat_obj <- readRDS(opt$`seurat-object`)
  
  # Perform pathway analysis
  pathways <- scPipeline::pathway_analysis(
    seurat_obj,
    ident.column = opt$`ident-column`,
    database = opt$database
  )
  
  results$pathways <- pathways
  results$summary <- list(
    total_pathways = nrow(pathways),
    significant_pathways = sum(pathways$p.adjust < 0.05, na.rm = TRUE),
    top_pathways = head(pathways[order(pathways$p.adjust), ], 20)
  )
  
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  
}, error = function(e) {
  results$status <- "error"
  results$error_message <- as.character(e)
  cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE))
  quit(status = 1)
})




