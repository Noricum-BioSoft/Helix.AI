#!/usr/bin/env Rscript
# DESeq2 differential expression analysis

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(DESeq2)
  library(dplyr)
})

# Parse command line arguments
option_list <- list(
  make_option(c("--count-matrix"), type="character", default=NULL,
              help="Path to count matrix CSV file"),
  make_option(c("--sample-metadata"), type="character", default=NULL,
              help="Path to sample metadata CSV file"),
  make_option(c("--design-formula"), type="character", default="~condition",
              help="Design formula for DESeq2"),
  make_option(c("--output-dir"), type="character", default=".",
              help="Output directory"),
  make_option(c("--alpha"), type="double", default=0.05,
              help="FDR threshold for significance")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Initialize results
results <- list(
  status = "success",
  output_files = c(),
  errors = c()
)

tryCatch({
  # Load count matrix
  count_matrix <- read.csv(opt$`count-matrix`, row.names=1, check.names=FALSE)
  
  # Load sample metadata
  sample_metadata <- read.csv(opt$`sample-metadata`, row.names=1, check.names=FALSE)
  
  # Ensure column order matches
  count_matrix <- count_matrix[, rownames(sample_metadata)]
  
  # Create DESeqDataSet
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = sample_metadata,
    design = as.formula(opt$`design-formula`)
  )
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  res <- results(dds, alpha=opt$alpha)
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  # Save results
  output_file <- file.path(opt$`output-dir`, "deseq2_results.csv")
  write.csv(res_df, output_file, row.names=FALSE)
  results$output_files <- c(results$output_files, output_file)
  
  # Summary statistics
  results$summary <- list(
    total_genes = nrow(res_df),
    significant_genes = sum(res_df$padj < opt$alpha, na.rm=TRUE),
    upregulated = sum(res_df$log2FoldChange > 0 & res_df$padj < opt$alpha, na.rm=TRUE),
    downregulated = sum(res_df$log2FoldChange < 0 & res_df$padj < opt$alpha, na.rm=TRUE)
  )
  
}, error = function(e) {
  results$status <<- "error"
  results$errors <<- c(results$errors, as.character(e))
})

# Output JSON
cat(toJSON(results, auto_unbox=TRUE))
