#!/usr/bin/env Rscript

# General DESeq2 Differential Expression Analysis Script
# Usage: Rscript deseq2_analysis.R <count_matrix> <metadata> <output_dir> <condition_col> <control_group>
# Arguments:
#   <count_matrix>: Gene expression count matrix file (rows: genes, columns: samples)
#   <metadata>: Sample metadata file (contains group information)
#   <output_dir>: Output directory for results
#   <condition_col>: Column name in metadata with group information
#   <control_group>: Name of the control group

# Load required packages
load_packages <- function() {
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  required_pkgs <- c("DESeq2", "ggplot2", "pheatmap", "dplyr", "tidyr")
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
}

# Parse command line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 5) {
    cat("Error: Incorrect number of arguments\n")
    cat("Correct usage: Rscript deseq2_analysis.R <count_matrix> <metadata> <output_dir> <condition_col> <control_group>\n")
    quit(status = 1)
  }
  
  list(
    count_file = args[1],
    meta_file = args[2],
    out_dir = args[3],
    condition_col = args[4],
    control = args[5]
  )
}

# Create output directories
create_output_dir <- function(out_dir) {
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  # Create subdirectory for plots
  dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)
}

# Load and prepare data
load_data <- function(count_file, meta_file) {
  # Read count matrix
  counts <- read.delim(count_file, row.names = 1, check.names = FALSE)
  
  # Read metadata
  meta <- read.delim(meta_file, row.names = 1)
  
  # Ensure sample names match between counts and metadata
  common_samples <- intersect(colnames(counts), rownames(meta))
  if (length(common_samples) == 0) {
    stop("Sample names do not match between count matrix and metadata")
  }
  
  list(
    counts = counts[, common_samples, drop = FALSE],
    meta = meta[common_samples, , drop = FALSE]
  )
}

# Main DESeq2 analysis pipeline
run_deseq2 <- function(counts, meta, condition_col, control_group) {
  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = as.formula(paste("~", condition_col))
  )
  
  # Set control group as reference
  dds[[condition_col]] <- relevel(dds[[condition_col]], ref = control_group)
  
  # Filter low-expression genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Perform differential expression analysis
  dds <- DESeq(dds)
  
  # Get all comparison results
  contrasts <- resultsNames(dds)
  res_list <- lapply(contrasts, function(contrast) {
    results(dds, name = contrast, alpha = 0.05)
  })
  names(res_list) <- contrasts
  
  list(dds = dds, results = res_list)
}

# Process and save results
process_results <- function(res_list, out_dir) {
  for (contrast in names(res_list)) {
    # Convert to data frame
    res_df <- as.data.frame(res_list[[contrast]]) %>%
      rownames_to_column("gene") %>%
      arrange(padj)
    
    # Save complete results
    write.csv(
      res_df,
      file.path(out_dir, paste0("deseq2_results_", gsub("/", "_vs_", contrast), ".csv")),
      row.names = FALSE
    )
    
    # Save significantly differentially expressed genes (padj < 0.05 & |log2FoldChange| > 1)
    sig_res <- filter(res_df, padj < 0.05, abs(log2FoldChange) > 1)
    write.csv(
      sig_res,
      file.path(out_dir, paste0("significant_genes_", gsub("/", "_vs_", contrast), ".csv")),
      row.names = FALSE
    )
  }
}

# Generate visualizations
generate_plots <- function(dds, res_list, out_dir, condition_col) {
  # Extract first comparison result for visualization
  contrast <- names(res_list)[1]
  res <- res_list[[1]]
  
  # 1. MA plot
  png(file.path(out_dir, "plots", paste0("ma_plot_", gsub("/", "_vs_", contrast), ".png")),
      width = 800, height = 700)
  plotMA(res, ylim = c(-2, 2), alpha = 0.05)
  dev.off()
  
  # 2. Volcano plot
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene") %>%
    mutate(
      significance = case_when(
        padj < 0.05 & abs(log2FoldChange) > 1 ~ "Significant",
        padj < 0.05 ~ "Significant (FC < 2)",
        abs(log2FoldChange) > 1 ~ "FC > 2 (ns)",
        TRUE ~ "Not significant"
      )
    )
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("Significant" = "red", "Significant (FC < 2)" = "blue", 
                                 "FC > 2 (ns)" = "green", "Not significant" = "gray")) +
    theme_minimal() +
    labs(title = "Differentially Expressed Genes", x = "log2(Fold Change)", y = "-log10(Adjusted p-value)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(
    file.path(out_dir, "plots", paste0("volcano_plot_", gsub("/", "_vs_", contrast), ".png")),
    plot = p, width = 10, height = 8, dpi = 300
  )
  
  # 3. Sample clustering heatmap (using rlog transformed data)
  rld <- rlog(dds, blind = FALSE)
  sample_dist <- dist(t(assay(rld)))
  sample_dist_matrix <- as.matrix(sample_dist)
  
  png(file.path(out_dir, "plots", "sample_clustering.png"), width = 900, height = 800)
  pheatmap(sample_dist_matrix,
           annotation_col = colData(dds)[, condition_col, drop = FALSE],
           main = "Sample Clustering Heatmap")
  dev.off()
}

# Main function
main <- function() {
  # Initialization
  load_packages()
  args <- parse_args()
  create_output_dir(args$out_dir)
  
  cat("Starting DESeq2 differential expression analysis...\n")
  cat("Input count matrix:", args$count_file, "\n")
  cat("Sample metadata:", args$meta_file, "\n")
  cat("Output directory:", args$out_dir, "\n")
  cat("Condition column:", args$condition_col, "\n")
  cat("Control group:", args$control, "\n\n")
  
  # Analysis pipeline
  data <- load_data(args$count_file, args$meta_file)
  deseq_result <- run_deseq2(data$counts, data$meta, args$condition_col, args$control)
  process_results(deseq_result$results, args$out_dir)
  generate_plots(deseq_result$dds, deseq_result$results, args$out_dir, args$condition_col)
  
  cat("\nAnalysis completed! Results saved to", args$out_dir, "\n")
}

# Run main function
main()
