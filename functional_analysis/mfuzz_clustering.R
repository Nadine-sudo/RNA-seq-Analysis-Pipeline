#!/usr/bin/env Rscript

# Generalized Mfuzz Gene Expression Clustering Pipeline
# Usage: Rscript mfuzz_clustering.R
# All parameters are defined in the 'User Parameters' section below

# ------------------------------------------------------------------------------
# User Parameters - MODIFY THESE VALUES ACCORDING TO YOUR ANALYSIS
# ------------------------------------------------------------------------------
params <- list(
  # Input/output settings
  input_dir = "path/to/expression_data",         # Directory containing expression file
  output_dir = "./mfuzz_results",          # Directory to save results
  expression_file = "expression_values.csv", # Name of expression file
  gene_id_column = "GeneSymbol",           # Column name with gene identifiers
  separator = ",",                         # Column separator in input file (e.g., "," or "\t")
  
  # Clustering parameters
  num_clusters = 5,                        # Number of clusters to generate
  na_threshold = 0.25,                     # Maximum proportion of NAs allowed (0-1)
  min_std = 0,                             # Minimum standard deviation for filtering
  na_fill_method = "knn",                  # Method for filling NAs ("knn" or "mean")
  fuzzifier_m = NULL,                      # Fuzzifier parameter (NULL = auto-estimate)
  
  # Plot settings
  plot_width = 1000,                       # Plot width in pixels
  plot_height = 800,                       # Plot height in pixels
  plot_resolution = 150,                   # Plot resolution in dpi
  plot_rows = 2,                           # Number of rows in cluster plot
  
  # Sample group definitions
  # Format: list of groups, each with "name" and "time_points"
  # Time points: named list where each entry has column names/indices for that time point
  sample_groups = list(
    list(
      name = "control",
      time_points = list(
        t0 = c("ctrl_t0_rep1", "ctrl_t0_rep2"),  # Columns for time point 0
        t1 = c("ctrl_t1_rep1", "ctrl_t1_rep2"),  # Columns for time point 1
        t2 = c("ctrl_t2_rep1", "ctrl_t2_rep2")   # Columns for time point 2
      )
    ),
    list(
      name = "treatment",
      time_points = list(
        t0 = c("trt_t0_rep1", "trt_t0_rep2"),    # Columns for time point 0
        t1 = c("trt_t1_rep1", "trt_t1_rep2"),    # Columns for time point 1
        t2 = c("trt_t2_rep1", "trt_t2_rep2")     # Columns for time point 2
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------
load_packages <- function() {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", dependencies = TRUE)
  }
  
  required_pkgs <- c("Mfuzz", "clusterProfiler", "org.Mm.eg.db", 
                    "plyr", "tidyverse", "rtracklayer")
  
  # Install missing packages
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg %in% c("Mfuzz", "clusterProfiler", "org.Mm.eg.db", "rtracklayer")) {
        BiocManager::install(pkg, update = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE, repos = "https://cran.r-project.org")
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# ------------------------------------------------------------------------------
# Validate parameters
# ------------------------------------------------------------------------------
validate_parameters <- function(params) {
  # Check required input parameters
  if (!dir.exists(params$input_dir)) {
    stop("Input directory does not exist: ", params$input_dir)
  }
  
  expr_file <- file.path(params$input_dir, params$expression_file)
  if (!file.exists(expr_file)) {
    stop("Expression file not found: ", expr_file)
  }
  
  # Validate sample groups
  if (length(params$sample_groups) == 0) {
    stop("No sample groups defined in parameters")
  }
  
  for (group in params$sample_groups) {
    if (is.null(group$name) || is.null(group$time_points)) {
      stop("Sample groups must contain 'name' and 'time_points'")
    }
    
    if (length(group$time_points) == 0) {
      stop("Group '", group$name, "' has no time points defined")
    }
  }
  
  # Validate NA fill method
  if (!params$na_fill_method %in% c("knn", "mean")) {
    stop("NA fill method must be either 'knn' or 'mean'")
  }
  
  # Validate numeric parameters
  if (params$num_clusters < 2) {
    stop("Number of clusters must be at least 2")
  }
  
  if (params$na_threshold < 0 || params$na_threshold > 1) {
    stop("NA threshold must be between 0 and 1")
  }
  
  if (params$min_std < 0) {
    stop("Minimum standard deviation cannot be negative")
  }
}

# ------------------------------------------------------------------------------
# Create output directories
# ------------------------------------------------------------------------------
setup_directories <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)
  dir.create(file.path(output_dir, "clusters"), showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# Load and preprocess data
# ------------------------------------------------------------------------------
load_data <- function(params) {
  # Define input file path
  expr_file <- file.path(params$input_dir, params$expression_file)
  
  # Load data with specified separator
  expr_data <- read.table(expr_file, header = TRUE, sep = params$separator, stringsAsFactors = FALSE)
  
  # Validate gene ID column exists
  if (!params$gene_id_column %in% colnames(expr_data)) {
    stop("Gene ID column not found in expression data: ", params$gene_id_column)
  }
  
  expr_data
}

# ------------------------------------------------------------------------------
# Prepare group-specific datasets with time points
# ------------------------------------------------------------------------------
prepare_group_data <- function(expr_data, params) {
  gene_id_col <- params$gene_id_column
  
  # Process each sample group defined in parameters
  lapply(params$sample_groups, function(group) {
    # Initialize with gene identifiers
    group_df <- expr_data[, gene_id_col, drop = FALSE]
    colnames(group_df) <- "Gene"  # Standardize column name
    
    # Calculate mean expression for each time point
    time_points <- group$time_points
    for (tp_name in names(time_points)) {
      # Get column indices or names
      cols <- time_points[[tp_name]]
      
      # Check if columns exist
      if (is.numeric(cols)) {
        invalid <- cols[cols > ncol(expr_data)]
      } else {
        invalid <- setdiff(cols, colnames(expr_data))
      }
      
      if (length(invalid) > 0) {
        stop("Invalid columns specified for time point '", tp_name, 
             "' in group '", group$name, "': ", paste(invalid, collapse = ", "))
      }
      
      # Calculate mean expression
      group_df[[tp_name]] <- apply(expr_data[, cols, drop = FALSE], 1, mean, na.rm = TRUE)
    }
    
    list(name = group$name, data = group_df)
  })
}

# ------------------------------------------------------------------------------
# Perform Mfuzz clustering
# ------------------------------------------------------------------------------
perform_clustering <- function(group_data, params) {
  group_name <- group_data$name
  df <- group_data$data
  cat("Processing", group_name, "data...\n")
  
  # Extract expression values (exclude gene ID column)
  expr_matrix <- data.matrix(df[, -1, drop = FALSE])
  
  # Create ExpressionSet object
  eset <- new("ExpressionSet", exprs = expr_matrix)
  
  # Filter genes with excessive missing values
  filtered_eset <- filter.NA(eset, thres = params$na_threshold)
  
  # Fill remaining missing values
  filled_eset <- fill.NA(filtered_eset, mode = params$na_fill_method)
  
  # Filter genes with low standard deviation
  final_eset <- filter.std(filled_eset, min.std = params$min_std)
  
  # Standardize data
  standardized_eset <- standardise(final_eset)
  
  # Estimate fuzzifier parameter m or use value from parameters
  m_value <- ifelse(is.null(params$fuzzifier_m), 
                   mestimate(standardized_eset), 
                   params$fuzzifier_m)
  
  # Perform clustering
  clusters <- mfuzz(standardized_eset, c = params$num_clusters, m = m_value)
  
  # Generate and save cluster plot
  plot_file <- file.path(params$output_dir, "plots", paste0(group_name, "_clusters.png"))
  png(plot_file, 
      width = params$plot_width, 
      height = params$plot_height, 
      res = params$plot_resolution)
  
  # Calculate plot layout
  rows <- params$plot_rows
  cols <- ceiling(params$num_clusters / rows)
  
  mfuzz.plot(
    standardized_eset, 
    clusters, 
    mfrow = c(rows, cols), 
    new.window = FALSE,
    time.labels = colnames(df[, -1]),
    main = paste(group_name, "Gene Expression Clusters")
  )
  dev.off()
  
  # Save cluster assignments
  cluster_results <- data.frame(
    Gene = df$Gene[match(rownames(exprs(standardized_eset)), rownames(expr_matrix))],
    Cluster = clusters$cluster,
    Membership = apply(clusters$membership, 1, max)
  )
  
  output_file <- file.path(params$output_dir, "clusters", paste0(group_name, "_clusters.txt"))
  write.table(
    cluster_results,
    output_file,
    quote = FALSE,
    row.names = FALSE,
    sep = "\t"
  )
  
  list(
    group = group_name,
    clusters = clusters,
    results = cluster_results
  )
}

# ------------------------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------------------------
main <- function() {
  # Initialize
  load_packages()
  
  # Validate parameters
  validate_parameters(params)
  
  # Setup output directories
  setup_directories(params$output_dir)
  
  # Print run parameters
  cat("Starting Mfuzz clustering analysis...\n")
  cat("Input directory: ", params$input_dir, "\n")
  cat("Output directory: ", params$output_dir, "\n")
  cat("Expression file: ", params$expression_file, "\n")
  cat("Gene ID column: ", params$gene_id_column, "\n")
  cat("Number of clusters: ", params$num_clusters, "\n")
  cat("Analyzing groups: ", paste(sapply(params$sample_groups, function(x) x$name), collapse = ", "), "\n\n")
  
  # Run analysis pipeline
  expr_data <- load_data(params)
  group_datasets <- prepare_group_data(expr_data, params)
  
  # Process each group
  lapply(group_datasets, function(group) {
    perform_clustering(group, params)
  })
  
  cat("\nAnalysis complete! Results saved to", params$output_dir, "\n")
}

# Run main function
main()

