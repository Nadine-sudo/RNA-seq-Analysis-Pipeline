#!/usr/bin/env Rscript

# Isoform Switch Analysis Pipeline - Part 2
# Function: Integrate external tool results, complete analysis, and generate final reports
# Prerequisite: Web tool analyses completed and results saved to specified locations

# --------------------------
# Configuration parameters (must match Part 1)
# --------------------------
config <- list(
  # Input/output paths
  output_dir = "./isoform_results",
  gtf_file = "path/to/stringtie_merged.gtf",
  
  # Analysis parameters
  genome_package = "BSgenome.Mmusculus.UCSC.mm10",
  org_db = "org.Mm.eg.db",
  gene_mean_cutoff = 0.1,
  consequence_dif_cutoff = 0.1,
  
  # Web tool configuration (matches Part 1)
  web_tools = list(
    cpc2 = list(
      output_file = "cpc2_results.txt"
    ),
    pfam = list(
      output_file = "pfam_results.txt"
    ),
    signalp = list(
      output_file = "signalp_results.txt"
    ),
    iupred2a = list(
      output_file = "iupred_results.txt"
    )
  ),
  
  # Consequence types to analyze
  consequences_of_interest = c(
    'intron_retention', 'coding_potential', 
    'NMD_status', 'domains_identified', 
    'ORF_seq_similarity'
  )
)

# --------------------------
# Load required packages
# --------------------------
load_packages <- function() {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", dependencies = TRUE)
  }
  
  required_pkgs <- c(
    "IsoformSwitchAnalyzeR", "BSgenome", "clusterProfiler", 
    "rtracklayer", "tidyverse", "VennDiagram"
  )
  
  for (pkg in required_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      BiocManager::install(pkg, update = FALSE)
      library(pkg, character.only = TRUE)
    }
  }
  
  library(config$genome_package, character.only = TRUE)
  library(config$org_db, character.only = TRUE)
}

# --------------------------
# Check completeness of external tool results
# --------------------------
check_external_results <- function(web_tools) {
  cat("Checking external tool result files...\n")
  ext_dir <- file.path(config$output_dir, "external_analysis")
  missing <- c()
  
  for (tool in names(web_tools)) {
    tool_info <- web_tools[[tool]]
    output_file <- file.path(ext_dir, tool_info$output_file)
    
    if (!file.exists(output_file)) {
      missing <- c(missing, output_file)
    }
  }
  
  if (length(missing) > 0) {
    cat("\nThe following external tool result files are missing:\n")
    for (file in missing) {
      cat("- ", file, "\n")
    }
    cat("Please complete these tool analyses according to the guide and save results to specified locations\n")
    cat("Guide location: ", file.path(config$output_dir, "web_tool_guide"), "\n")
    return(FALSE)
  } else {
    cat("All external tool result files are ready\n")
    return(TRUE)
  }
}

# --------------------------
# Main analysis pipeline - Part 2
# --------------------------
main <- function() {
  # Initialization
  load_packages()
  output_dir <- config$output_dir
  ext_analysis_dir <- file.path(output_dir, "external_analysis")
  
  # Check if output directory exists
  if (!dir.exists(output_dir)) {
    stop("Output directory does not exist. Please run Part 1 first: ", output_dir)
  }
  
  # Load analysis object from Part 1
  cat("Loading Part 1 analysis results...\n")
  analysis_file <- file.path(output_dir, "aSwitchListAnalyzed.RData")
  if (!file.exists(analysis_file)) {
    stop("Part 1 analysis results not found. Ensure Part 1 completed successfully: ", analysis_file)
  }
  load(analysis_file)  # Loads the aSwitchListAnalyzed object
  
  # Check external tool results
  if (!check_external_results(config$web_tools)) {
    cat("Please complete missing result files and rerun this script\n")
    quit(status = 1)
  }
  
  # Step 1: Integrate external tool results
  cat("\nStep 1/3: Integrating external tool analysis results...\n")
  ext_results <- list()
  for (tool in names(config$web_tools)) {
    ext_results[[tool]] <- file.path(ext_analysis_dir, config$web_tools[[tool]]$output_file)
  }
  
  # Integrate CPC2 results
  if (file.exists(ext_results$cpc2)) {
    aSwitchListAnalyzed <- analyzeCPC2(
      aSwitchListAnalyzed,
      pathToCPC2resultFile = ext_results$cpc2,
      removeNoncodinORFs = TRUE
    )
  }
  
  # Integrate PFAM results
  if (file.exists(ext_results$pfam)) {
    aSwitchListAnalyzed <- analyzePFAM(
      aSwitchListAnalyzed,
      pathToPFAMresultFile = ext_results$pfam,
      showProgress = FALSE
    )
  }
  
  # Integrate SignalP results
  if (file.exists(ext_results$signalp)) {
    aSwitchListAnalyzed <- analyzeSignalP(
      aSwitchListAnalyzed,
      pathToSignalPresultFile = ext_results$signalp
    )
  }
  
  # Integrate IUPred2A results
  if (file.exists(ext_results$iupred2a)) {
    aSwitchListAnalyzed <- analyzeIUPred2A(
      aSwitchListAnalyzed,
      pathToIUPred2AresultFile = ext_results$iupred2a,
      showProgress = FALSE
    )
  }
  
  # Step 2: Analyze alternative splicing and switch consequences
  cat("Step 2/3: Analyzing alternative splicing and switch consequences...\n")
  
  # Analyze alternative splicing
  aSwitchListAnalyzed <- analyzeAlternativeSplicing(
    switchAnalyzeRlist = aSwitchListAnalyzed,
    quiet = TRUE
  )
  
  # Analyze intron retention
  aSwitchListAnalyzed <- analyzeIntronRetention(aSwitchListAnalyzed)
  cat("Intron retention summary:\n")
  print(table(aSwitchListAnalyzed$AlternativeSplicingAnalysis$IR))
  
  # Predict switch consequences
  aSwitchListAnalyzed <- analyzeSwitchConsequences(
    aSwitchListAnalyzed,
    consequencesToAnalyze = config$consequences_of_interest,
    dIFcutoff = config$consequence_dif_cutoff,
    showProgress = FALSE
  )
  
  # Add gene annotations (gene symbols)
  tmp <- as.data.frame(aSwitchListAnalyzed$isoformFeatures)
  gene_id <- gsub('\\..*', '', tmp$gene_id)
  gene.df <- bitr(
    gene_id, 
    fromType = "REFSEQ", 
    toType = "SYMBOL",
    OrgDb = org.Mm.eg.db,
    drop = FALSE
  )
  
  tmp$gene <- gsub('\\..*', '', tmp$gene_id)
  tmp$gene_name <- gene.df[match(tmp$gene, gene.df$REFSEQ), "SYMBOL"]
  aSwitchListAnalyzed$isoformFeatures <- tmp
  
  # Step 3: Generate final results and visualizations
  cat("Step 3/3: Generating final results and visualizations...\n")
  
  # Create subset for visualization
  aSwitchListAnalyzedSubset <- subsetSwitchAnalyzeRlist(
    aSwitchListAnalyzed,
    subset = aSwitchListAnalyzed$isoformFeatures$gene_overall_mean > config$gene_mean_cutoff
  )
  save(aSwitchListAnalyzedSubset, file = file.path(output_dir, "aSwitchListAnalyzedSubset.RData"))
  
  # Extract top switching events
  top_switches_q <- extractTopSwitches(
    aSwitchListAnalyzedSubset,
    filterForConsequences = TRUE,
    n = 10,
    sortByQvals = TRUE
  )
  write.csv(
    top_switches_q,
    file.path(output_dir, "results", "top_switches_by_qvalue.csv"),
    row.names = FALSE
  )
  
  # Generate switch visualization plots
  switchPlotTopSwitches(
    switchAnalyzeRlist = aSwitchListAnalyzed,
    n = 20,
    filterForConsequences = TRUE,
    splitFunctionalConsequences = TRUE,
    pathToOutput = file.path(output_dir, "plots"),
    fileType = "png"
  )
  
  # Save final analysis results
  save(aSwitchListAnalyzed, file = file.path(output_dir, "aSwitchListAnalyzed_final.RData"))
  
  cat("\nPart 2 analysis completed!\n")
  cat("All results saved to: ", output_dir, "\n")
  cat("Key result files include:\n")
  cat("- Final analysis object: aSwitchListAnalyzed_final.RData\n")
  cat("- Top switching events: results/top_switches_by_qvalue.csv\n")
  cat("- Visualization plots: files in the plots/ directory\n")
}

# Run main function
main()
