#!/usr/bin/env Rscript

# Isoform Switch Analysis Pipeline - Part 1
# Function: Data preprocessing, isoform identification, and sequence extraction
# Output: Sequence files and intermediate results for downstream analysis

# --------------------------
# Configuration parameters
# --------------------------
config <- list(
  # Input/output paths
  output_dir = "./isoform_results",
  expression_dir = "path/to/expression_data/",
  gtf_file = "path/to/stringtie_merged.gtf",
  
  # Analysis parameters
  read_length = 150,
  gene_expr_cutoff = 10,
  isoform_expr_cutoff = 3,
  genome_package = "BSgenome.Mmusculus.UCSC.mm10",
  org_db = "org.Mm.eg.db"
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
# Create directories
# --------------------------
setup_directories <- function(output_dir) {
  dirs <- c(
    output_dir,
    file.path(output_dir, "plots"),
    file.path(output_dir, "results"),
    file.path(output_dir, "sequences"),
    file.path(output_dir, "external_analysis")  # For web tool results
  )
  sapply(dirs, function(d) dir.create(d, showWarnings = FALSE, recursive = TRUE))
  return(output_dir)
}

# --------------------------
# Main analysis pipeline - Part 1
# --------------------------
main <- function() {
  # Initialization
  load_packages()
  output_dir <- setup_directories(config$output_dir)
  sequence_dir <- file.path(output_dir, "sequences")
  
  # Step 1: Import expression data
  cat("Step 1/6: Importing expression data...\n")
  stringtieQuant <- importIsoformExpression(
    parentDir = config$expression_dir,
    addIsofomIdAsColumn = FALSE,
    readLength = config$read_length
  )
  
  # Filter isoform IDs (keep only those starting with M)
  abundance <- stringtieQuant[[1]]
  counts <- stringtieQuant[[2]]
  keep <- grepl("^M", rownames(abundance))
  stringtieQuant[[1]] <- abundance[keep, ]
  stringtieQuant[[2]] <- counts[keep, ]
  save(stringtieQuant, file = file.path(output_dir, "stringtieQuant.RData"))
  
  # Step 2: Create analysis object
  cat("Step 2/6: Creating switchAnalyzeRlist...\n")
  myDesign <- read.delim("brain/brain.csv", sep = ",", header = TRUE)
  aSwitchList <- importRdata(
    isoformRepExpression = stringtieQuant$abundance,
    isoformCountMatrix   = stringtieQuant$counts,
    designMatrix         = myDesign,
    isoformExonAnnoation = config$gtf_file
  )
  save(aSwitchList, file = file.path(output_dir, "aSwitchList.RData"))
  
  # Step 3: Filter low expression
  cat("Step 3/6: Filtering low expression isoforms...\n")
  aSwitchListFiltered <- preFilter(
    switchAnalyzeRlist = aSwitchList,
    geneExpressionCutoff = config$gene_expr_cutoff,
    isoformExpressionCutoff = config$isoform_expr_cutoff,
    removeSingleIsoformGenes = TRUE
  )
  save(aSwitchListFiltered, file = file.path(output_dir, "aSwitchListFiltered.RData"))
  
  # Step 4: Identify isoform switches
  cat("Step 4/6: Identifying isoform switch events...\n")
  aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
    switchAnalyzeRlist = aSwitchListFiltered,
    reduceToSwitchingGenes = TRUE
  )
  save(aSwitchListAnalyzed, file = file.path(output_dir, "aSwitchListAnalyzed.RData"))
  
  # Step 5: ORF analysis
  cat("Step 5/6: Analyzing open reading frames...\n")
  aSwitchListAnalyzed <- addORFfromGTF(aSwitchListAnalyzed, pathToGTF = config$gtf_file)
  genome_obj <- get(config$genome_package)
  aSwitchListAnalyzed <- analyzeNovelIsoformORF(
    aSwitchListAnalyzed,
    analysisAllIsoformsWithoutORF = TRUE,
    genomeObject = genome_obj
  )
  
  # Step 6: Extract sequences
  cat("Step 6/6: Extracting nucleotide and amino acid sequences...\n")
  aSwitchListAnalyzed <- extractSequence(
    aSwitchListAnalyzed,
    pathToOutput = sequence_dir,
    writeToFile = TRUE,
    genomeObject = genome_obj,
    alpha = 0,
    dIFcutoff = 0
  )
  save(aSwitchListAnalyzed, file = file.path(output_dir, "aSwitchListAnalyzed.RData"))
  
  # ----------------------------------------------------------------------------
  # EXTERNAL ANALYSIS INSTRUCTIONS (save these results for Part 2)
  # ----------------------------------------------------------------------------
  cat("\nPart 1 analysis completed!\n")
  cat("Before running Part 2, perform these external analyses:\n\n")
  
  # CPC2 - Coding potential analysis
  cat("# 1. Coding Potential Analysis with CPC2\n")
  cat("# Website: https://cpc2.gao-lab.org/\n")
  cat("# Input file: ", file.path(sequence_dir, "isoform_switch_analyzer_nt_sequences.fasta\n"))
  cat("# Save result as: ", file.path(output_dir, "external_analysis", "cpc2_results.txt\n\n"))
  
  # PFAM - Domain analysis
  cat("# 2. Protein Domain Analysis with PFAM\n")
  cat("# Website: https://www.ebi.ac.uk/interpro/\n")
  cat("# Input file: ", file.path(sequence_dir, "isoform_switch_analyzer_aa_sequences.fasta\n"))
  cat("# Save result as: ", file.path(output_dir, "external_analysis", "pfam_results.txt\n\n"))
  
  # SignalP - Signal peptide prediction
  cat("# 3. Signal Peptide Prediction with SignalP\n")
  cat("# Website: https://services.healthtech.dtu.dk/services/SignalP-6.0/\n")
  cat("# Input file: ", file.path(sequence_dir, "isoform_switch_analyzer_aa_sequences.fasta\n"))
  cat("# Parameters: Organism = 'Eukarya, Output format = Short output'\n")
  cat("# Save result as: ", file.path(output_dir, "external_analysis", "signalp_results.txt\n\n"))
  
  # IUPred2A - Intrinsic disorder prediction
  cat("# 4. Intrinsic Disorder Prediction with IUPred2A\n")
  cat("# Website: https://iupred2a.elte.hu/\n")
  cat("# Input file: ", file.path(sequence_dir, "isoform_switch_analyzer_aa_sequences.fasta\n"))
  cat("# Save result as: ", file.path(output_dir, "external_analysis", "iupred_results.txt\n\n"))
  
  cat("After completing these analyses, run Part 2: Rscript isoform_analysis_part2.R\n")
}

# Run main function
main()


  
