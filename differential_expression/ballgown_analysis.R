#!/usr/bin/env Rscript
# RNA-seq Data Analysis Pipeline
# Transcriptome analysis using ballgown, including expression analysis, 
# differential expression analysis, and functional enrichment analysis

# --------------------------
# 1. Install and load required packages
# --------------------------
install_packages <- function() {
  # Install Bioconductor manager
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", dependencies = TRUE)
  }
  
  # Define required packages
  bioc_pkgs <- c(
    "ballgown", "MatrixGenerics", "clusterProfiler", 
    "org.Mm.eg.db", "ggtree"
  )
  
  cran_pkgs <- c(
    "devtools", "tidyverse", "pheatmap", 
    "UpSetR", "corrplot", "ggrepel"
  )
  
  # Install Bioconductor packages
  for (pkg in bioc_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg == "ggtree") {
        devtools::install_github('YuLab-SMU/ggtree')
      } else {
        BiocManager::install(pkg, update = FALSE)
      }
      library(pkg, character.only = TRUE)
    }
  }
  
  # Install CRAN packages
  for (pkg in cran_pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE)
      library(pkg, character.only = TRUE)
    }
  }
}

# --------------------------
# 2. Configure analysis parameters
# --------------------------
config <- list(
  base_dir = "/path/to/working/",   # Base working directory
  data_dir = "data",                # Raw data directory (relative to base_dir)
  pheno_file = "phenodata.csv",     # Sample information file
  output_prefix = "results_",       # Prefix for output files
  min_nonzero = 1,                  # Minimum value for log transformation
  fpkm_sum_cutoff = 5,              # Threshold for filtering low expression by FPKM sum
  fc_threshold = 1,                 # Log2 fold change threshold for differential expression
  pval_threshold = 0.05,            # P-value threshold for differential expression
  go_display_number = 10,           # Number of GO terms to display
  organism = "mmu"                  # Species code for KEGG analysis (mouse: mmu)
)

# --------------------------
# 3. Data preparation and loading
# --------------------------
load_data <- function(config) {
  # Create output directory
  if (!dir.exists(file.path(config$base_dir, "output"))) {
    dir.create(file.path(config$base_dir, "output"), recursive = TRUE)
  }
  
  # Load sample information
  pheno_path <- file.path(config$base_dir, config$pheno_file)
  if (!file.exists(pheno_path)) {
    stop("Sample information file not found: ", pheno_path)
  }
  pheno_data <- read.table(pheno_path, header = TRUE, sep = ',')
  
  # Create ballgown object
  cat("Creating ballgown object...\n")
  bg <- ballgown(
    dataDir = file.path(config$base_dir, config$data_dir),
    samplePattern = "",
    pData = pheno_data
  )
  
  list(bg = bg, pheno_data = pheno_data)
}

# --------------------------
# 4. Expression data extraction and quality control
# --------------------------
perform_qc <- function(bg, pheno_data, config) {
  # Extract expression data
  cat("Extracting expression data...\n")
  trans_data <- texpr(bg, 'all')       # Transcript-level data
  gene_expression <- gexpr(bg)         # Gene-level data
  
  # 1. Distribution of transcript counts per gene
  counts <- data.frame(table(trans_data[, "gene_id"]))
  colnames(counts) <- c("gene_id", "numbers")
  
  p1 <- ggplot(data = counts) + 
    geom_histogram(aes(x = numbers), fill = "#66C3A5") + 
    labs(title = "Distribution of transcript count per gene") +
    xlab("Transcripts per gene") +
    theme_bw()
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "transcript_count_distribution.png")),
    plot = p1, width = 8, height = 6
  )
  
  # 2. Distribution of transcript lengths
  p2 <- ggplot() + 
    geom_histogram(aes(x = trans_data$length), fill = "#8DA1CB") + 
    labs(title = "Distribution of transcript lengths") +
    xlab("Transcript length (bp)") +
    theme_bw()
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "transcript_length_distribution.png")),
    plot = p2, width = 8, height = 6
  )
  
  # 3. Visualize transcript structure for specific gene
  if ("NM_001001176.2" %in% trans_data$gene_id) {
    pdf(file.path(config$base_dir, "output", paste0(config$output_prefix, "transcript_structure_NM_001001176.2.pdf")))
    plotTranscripts(
      gene = 'NM_001001176.2', 
      gown = bg, 
      samples = '1m1_l', 
      meas = 'FPKM', 
      colorby = 'transcript', 
      main = 'Transcripts from gene NM_001001176.2 (1m1_l, FPKM)'
    )
    dev.off()
  }
  
  # 4. Sample reproducibility check (scatter plot)
  fpkm <- texpr(bg, meas = "FPKM")
  colnames(fpkm) <- gsub("FPKM.", "", pheno_data$id)
  
  if ("c1_l" %in% colnames(fpkm) && "c2_l" %in% colnames(fpkm)) {
    c1_l <- log2(fpkm[, "c1_l"] + config$min_nonzero)
    c2_l <- log2(fpkm[, "c2_l"] + config$min_nonzero)
    
    p3 <- ggplot() + 
      geom_hex(aes(x = c1_l, y = c2_l)) + 
      geom_abline(intercept = 0, slope = 1, size = 1, color = "red") +
      xlab("c1_l (log2 FPKM)") + 
      ylab("c2_l (log2 FPKM)") + 
      theme_bw()
    
    ggsave(
      file.path(config$base_dir, "output", paste0(config$output_prefix, "sample_correlation_scatter.png")),
      plot = p3, width = 8, height = 6
    )
  }
  
  # 5. Filter low expression transcripts
  fpkm_df <- data.frame(fpkm)
  fpkm_df[,"Row_Sum"] <- rowSums(fpkm_df)
  judge_sum <- fpkm_df[,"Row_Sum"] > config$fpkm_sum_cutoff
  ballgown_filtted <- subset(bg, "judge_sum", genomesubset = TRUE)
  index <- which(fpkm_df[,"Row_Sum"] > config$fpkm_sum_cutoff)
  data_mat <- fpkm_df[index, -ncol(fpkm_df)]
  
  # 6. Transcript sharing across samples (UpSet plot)
  judge_data <- function(x) { ifelse(x == 0, 0, 1) }
  upset_data <- data.frame(apply(data_mat, 2, judge_data))
  col_names <- colnames(upset_data)
  
  pdf(file.path(config$base_dir, "output", paste0(config$output_prefix, "transcript_upset_plot.pdf")))
  upset(upset_data, sets = col_names, order.by = c("freq"), main.bar.color = "#009688")
  dev.off()
  
  # 7. Intersample correlation heatmap
  mat <- cor(data_mat)
  
  pdf(file.path(config$base_dir, "output", paste0(config$output_prefix, "sample_correlation_heatmap.pdf")))
  corrplot(corr = mat, is.corr = FALSE, type = "upper", diag = FALSE)
  dev.off()
  
  # 8. MDS analysis
  dist_mat <- dist(t(data_mat))
  mds_data <- cmdscale(dist_mat)
  
  mds_data <- as.data.frame(scale(mds_data))
  colnames(mds_data) <- c("MDS 1", "MDS 2")
  mds_data$name <- rownames(mds_data)
  
  p4 <- ggplot(data = mds_data) + 
    geom_point(aes(x = `MDS 1`, y = `MDS 2`), size = 3) +
    geom_text(aes(x = `MDS 1` + 0.2, y = `MDS 2`, label = name), size = 2) + 
    theme_bw()
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "mds_plot.png")),
    plot = p4, width = 8, height = 6
  )
  
  # 9. PCA analysis
  pca_fit <- prcomp(t(data_mat))
  components <- summary(pca_fit)$importance[2, c(1, 2)]
  plot_data <- as.data.frame(pca_fit$x)
  plot_data$name <- row.names(plot_data)
  
  p5 <- ggplot(data = plot_data) + 
    geom_point(aes(x = PC1, y = PC2), size = 2) +
    ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = name)) +
    xlab(label = paste0("PC 1 (", round(components[1], 4)*100, "%)")) + 
    ylab(label = paste0("PC 2 (", round(components[2], 4)*100, "%)")) +
    theme_bw()
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "pca_plot.png")),
    plot = p5, width = 8, height = 6
  )
  
  list(
    ballgown_filtted = ballgown_filtted,
    data_mat = data_mat,
    fpkm = fpkm
  )
}

# --------------------------
# 5. Differential expression analysis
# --------------------------
perform_de_analysis <- function(ballgown_filtted, config) {
  cat("Performing differential expression analysis...\n")
  
  # Transcript-level differential analysis
  results_transcripts <- stattest(
    ballgown_filtted,
    feature = "transcript",
    covariate = "group",
    getFC = TRUE,
    meas = "FPKM"
  )
  
  # Merge gene IDs
  results_transcripts <- data.frame(
    REFSEQ = geneIDs(ballgown_filtted),
    results_transcripts
  )
  
  # Label differentially expressed transcripts
  results_transcripts$threshold <- ifelse(
    results_transcripts$pval < config$pval_threshold & 
    abs(log2(results_transcripts$fc)) > config$fc_threshold,
    ifelse(log2(results_transcripts$fc) > config$fc_threshold, 'up', 'down'),
    'no'
  )
  
  # Create volcano plot
  p_volcano <- ggplot(results_transcripts, aes(log2(fc), -log10(pval))) +
    geom_point(aes(color = threshold)) +
    scale_color_manual(values = c("#303F9F", "#757575", "#FF5252")) +
    theme_bw() +
    xlab(expression(log[2] (fold-change))) + 
    ylab(expression(-log[10] (p-value))) +
    ggtitle("Volcano plot of differential expression")
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "volcano_plot.png")),
    plot = p_volcano, width = 8, height = 6
  )
  
  # Convert gene IDs to gene symbols
  gene <- results_transcripts$REFSEQ
  tmp <- gsub("\\..*", "", gene)
  gene.df <- bitr(
    tmp, 
    fromType = "REFSEQ", 
    toType = c("SYMBOL", "ENTREZID"),
    OrgDb = org.Mm.eg.db,
    drop = FALSE
  )
  
  # Merge annotation information
  results_transcripts$REFSEQ <- gsub("\\..*", "", results_transcripts$REFSEQ)
  result <- merge(gene.df, results_transcripts, by = "REFSEQ")
  
  # Filter up and down regulated genes
  result_up <- result[which(result$fc >= 2 & result$pval < config$pval_threshold), ]
  result_down <- result[which(result$fc <= 0.5 & result$pval < config$pval_threshold), ]
  result_total <- rbind(result_up, result_down)
  
  # Save differential expression results
  write.csv(
    result_total, 
    file.path(config$base_dir, "output", paste0(config$output_prefix, "result_total.csv")), 
    row.names = FALSE
  )
  write.csv(
    result_up, 
    file.path(config$base_dir, "output", paste0(config$output_prefix, "result_up.csv")), 
    row.names = FALSE
  )
  write.csv(
    result_down, 
    file.path(config$base_dir, "output", paste0(config$output_prefix, "result_down.csv")), 
    row.names = FALSE
  )
  
  list(
    de_results = result_total,
    entrez_ids = result_total$ENTREZID[!is.na(result_total$ENTREZID)]
  )
}

# --------------------------
# 6. Functional enrichment analysis
# --------------------------
perform_enrichment_analysis <- function(entrez_ids, config) {
  if (length(entrez_ids) == 0) {
    warning("No valid ENTREZ IDs obtained, cannot perform enrichment analysis")
    return(NULL)
  }
  
  cat("Performing functional enrichment analysis...\n")
  
  # GO enrichment analysis
  ego_BP <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    minGSSize = 1,
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_CC <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "CC",
    pAdjustMethod = "BH",
    minGSSize = 1,
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  ego_MF <- enrichGO(
    gene = entrez_ids,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = "MF",
    pAdjustMethod = "BH",
    minGSSize = 1,
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05,
    readable = TRUE
  )
  
  # Save GO enrichment results
  write.csv(
    as.data.frame(ego_BP),
    file.path(config$base_dir, "output", paste0(config$output_prefix, "go_bp.csv")),
    row.names = TRUE
  )
  write.csv(
    as.data.frame(ego_CC),
    file.path(config$base_dir, "output", paste0(config$output_prefix, "go_cc.csv")),
    row.names = TRUE
  )
  write.csv(
    as.data.frame(ego_MF),
    file.path(config$base_dir, "output", paste0(config$output_prefix, "go_mf.csv")),
    row.names = TRUE
  )
  
  # Create GO enrichment bar plot
  display_number <- rep(config$go_display_number, 3)
  ego_partial_result_BP <- as.data.frame(ego_BP)[1:display_number[1], ]
  ego_partial_result_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
  ego_partial_result_MF <- as.data.frame(ego_MF)[1:display_number[3], ]
  
  go_enrich_df <- data.frame(
    ID = c(ego_partial_result_BP$ID, ego_partial_result_CC$ID, ego_partial_result_MF$ID),
    Description = c(ego_partial_result_BP$Description, ego_partial_result_CC$Description, ego_partial_result_MF$Description),
    GeneNumber = c(ego_partial_result_BP$Count, ego_partial_result_CC$Count, ego_partial_result_MF$Count),
    type = factor(
      c(rep("biological process", display_number[1]),
        rep("cellular component", display_number[2]),
        rep("molecular function", display_number[3])),
      levels = c("biological process", "cellular component", "molecular function")
    )
  )
  
  # Simplify pathway names
  for (i in 1:nrow(go_enrich_df)) {
    description_splite <- strsplit(go_enrich_df$Description[i], split = " ")
    description_collapse <- paste(description_splite[[1]][1:5], collapse = " ")
    go_enrich_df$Description[i] <- description_collapse
    go_enrich_df$Description <- gsub(pattern = "NA", "", go_enrich_df$Description)
  }
  
  # Create GO enrichment plot
  go_enrich_df$type_order <- factor(
    rev(as.integer(rownames(go_enrich_df))),
    labels = rev(go_enrich_df$Description)
  )
  COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
  
  p_go <- ggplot(data = go_enrich_df, aes(x = type_order, y = GeneNumber, fill = type)) + 
    geom_bar(stat = "identity", width = 0.8) + 
    scale_fill_manual(values = COLS) + 
    coord_flip() + 
    xlab("GO term") + 
    ylab("Gene Number") + 
    labs(title = "The Most Enriched GO Terms") +
    theme_bw()
  
  ggsave(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "go_enrichment.png")),
    plot = p_go, width = 10, height = 8
  )
  
  # KEGG enrichment analysis (optional, not executed by default)
  if (FALSE) {
    kk <- enrichKEGG(
      gene = entrez_ids,
      keyType = "kegg",
      organism = config$organism,
      qvalueCutoff = 0.05,
      pvalueCutoff = 0.05
    )
    
    if (!is.null(kk)) {
      hh <- as.data.frame(kk)
      write.csv(
        hh,
        file.path(config$base_dir, "output", paste0(config$output_prefix, "kegg_pathways.csv")),
        row.names = TRUE
      )
      
      # Create KEGG enrichment plot
      rownames(hh) <- 1:nrow(hh)
      hh$order <- factor(rev(as.integer(rownames(hh))), labels = rev(hh$Description))
      
      p_kegg <- ggplot(hh, aes(y = order, x = Count)) +
        geom_point(aes(size = Count, color = -1 * p.adjust)) +
        scale_color_gradient(low = "green", high = "red") +
        labs(
          color = expression(p.adjust),
          size = "Count",
          x = "Gene Number",
          y = "Pathways",
          title = "KEGG Pathway Enrichment"
        ) +
        theme_bw()
      
      ggsave(
        file.path(config$base_dir, "output", paste0(config$output_prefix, "kegg_enrichment.png")),
        plot = p_kegg, width = 10, height = 8
      )
    }
  }
  
  # Save gene IDs for online KEGG analysis
  write.csv(
    data.frame(ENTREZID = entrez_ids),
    file.path(config$base_dir, "output", paste0(config$output_prefix, "gene_ids_for_kegg.csv")),
    row.names = FALSE
  )
  
  cat("KEGG analysis note: You can use the online tool KOBAS (http://kobas.cbi.pku.edu.cn/) for pathway enrichment analysis\n")
}

# --------------------------
# 7. Expression heatmap
# --------------------------
plot_expression_heatmap <- function(de_results, data_mat, config) {
  if (nrow(de_results) == 0) {
    warning("No differentially expressed genes, cannot plot heatmap")
    return(NULL)
  }
  
  cat("Plotting differentially expressed gene heatmap...\n")
  
  # Extract expression levels of differentially expressed genes
  id <- de_results$id
  pre_data <- subset(data_mat, rownames(data_mat) %in% id)
  
  # Create heatmap
  pdf(
    file.path(config$base_dir, "output", paste0(config$output_prefix, "de_heatmap.pdf")),
    width = 10, height = 8
  )
  pheatmap(
    mat = log10(pre_data + 1),
    show_rownames = FALSE,
    main = "Differentially Expressed Genes"
  )
  dev.off()
}

# --------------------------
# Main function: Run complete analysis pipeline
# --------------------------
main <- function() {
  # Record start time
  start_time <- Sys.time()
  cat("RNA-seq data analysis pipeline started at:", as.character(start_time), "\n\n")
  
  # 1. Install and load packages
  cat("Step 1/7: Installing and loading required packages...\n")
  install_packages()
  
  # 2. Load data
  cat("\nStep 2/7: Loading data...\n")
  data_objects <- load_data(config)
  bg <- data_objects$bg
  pheno_data <- data_objects$pheno_data
  
  # 3. Quality control and expression analysis
  cat("\nStep 3/7: Performing quality control and expression analysis...\n")
  qc_objects <- perform_qc(bg, pheno_data, config)
  ballgown_filtted <- qc_objects$ballgown_filtted
  data_mat <- qc_objects$data_mat
  
  # 4. Differential expression analysis
  cat("\nStep 4/7: Performing differential expression analysis...\n")
  de_objects <- perform_de_analysis(ballgown_filtted, config)
  de_results <- de_objects$de_results
  entrez_ids <- de_objects$entrez_ids
  
  # 5. Functional enrichment analysis
  cat("\nStep 5/7: Performing functional enrichment analysis...\n")
  perform_enrichment_analysis(entrez_ids, config)
  
  # 6. Plot expression heatmap
  cat("\nStep 6/7: Plotting differentially expressed gene heatmap...\n")
  plot_expression_heatmap(de_results, data_mat, config)
  
  # 7. Complete analysis
  end_time <- Sys.time()
  cat("\nStep 7/7: Analysis complete!\n")
  cat("Analysis ended at:", as.character(end_time), "\n")
  cat("Total runtime:", as.character(end_time - start_time), "\n")
  cat("All results saved to:", file.path(config$base_dir, "output"), "\n")
}

# Run main function
main()

