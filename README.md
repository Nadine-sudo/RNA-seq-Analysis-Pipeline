# RNA-seq Analysis Pipeline  
A complete, modular RNA-seq data analysis pipeline covering **upstream preprocessing**, **differential expression (DE) analysis**, **functional clustering**, and **isoform characterization**. Designed for reproducibility, flexibility, and ease of use—adaptable to diverse species, experimental designs (e.g., time-series, case-control), and research goals.


## 📊 Pipeline Overview  
This pipeline streamlines end-to-end RNA-seq analysis, from raw sequencing data to functional insights. Below is a high-level workflow schematic:  

| Analysis Stage         | Corresponding Directory       | Core Scripts                          | Key Outputs                              |
|------------------------|--------------------------------|---------------------------------------|------------------------------------------|
| 1. Upstream Processing | `upstream_analysis/`           | `RNAseq_pipeline.sh`                  | Cleaned FASTQs, alignment BAMs, StringTie quantification |
| 2. Data Postprocessing | `postprocessing/`              | `RNAseq_postprocessing.sh` | Gene/transcript count matrices (for DE analysis) |
| 3. Differential Expression | `differential_expression/`   | `ballgown_analysis.R` <br> `deseq2_analysis.R` | DE gene/transcript lists, volcano plots, MA plots |
| 4. Functional Clustering | `functional_analysis/`        | `mfuzz_clustering.R`                  | Expression pattern clusters, cluster-specific gene lists |
| 5. Isoform Analysis    | `isoform_analysis/`            | `isoform_analysis_part1.R`/`part2.R`  | Isoform switch lists, domain annotations, coding potential results |


## 📂 Repository Structure  
'''plaintext
RNA-seq-Analysis-Pipeline/
├── README.md                      # Main documentation (you're here)
├── upstream_analysis/             # Step 1: Raw data → Quantification
│   ├── RNAseq_pipeline.sh         # Upstream main script (QC → alignment → quantification)
│   └── rnaseq_pipeline.config.sh  # Upstream parameter config (paths, threads, thresholds)                    
├── postprocessing/                # Step 2: Quantification → DE input
│   ├── RNAseq_postprocessing.sh   # Postprocessing main script (format conversion)
│   ├── rnaseq_postprocessing.config.sh  # Postprocessing parameter config
│   ├── prepDE.py                  # Convert StringTie outputs to count matrices
│   ├── getTPM.py                  # Convert StringTie outputs to TPM matrices
│   └── getFPKM.py                 # Convert StringTie outputs to FPKM matrices
├── differential_expression/       # Step 3: DE analysis (choose one)
│   ├── ballgown_analysis.R        # Transcript-level DE analysis
│   └── deseq2_analysis.R          # Gene-level DE analysis
├── functional_analysis/           # Step 4: Functional clustering (optional)
│   └── mfuzz_clustering.R         # Cluster DE genes by expression patterns
└── isoform_analysis/              # Step 5: Isoform analysis (optional)
├── isoform_analysis_part1.R   # Preprocessing + sequence extraction
└── isoform_analysis_part2.R   # External tool integration + final analysis
'''

## 🚀 Quick Start  
Follow these steps to run the pipeline with your data.


### 1. Prerequisites  
Install required software and packages before starting:  

#### System Tools  
- Fastp (raw data QC)   
- HISAT2 (genome alignment)  
- StringTie (transcript assembly & quantification)  
- Samtools (BAM file processing)  
- Python 2.7+  
- R 4.0+  

**Installation Example**:  
```bash
1. Install R Packages
Required for DE analysis, clustering, and isoform analysis:
# Run this in an R console
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c(
  "ballgown", "DESeq2", "IsoformSwitchAnalyzeR", "clusterProfiler",
  "org.Hs.eg.db", "org.Mm.eg.db", "tidyverse", "ggplot2", "MFuzz",
  "ggrepel", "pheatmap", "corrplot"
))
install.packages(c("devtools", "VennDiagram"))
devtools::install_github("YuLab-SMU/ggtree")

Note: Install species-specific annotation packages (e.g., org.Hs.eg.db for human) if analyzing non-mouse data.
2. Clone the Repository
git clone https://github.com/your-username/RNA-seq-Analysis-Pipeline.git
cd RNA-seq-Analysis-Pipeline

3. Step 1: Upstream Processing (Raw Data → Quantification)
Handles raw data QC, genome alignment, and transcript quantification with StringTie.
Configure Parameters
Edit the upstream config file to match your dataset and environment:
nano upstream_analysis/rnaseq_pipeline.config.sh

Key Parameters to Set:
RAW_DATA_DIR: Path to raw FASTQ files (e.g., ~/data/RNAseq/raw/)
CLEAN_DATA_DIR: Path to save cleaned FASTQs (e.g., upstream_analysis/outputs/clean_fastq/)
GENOME_INDEX: Path to HISAT2 genome index (e.g., ~/genomes/human/hg38/hisat2_index/)
GTF_FILE: Path to genome annotation GTF (e.g., ~/genomes/human/hg38/annotation/Homo_sapiens.GRCh38.110.gtf)
THREADS: Number of CPU threads (e.g., 16; match to your system’s capacity)
Run the Upstream Script
cd upstream_analysis
bash RNAseq_pipeline.sh rnaseq_pipeline.config.sh

Output:
Cleaned FASTQs: outputs/clean_fastq/
Alignment BAMs: outputs/bam/
StringTie quantification: outputs/stringtie/ (TPM/FPKM files + assembled GTFs)
Logs: logs/ (for troubleshooting failed steps)
4. Step 2: Data Postprocessing (Quantification → DE Input)
Converts StringTie outputs to standardized count matrices (gene/transcript level) for DE analysis.
Configure Parameters
Edit the postprocessing config file:
cd ../postprocessing
nano rnaseq_postprocessing.config.sh

Key Parameters to Set:
QUANT_DIR: Path to StringTie outputs (from Step 3: ../upstream_analysis/outputs/stringtie/)
OUTPUT_DIR: Path to save count matrices (e.g., postprocessing/outputs/)
SAMPLE_METADATA: Path to your sample metadata CSV (format: sample_id,group; e.g., sample1,control\nsample2,treatment)
Run the Postprocessing Script
bash RNAseq_postprocessing.sh rnaseq_postprocessing.config.sh

Output:
gene_count_matrix.csv: Gene-level raw count matrix 
transcript_count_matrix.csv: Transcript-level raw count matrix 
normalized_tpm_matrix.csv: Gene-level TPM matrix 
5. Step 3: Differential Expression Analysis (Choose One)
Select the script based on your preference.
Option A: DESeq2
For standard differential expression (gene/transcript):
cd ../differential_expression
Rscript deseq2_analysis.R \
  --count_matrix ../postprocessing/outputs/gene_count_matrix.csv \
  --metadata ../postprocessing/sample_metadata.csv \
  --contrast "treatment,control" \  # Format: "test_group,control_group"
  --p_cutoff 0.05 \
  --fc_cutoff 1 \  # log2(fold-change) threshold
  --output_dir results/

Option B: Ballgown
Rscript ballgown_analysis.R \
  --count_matrix ../upstream_analysis/outputs/stringtie/ \
  --metadata ../postprocessing/sample_metadata.csv \
  --contrast "treatment,control" \
  --p_cutoff 0.05 \
  --fc_cutoff 1 \
  --output_dir results/

Output (for both options):
de_results.csv: Filtered DE gene/transcript list (log2FC, p-value, adjusted p-value)
volcano_plot.png: Volcano plot of DE results
ma_plot.png: MA plot (mean expression vs. log2FC)
de_heatmap.png: Heatmap of top 50 DE genes/transcripts
6. Step 4: Functional Clustering (Optional)
Use mfuzz_clustering.R to cluster DE genes by expression patterns (ideal for time-series or multi-condition data).
cd ../functional_analysis
Rscript mfuzz_clustering.R \
  --de_list ../differential_expression/results/de_results.csv \
  --expr_matrix ../postprocessing/outputs/normalized_tpm_matrix.csv \
  --num_clusters 5 \  # Adjust based on your data
  --output_dir results/

Output:
cluster_results.csv: Gene-to-cluster mapping
mfuzz_clusters.png: Visualization of expression patterns per cluster
cluster_stats.csv: Number of genes and enriched terms per cluster
7. Step 5: Isoform Analysis (Optional)
Characterizes isoform switches and functional changes (coding potential, domain loss/gain).
Part 1: Preprocessing & Sequence Extraction
First, extract isoform sequences and generate external tool guidance:
cd ../isoform_analysis
Rscript isoform_analysis_part1.R \
  --gtf_file ../upstream_analysis/outputs/stringtie/merged.gtf \
  --quant_dir ../upstream_analysis/outputs/stringtie/ \
  --output_dir outputs/

Action Required: After running Part 1, use the generated guides (in outputs/web_tool_guide/) to run these external tools:
CPC2: Coding potential prediction (https://cpc2.gao-lab.org/)
PFAM: Protein domain annotation (https://www.ebi.ac.uk/interpro/)
SignalP: Signal peptide prediction (https://services.healthtech.dtu.dk/services/SignalP-6.0/)
IUPred2A: Intrinsic disorder prediction (https://iupred2a.elte.hu/)
Save all external tool results to isoform_analysis/external_results/ (use exact filenames: cpc2_results.txt, pfam_results.txt, etc.).
Part 2: Integrate Results & Final Analysis
Rscript isoform_analysis_part2.R \
  --input_dir outputs/ \
  --external_dir external_results/ \
  --output_dir final_results/

Output:
isoform_switches.csv: List of isoform switches with functional consequences
coding_potential_results.csv: CPC2-derived coding potential scores
domain_annotation.csv: PFAM domain predictions
isoform_visualization.png: Plots of top isoform switches
⚙️ Configuration Details
All critical parameters are stored in .config.sh files or in the head of Rscripts (no need to edit core scripts). 
📁 Key Output Descriptions
Directory
File Name
Purpose
upstream_analysis/
outputs/stringtie/
StringTie quantification (TPM/FPKM)
postprocessing/
gene_count_matrix.csv
Input for DESeq2
differential_expression/
de_results.csv
Filtered DE genes/transcripts
functional_analysis/
cluster_results.csv
Gene-to-cluster mapping
isoform_analysis/
isoform_switches.csv
Isoform switches with functional impacts

❗ Troubleshooting
Upstream script fails: Check upstream_analysis/logs/ for error messages. Common issues: missing genome index, insufficient threads, or malformed FASTQ files.
DE analysis errors: Ensure count matrices have no missing values. Use the --filter_low_counts flag in RNAseq_postprocessing.sh to remove low-expression genes.
Isoform analysis issues: Verify external tool results are saved to external_results/ with exact filenames (e.g., cpc2_results.txt).
🤝 Contributing
We welcome contributions to improve the pipeline! To contribute:
Fork the repository
Create a feature branch (git checkout -b feature/new-tool)
Commit your changes (git commit -m "Add Salmon quantification support")
Push to the branch (git push origin feature/new-tool)
Open a Pull Request
📞 Support
For questions or issues, open a GitHub Issue or contact the maintainer at [smart_lotus@163.com].





