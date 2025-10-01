#!/usr/bin/env bash
# RNA-seq Pipeline Configuration File
# All parameters should be set here before running the pipeline script
# ------------------------------------------------------------------------------

# 1. Paths to required software tools
#    (Update these to match your system's installation paths)
FASTP="/path/to/fastp"                  # Path to fastp (quality control tool)
HISAT2="/path/to/hisat2"                # Path to HISAT2 (alignment tool)
SAMTOOLS="/path/to/samtools"            # Path to SAMtools (BAM processing tool)
STRINGTIE="/path/to/stringtie"          # Path to StringTie (transcript assembly tool)

# 2. Input sample information
INPUT1="/path/to/sample_R1.fastq.gz"    # Full path to read 1 FASTQ file (gzipped or uncompressed)
INPUT2="/path/to/sample_R2.fastq.gz"    # Full path to read 2 FASTQ file (gzipped or uncompressed)
SampleNAME="SampleID"                   # Unique sample name (used in output filenames)
num=8                                   # Number of threads to use (adjust based on system resources)

# 3. Reference genome information
#    (Update these to match your reference genome indices)
refGenome="/path/to/hisat2/index/genome"  # Prefix for HISAT2 genome index files (e.g., if index files are hg38.1.ht2, use "hg38")
refgtf="/path/to/annotation/genome.gtf"   # Full path to GTF annotation file for the reference genome
