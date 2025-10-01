#!/usr/bin/env bash

usage() {
    NAME=$(basename "$0")
    cat <<EOF
Usage:
  ${NAME} [output_dir]
Wrapper script for RNA-seq analysis pipeline (fastp + HISAT2 + StringTie).

In order to configure the pipeline options (tool paths, genome indices etc.)
please copy and edit a file rnaseq_pipeline.config.sh which must be
placed in the current working directory where this script is launched.

Output directories "fastp", "bams", "stringties", and "logs" will be created in the
current working directory or, if provided, in the given <output_dir>
(which will be created if it does not exist).

EOF
}

# Set default output directory
OUTDIR="."
if [[ "$1" ]]; then
    if [[ "$1" == "-h" || "$1" == "--help" ]]; then
        usage
        exit 1
    fi
    OUTDIR="$1"
fi

# Convert OUTDIR to absolute path to avoid issues with directory changes
OUTDIR_ABS=$(realpath "$OUTDIR")

# Load configuration variables
if [[ ! -f ./rnaseq_pipeline.config.sh ]]; then
    usage
    echo "Error: configuration file (rnaseq_pipeline.config.sh) missing!"
    exit 1
fi
source ./rnaseq_pipeline.config.sh

# Check required programs
errprog=""
if [[ ! -x "$FASTP" ]]; then
    errprog="fastp"
fi
if [[ ! -x "$HISAT2" ]]; then
    errprog="hisat2"
fi
if [[ ! -x "$SAMTOOLS" ]]; then
    errprog="samtools"
fi
if [[ ! -x "$STRINGTIE" ]]; then
    errprog="stringtie"
fi

if [[ -n "$errprog" ]]; then
    echo "ERROR: $errprog program not found, please edit the configuration script."
    exit 1
fi

# Check required parameters from config
required_params=("INPUT1" "INPUT2" "SampleNAME" "num" "refGenome" "refgtf")
for param in "${required_params[@]}"; do
    if [[ -z "${!param}" ]]; then
        echo "ERROR: Missing required parameter in configuration file: $param"
        exit 1
    fi
done

# Check input files exist
if [[ ! -f "$INPUT1" ]]; then
    echo "ERROR: Input file not found: $INPUT1"
    exit 1
fi
if [[ ! -f "$INPUT2" ]]; then
    echo "ERROR: Input file not found: $INPUT2"
    exit 1
fi

# Validate genome indices and GTF file
if [[ ! -f "${refGenome}.1.ht2" ]]; then
    echo "ERROR: HISAT2 index file not found: ${refGenome}.1.ht2"
    exit 1
fi
if [[ ! -f "$refgtf" ]]; then
    echo "ERROR: GTF file not found: $refgtf"
    exit 1
fi

# Create output directories
set -e
mkdir -p "$OUTDIR_ABS"
mkdir -p "$OUTDIR_ABS/fastp" "$OUTDIR_ABS/bams" "$OUTDIR_ABS/stringties" "$OUTDIR_ABS/logs"

# Set directory paths
FASTPDIR="$OUTDIR_ABS/fastp"
BAMDIR="$OUTDIR_ABS/bams"
STRINGTIEDIR="$OUTDIR_ABS/stringties"
LOGFILE="$OUTDIR_ABS/logs"  

# Main pipeline function
pipeline() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] #> START: $(basename "$0") $*"

    # Get absolute paths of input files
    INPUT1_ABS=$(realpath "$INPUT1")
    INPUT2_ABS=$(realpath "$INPUT2")

    # Step 1: Quality control with fastp
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Processing sample: $SampleNAME"
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Quality control (fastp)"
    "$FASTP" -w "$num" \
        -i "$INPUT1_ABS" -I "$INPUT2_ABS" \
        -o "$FASTPDIR/out.$(basename "$INPUT1")" \
        -O "$FASTPDIR/out.$(basename "$INPUT2")" \
        --html "$FASTPDIR/$SampleNAME.html" \
        --json "$FASTPDIR/$SampleNAME.json" \
        > "$LOGFILE/fastp.log" 2>&1

    # Step 2: Alignment with HISAT2
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Alignment to genome (HISAT2)"
    "$HISAT2" -p "$num" -t --dta \
        -x "$refGenome" \
        -1 "$FASTPDIR/out.$(basename "$INPUT1")" \
        -2 "$FASTPDIR/out.$(basename "$INPUT2")" \
        -S "$BAMDIR/$SampleNAME.sam" \
        > "$LOGFILE/hisat2.log" 2>&1

    # Step 3: SAM to BAM conversion and sorting
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * BAM processing (SAMtools)"
    "$SAMTOOLS" sort -@ "$num" \
        -o "$BAMDIR/$SampleNAME.sorted.bam" \
        "$BAMDIR/$SampleNAME.sam" \
        > "$LOGFILE/samtools.log" 2>&1
    rm "$BAMDIR/$SampleNAME.sam"  # Clean up intermediate file

    # Step 4: Generate flagstat report
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Generating flagstat report"
    "$SAMTOOLS" flagstat -@ "$num" \
        "$BAMDIR/$SampleNAME.sorted.bam" \
        > "$BAMDIR/$SampleNAME.sorted.bam.flagstat"

    # Step 5: Transcript assembly and quantification with StringTie
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Transcript assembly (StringTie)"
    "$STRINGTIE" -p "$num" \
        -G "$refgtf" \
        -o "$STRINGTIEDIR/$SampleNAME.stringtie.gtf" \
        -B -e -A "$STRINGTIEDIR/$SampleNAME.stringtie.abund.txt" \
        "$BAMDIR/$SampleNAME.sorted.bam" \
        > "$LOGFILE/stringtie.log" 2>&1

    echo "[$(date +"%Y-%m-%d %H:%M:%S")] #> DONE. Results in $OUTDIR_ABS"
}

# Run pipeline and log output
pipeline "$@" 2>&1 | tee "$LOGFILE/run.log"
