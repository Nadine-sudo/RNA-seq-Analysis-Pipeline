#!/usr/bin/env bash

usage() {
    NAME=$(basename "$0")
    cat <<EOF
Usage:
  ${NAME} [output_dir]
Post-processing pipeline for RNA-seq analysis (StringTie merge + expression quantification).

In order to configure the pipeline options (tool paths, genome indices etc.)
please copy and edit a file rnaseq_postprocessing.config.sh which must be
placed in the current working directory where this script is launched.

Output files and subdirectories will be created in the
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

# Convert OUTDIR to absolute path
OUTDIR_ABS=$(realpath "$OUTDIR")

# Load configuration variables
if [[ ! -f ./rnaseq_postprocessing.config.sh ]]; then
    usage
    echo "Error: configuration file (rnaseq_postprocessing.config.sh) missing!"
    exit 1
fi
source ./rnaseq_postprocessing.config.sh

# Check required programs
errprog=""
if [[ ! -x "$STRINGTIE" ]]; then
    errprog="stringtie"
fi
if [[ ! -x "$GFFCOMPARE" ]]; then
    errprog="gffcompare"
fi
if ! command -v "$PYTHON2" &> /dev/null; then
    errprog="python2"
fi

if [[ -n "$errprog" ]]; then
    echo "ERROR: $errprog program not found, please edit the configuration script."
    exit 1
fi

# Check required parameters from config
required_params=("INPUT_DIR" "GENOME" "num" "refgtf" "PREPDE_SCRIPT" "GETFPKM_SCRIPT" "GETTPM_SCRIPT")
for param in "${required_params[@]}"; do
    if [[ -z "${!param}" ]]; then
        echo "ERROR: Missing required parameter in configuration file: $param"
        exit 1
    fi
done

# Validate input directory
if [[ ! -d "$INPUT_DIR" ]]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Validate reference GTF
if [[ ! -f "$refgtf" ]]; then
    echo "ERROR: GTF file not found: $refgtf"
    exit 1
fi

# Validate Python scripts
if [[ ! -f "$PREPDE_SCRIPT" ]]; then
    echo "ERROR: prepDE.py script not found: $PREPDE_SCRIPT"
    exit 1
fi
if [[ ! -f "$GETFPKM_SCRIPT" ]]; then
    echo "ERROR: getFPKM.py script not found: $GETFPKM_SCRIPT"
    exit 1
fi
if [[ ! -f "$GETTPM_SCRIPT" ]]; then
    echo "ERROR: getTPM.py script not found: $GETTPM_SCRIPT"
    exit 1
fi

# Create output directories
set -e
mkdir -p "$OUTDIR_ABS"
mkdir -p "$OUTDIR_ABS/logs" "$OUTDIR_ABS/merged_gtf"

# Set directory paths
LOGFILE="$OUTDIR_ABS/logs"
MERGED_DIR="$OUTDIR_ABS/merged_gtf"

# Main pipeline function
pipeline() {
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] #> START: $(basename "$0") $*"
    echo "[$(date +"%Y-%m-%d %H:%M:%S")] Using input directory: $INPUT_DIR"

    # Step 1: StringTie merge
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Starting StringTie merge"
    ls "$INPUT_DIR"/*/stringties/*gtf > "$MERGED_DIR/mergelist.txt"
    
    if [[ ! -s "$MERGED_DIR/mergelist.txt" ]]; then
        echo "ERROR: mergelist.txt is empty - no GTF files found in input directory"
        exit 1
    fi

    "$STRINGTIE" --merge -p "$num" \
        -G "$refgtf" \
        -o "$MERGED_DIR/stringtie_merged.gtf" \
        "$MERGED_DIR/mergelist.txt" \
        > "$LOGFILE/stringtie_merge.log" 2>&1

    # Step 2: GFFCompare
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Running GFFCompare"
    "$GFFCOMPARE" -r "$refgtf" \
        -o "$MERGED_DIR/stringtieCmp" \
        "$MERGED_DIR/stringtie_merged.gtf" \
        > "$LOGFILE/gffcompare.log" 2>&1

    # Step 3: StringTie with merged GTF
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Running StringTie with merged annotation"
    for sample in $(ls -d "$INPUT_DIR"/*/ | xargs -n1 basename); do
        echo "[$(date +"%Y-%m-%d %H:%M:%S")]     - Processing sample: $sample"
        sample_dir="$INPUT_DIR/$sample"
        output_dir="$sample_dir/stringtiesExp"
        
        mkdir -p "$output_dir"
        
        "$STRINGTIE" -p "$num" \
            -G "$MERGED_DIR/stringtie_merged.gtf" \
            -o "$output_dir/$sample.stringtie.exp.gtf" \
            -e -B -A "$output_dir/$sample.stringtie.exp.abund.txt" \
            "$sample_dir/bams/$sample.sorted.bam" \
            > "$LOGFILE/stringtie.log" 2>&1
    done

    # Step 4: Create sample list
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Creating sample list"
    > "$OUTDIR_ABS/sample_lst.txt"  # Ensure file is empty before writing
    
    for sample in $(ls -d "$INPUT_DIR"/*/ | xargs -n1 basename); do
        gtf_path="$INPUT_DIR/$sample/stringtiesExp/$sample.stringtie.exp.gtf"
        if [[ -f "$gtf_path" ]]; then
            echo "$sample $gtf_path" >> "$OUTDIR_ABS/sample_lst.txt"
        else
            echo "WARNING: GTF file not found for $sample: $gtf_path"
        fi
    done

    # Step 5: Generate expression tables
    echo "[$(date +"%Y-%m-%d %H:%M:%S")]   * Generating expression tables"
    "$PYTHON2" "$PREPDE_SCRIPT" -i "$OUTDIR_ABS/sample_lst.txt" \
        > "$LOGFILE/prepDE.log" 2>&1
        
    "$PYTHON2" "$GETFPKM_SCRIPT" -i "$OUTDIR_ABS/sample_lst.txt" \
        > "$LOGFILE/getFPKM.log" 2>&1
        
    "$PYTHON2" "$GETTPM_SCRIPT" -i "$OUTDIR_ABS/sample_lst.txt" \
        > "$LOGFILE/getTPM.log" 2>&1

    echo "[$(date +"%Y-%m-%d %H:%M:%S")] #> DONE. Results in $OUTDIR_ABS"
}

# Run pipeline and log output
pipeline "$@" 2>&1 | tee "$LOGFILE/postprocessing.log"
