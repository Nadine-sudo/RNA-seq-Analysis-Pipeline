#!/usr/bin/env bash
# RNA-seq Post-processing Pipeline Configuration File
# All parameters should be set here before running the postprocessing script
# ------------------------------------------------------------------------------

# 1. Paths to required software tools
#    (Update these to match your system's installation paths)
STRINGTIE="/path/to/stringtie"          # Path to StringTie (transcript assembly tool)
GFFCOMPARE="/path/to/gffcompare"        # Path to gffcompare (annotation comparison tool)
PYTHON2="/path/to/python2"              # Path to Python 2 executable

# 2. Analysis parameters
INPUT_DIR="/path/to/rnaseq_results"     # Directory containing sample subdirectories from main pipeline
num=8                                   # Number of threads to use

# 3. Reference genome information
refgtf="/path/to/data/reference_gtf"  # Path to reference GTF file

# 4. Paths to Python scripts
PREPDE_SCRIPT="/path/to/prepDE.py"    # Path to prepDE.py script
GETFPKM_SCRIPT="/path/to/getFPKM.py"  # Path to getFPKM.py script
GETTPM_SCRIPT="/path/to/getTPM.py"    # Path to getTPM.py script
