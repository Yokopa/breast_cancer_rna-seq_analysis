#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/reads/fastqc_output/fastqc.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/reads/fastqc_output/fastqc.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=10:00:00
#SBATCH --partition=pall

#---------
# FASTQC
# --------

# Load the module
module add UHTS/Quality_control/fastqc/0.11.9

# Specify the directory containing read samples
INPUT_DIR="/data/users/ascarpellini/rnaseq_course/reads"
OUTPUT_DIR="/data/users/ascarpellini/rnaseq_course/reads/fastqc_results"

for file in $INPUT_DIR/*.fastq.gz; do
    fastqc -o $OUTPUT_DIR $file
done
