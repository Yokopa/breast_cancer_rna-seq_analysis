#!/bin/bash

#SBATCH --job-name=sort_bam
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/sort_bam.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/sort_bam.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --time=16:00:00
#SBATCH --partition=pall
#SBATCH --array=1-12

#----------------------------------
# SAMTOOLS (SORT BAM BY COORDINATES)
# ---------------------------------

module add UHTS/Analysis/samtools/1.10

# Path to the TSV file
BAM_LIST="BAM_list.tsv"

# Extract sample information from the TSV file using awk
sample_name=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $BAM_LIST)
bam_path=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $BAM_LIST)

# Define output path
OUTPUT_SORTED_BAM="/data/users/ascarpellini/rnaseq_course/mapping/samtools_sorted_bam/${sample_name}_sorted.bam"

# Define the path of the directory where temporary files will be stored during sorting
TEMP="/data/users/ascarpellini/rnaseq_course/mapping/temp_sorting"

# Sort the bam files by genomic coordinates
echo "Sorting BAM file ${sample_name} at $(date)"
samtools sort -@ 4 -o $OUTPUT_SORTED_BAM -T $TEMP $bam_path