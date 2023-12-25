#!/bin/bash
#SBATCH --job-name=sam_to_bam
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/sam_to_bam.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/sam_to_bam.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=pall
#SBATCH --array=1-12

#----------------------------------
# SAMTOOLS (CONVERT SAM TO BAM)
# ---------------------------------

module add UHTS/Analysis/samtools/1.10

# Path to the TSV file
SAM_LIST="SAM_list.tsv"

# Extract sample information from the TSV file using awk
sample_name=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAM_LIST)
sam_path=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAM_LIST)

# Define output path
OUTPUT_BAM="/data/users/ascarpellini/rnaseq_course/mapping/samtools_bam_files/${sample_name}_mapped.bam"

# Convert SAM to BAM
echo "Starting SAM to BAM conversion for $sample_name at $(date)"
samtools view -hbS $sam_path > $OUTPUT_BAM
echo "Finished SAM to BAM conversion for $sample_name at $(date)"