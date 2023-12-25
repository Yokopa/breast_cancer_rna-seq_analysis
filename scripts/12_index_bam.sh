#!/bin/bash
#SBATCH --job-name=index_bam
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/index_bam.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/index_bam.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=1:00:00
#SBATCH --partition=pall
#SBATCH --array=1-12

#-----------------------------------
# SAMTOOLS (INDEX SORTED BAM FILES)
# ----------------------------------

module add UHTS/Analysis/samtools/1.10

# Path to the TSV file containing sorted BAM file information
BAM_LIST="sorted_BAM_list.tsv"

# Extract information from the TSV file using awk
sample_name=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $BAM_LIST)
sorted_bam_path=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $BAM_LIST)

# Index the sorted BAM file
echo "Indexing sorted BAM file ${sample_name} at $(date)"
samtools index $sorted_bam_path
echo "Finished indexing sorted BAM file ${sample_name} at $(date)"