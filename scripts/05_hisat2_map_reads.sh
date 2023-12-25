#!/bin/bash
#SBATCH --job-name=hisat2_mapping
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_output/mapping_hisat2.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_output/mapping_hisat2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=13:00:00
#SBATCH --partition=pall
#SBATCH --array=1-12

#---------------------
# HISAT 2 (MAPPING)
# --------------------

module add UHTS/Aligner/hisat/2.2.1

# Path to the TSV file
SAMPLE_LIST="sample_list.tsv"
HISAT2_INDEX="/data/users/ascarpellini/rnaseq_course/mapping/hisat2_index/hisat2_human_genome_index"

# Extract sample information from the TSV file using awk
sample_name=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $1; exit}' $SAMPLE_LIST)
read1_path=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $2; exit}' $SAMPLE_LIST)
read2_path=$(awk -v line=$SLURM_ARRAY_TASK_ID 'NR==line{print $3; exit}' $SAMPLE_LIST)

# Output SAM file
OUTPUT_SAM="/data/users/ascarpellini/rnaseq_course/mapping/hisat2_sam_files/${sample_name}_mapped.sam"

# Print sample information for debugging
# echo "Sample Name: $sample_name"
# echo "Read 1 Path: $read1_path"
# echo "Read 2 Path: $read2_path"

# hisat2 mapping command
hisat2 -x $HISAT2_INDEX -1 $read1_path -2 $read2_path -S $OUTPUT_SAM