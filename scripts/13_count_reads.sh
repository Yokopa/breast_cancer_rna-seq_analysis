#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/counts/subread_output/featureCounts.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/counts/subread_output/featureCounts.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=16:00:00
#SBATCH --partition=pall

#-------------------------------------
# FEATURECOUNTS (COUNT READS PER GENE)
# ------------------------------------

# Load required module
module load UHTS/Analysis/subread/2.0.1

# Paths
SORTED_BAM="/data/users/ascarpellini/rnaseq_course/mapping/samtools_sorted_bam"
ANNOTATION_FILE="/data/users/ascarpellini/rnaseq_course/reference/Homo_sapiens.GRCh38.110.gtf.gz"
OUTPUT_DIR="/data/users/ascarpellini/rnaseq_course/counts"

# Run FeatureCounts
featureCounts -p -t exon -g gene_id -a $ANNOTATION_FILE -o $OUTPUT_DIR/featurecounts.txt $SORTED_BAM/*.bam