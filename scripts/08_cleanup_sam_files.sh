#!/bin/bash
#SBATCH --job-name=cleanup_sam_files
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_sam_files/cleanup_sam_files.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_sam_files/cleanup_sam_files.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G 
#SBATCH --time=0:30:00 

#---------------------
# CLEANUP SAM FILES
# --------------------

# Path of the directory where the sam files are stored
SAM_DIR="/data/users/ascarpellini/rnaseq_course/mapping/hisat2_sam_files"

# Remove all SAM files in the directory
rm $SAM_DIR/*.sam

# Write a txt file explaining that the sam files have been removed
touch $SAM_DIR/README.md 
echo -e "These 12 SAM files below have been removed after their conversion to BAM format.\nYou can find the BAM files in: '/data/users/ascarpellini/rnaseq_course/mapping/samtools_bam_files'.\nHER21_mapped.sam\tNonTNBC1_mapped.sam\tNormal1_mapped.sam\tTNBC1_mapped.sam\nHER22_mapped.sam\tNonTNBC2_mapped.sam\tNormal2_mapped.sam\tTNBC2_mapped.sam\nHER23_mapped.sam\tNonTNBC3_mapped.sam\tNormal3_mapped.sam\tTNBC3_mapped.sam" > $SAM_DIR/README.md
