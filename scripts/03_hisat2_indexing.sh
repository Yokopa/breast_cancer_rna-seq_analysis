#!/usr/bin/env bash
#SBATCH --job-name=indexing_hisat2
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_output/indexing_hisat2.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/hisat2_output/indexing_hisat2.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=3:00:00
#SBATCH --partition=pall

#---------------------
# HISAT 2 (INDEXING)
# --------------------

module add UHTS/Aligner/hisat/2.2.1

# Input reference genome FASTA file
GENOME_FILE="/data/users/ascarpellini/rnaseq_course/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
INDEX_BASENAME="/data/users/ascarpellini/rnaseq_course/mapping/hisat2_human_genome_index"

# Hisat2 indexing command
hisat2-build -p 16 $GENOME_FILE $INDEX_BASENAME

