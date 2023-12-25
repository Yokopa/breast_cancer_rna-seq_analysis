#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=600M
#SBATCH --time=00:10:00

#----------------------------------
# REFERENCE GENOME AND ANNOTATION
# ---------------------------------

# Download latest reference genome 
wget -P /data/users/ascarpellini/rnaseq_course/reference https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# Download the associated annotation
wget -P /data/users/ascarpellini/rnaseq_course/reference https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz


