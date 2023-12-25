#!/bin/bash
#SBATCH --job-name=generate_tsv
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/reads/generate_tsv.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/reads/generate_tsv.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --partition=pall

#---------------------
# SAMPLE LIST
# --------------------

# Specify the output TSV file
TSV_FILE="sample_list.tsv"
DATA_DIR="/data/users/ascarpellini/rnaseq_course/reads"

## Write the header to the TSV file
## echo -e "Sample_Name\tRead1_Path\tRead2_Path" > $TSV_FILE

# Loop through the data directory and write the TSV file
for file in $DATA_DIR/*.fastq.gz; do
     # Extract sample name from the filename
    sample_name=$(basename "$file" | cut -d '_' -f 1)

    # Determine if it's R1 or R2
    if [[ $file == *"_R1.fastq.gz" ]]; then
        read1_path=$file
    elif [[ $file == *"_R2.fastq.gz" ]]; then
        read2_path=$file
    
        # Append the sample information to the TSV file
        echo -e "$sample_name\t$read1_path\t$read2_path" >> $TSV_FILE
    fi
done