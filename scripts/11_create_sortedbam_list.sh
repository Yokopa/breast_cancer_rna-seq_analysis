#!/bin/bash
#SBATCH --job-name=generate_sorted_BAM_tsv
#SBATCH --output=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/generate_tsv.out
#SBATCH --error=/data/users/ascarpellini/rnaseq_course/mapping/samtools_output/generate_tsv.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=1:00:00
#SBATCH --partition=pall

#---------------------------------------
# CREATE TSV FILE WITH SORTED BAM FILES
# --------------------------------------

# Specify the output TSV file
TSV_FILE="sorted_BAM_list.tsv"
DATA_DIR="/data/users/ascarpellini/rnaseq_course/mapping/samtools_sorted_bam"

# Check if TSV file already exists
if [ -e "$TSV_FILE" ]; then
    echo "Error: TSV file '$TSV_FILE' already exists. Please choose a different filename."
    exit 1
fi

## Write the header to the TSV file
## echo -e "Sample_Name\sorted_BAM_files" > $TSV_FILE

# Loop through the data directory and write the TSV file
for file in $DATA_DIR/*; do
     # Extract sample name from the filename
    sample_name=$(basename "$file" | cut -d '_' -f 1)

    sorted_BAM_path=$file
    
    # Append the path information to the TSV file
    echo -e "$sample_name\t$sorted_BAM_path" >> $TSV_FILE

done