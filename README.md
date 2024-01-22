# RNA-Seq analysis of breast cancer subtypes (GSE52194 dataset).
Welcome to my Breast Cancer project for the RNA-sequencing course (467713-HS2023) at the University of Bern! 
This project involves the RNA-seq analysis of breast cancer subtypes using the GSE52194 dataset to investigate differentially expressed genes.

To reproduce this project, follow these steps. All the scripts I used are located in the 'scripts' directory of this repository.
Please be aware that you need to modify the paths to match your work environment.

### 1. Dataset

### 2. Quality control
To ensure the raw data's suitability for subsequent analysis steps, perform a quality control using FastQC (version 0.11.9) with the bash script 01_run_fastqc.sh.
```
sbatch 01_run_fastqc.sh
```
Download the HTML output files to your remote computer and check the generated reports for insights into the data quality and potential issues.

### 3. Map reads to the reference genome
description

### 4. Count the number of reads per gene

### 5. Exploratory data analysis

### 6. Differential expression analysis

### 7. Overrepresentation analysis
