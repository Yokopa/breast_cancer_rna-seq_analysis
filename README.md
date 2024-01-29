# RNA-Seq analysis of breast cancer subtypes (GSE52194 dataset)
Welcome to my Breast Cancer project for the RNA-sequencing course (467713-HS2023) at the University of Bern! 
This project involves the RNA-seq analysis of breast cancer subtypes using the GSE52194 dataset to investigate differentially expressed genes.

To reproduce this project, follow these steps. All the scripts I used are located in the 'scripts' directory of this repository. Scripts from 1 to 13 were run in the IBU cluster, while script 14 was run locally using RStudio.
Please be aware that you need to modify the paths to match your work environment.

### 1. Dataset
The dataset includes 3 replicates each from 3 different subtypes of human breast tumors and 3 healthy samples. The fastq files can be found in the IBU clusterat the path /data/courses/rnaseq_course/breastcancer_de/reads. The samples are a subset from Eswaran et al. 2012 and the fastq files were downloaded through the Gene Expression Omnibus (GEO), accession GSE52194. The library preparation protocol did not preserve information on the transcribed strand (i.e. unstranded). The libraries were sequenced on an Illumina HiSeq 2000 in paired-end mode.

### 2. Quality control
To ensure the raw data's suitability for subsequent analysis steps, perform a quality control using FastQC (version 0.11.9) with the bash script `01_run_fastqc.sh`.
```bash
sbatch 01_run_fastqc.sh
```
Download the HTML output files to your remote computer and check the generated reports for insights into the data quality and potential issues.

### 3. Map reads to the reference genome
Download the latest reference genome sequence and associated annotation for your species from the Ensembl ftp site. For the sequence, use the file named species.assembly.dna.primary_assembly.fa.gz (under DNA(FASTA)). For the annotation, use the gtf (under gene sets) named species.assembly.build.gtf.gz. 
```bash
sbatch 02_download_ref_genome.sh
```
Make sure you check that the downloaded files are intact by computing checksums sum <yourfile> and comparing them to the values in the CHECKSUMS file on the ftp server!

Produce all required index files for read alignment using Hisat2 (version 2.2.1) with the following bash script:
```bash
sbatch 03_hisat2_indexing.sh
```
With Hisat2, map the reads to the reference genome for each sample separately. To do so for all the samples with an array job, first create a tsv file with all samples:
```bash
sbatch 04_create_sample_list.sh
```
Then, run the following code:
```bash
sbatch 05_hisat2_map_reads.sh
```

Convert the resulting sam files to bam format using Samtools (version 1.10) submitting an array job:
```bash
sbatch 06_create_sam_list.sh
sbatch 07_samtools_convert.sh
```
Delete the sam files:
```bash
sbatch 08_cleanup_sam_files.sh
```

Sort the resulting bam files by genomic coordinates using Samtools submitting an array job:
```bash
sbatch 09_create_bam_list.sh
sbatch 10_samtools_sort.sh
```

Index the coordinate sorted bam files using Samtools with an array job:
```bash
sbatch 11_create_sortedbam_list.sh
sbatch 12_index_bam.sh
```
### 4. Count the number of reads per gene
Use all of your bam files as input for featureCounts to produce a table of counts containing the number of reads per gene in each sample. This step will also require the annotation file you downloaded together with the reference sequence. Submit the script number 13:
```bash
sbatch 13_count_reads.sh
```
### 6. Differential expression analysis
Download the counts from featureCounts to your local computer, and run the Differential expression analysis on RStudio using the followign script: `14_run_DESeq2_enrichGO.R `.

Once you complete the differential expression analysis, you can extract the results for the pairwise contrasts of your interest!

### 7. Overrepresentation analysis
You will find also the code for the Overrepresentation analysis in the same script used in the previous step.
You can specify the Gene Ontology (GO) subontology for the Overrepresentation analysis using the `ont` parameter in the `enrichGO` function within the script `14_run_DESeq2_enrichGO.R`. 
The available options for the `ont` parameter are:

- "BP" for Biological Processes
- "MF" for Molecular Functions
- "CC" for Cellular Components
- "ALL" for all three subontologies

For this project, I set `ont=BP` to focus on Biological Processes. However, you have the flexibility to choose other options based on your specific analysis goals.
