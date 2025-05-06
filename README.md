## Project BINF05555: ChIP-seq Data Analysis Pipeline Documentation

**Date:** May 5, 2025

### 1. Overview

This document describes the bioinformatics pipeline scripted in `#BINF05555 PROJECT CODE`. The pipeline is designed to process Chromatin Immunoprecipitation sequencing (ChIP-seq) data for two cell lines, K562 and MCF7, targeting the POLR2A protein, along with their respective control samples.

The workflow includes the following major steps:
1.  **Setup and Dependency Installation:** Downloading the reference genome index and installing necessary bioinformatics tools.
2.  **Data Acquisition:** Downloading raw ChIP-seq read data (FASTQ files) from the ENCODE project.
3.  **Quality Control:** Assessing the quality of the raw reads using FastQC.
4.  **Read Alignment:** Aligning the reads to the hg19 human reference genome using Bowtie2.
5.  **Post-Alignment Processing:**
    * Converting SAM to BAM format using Samtools.
    * Filtering BAM files for uniquely aligned reads with a minimum quality score using Samtools.
    * Sorting and indexing the filtered BAM files using Samtools.
6.  **Gene Expression Analysis (Cufflinks & Cuffdiff):**
    * Quantifying gene expression for each sample using Cufflinks.
    * Identifying differential gene expression between conditions using Cuffdiff.
7.  **Peak Calling:** Identifying regions of significant protein binding (peaks) using MACS2.
8.  **Peak Sequence Extraction:** Extracting the DNA sequences corresponding to the called peaks using Bedtools.
9.  **Variant Calling:** Identifying genetic variants (SNPs) from the aligned reads using Bcftools.
10. **Timing:** Recording the start, end, and total elapsed time for the pipeline execution.

### 2. Prerequisites

* **Operating System:** A Linux-based system with a command-line interface.
* **Internet Connection:** Required for downloading the genome index, tools, and raw data.
* **Conda Environment:** `conda` should be installed and configured for managing software packages.
* **Permissions:** Sufficient permissions to download files, install software, and write output files in the working directory.
* **Reference Genome Annotation File:** A `hg19_ncbiRefSeq_annotation.gtf` file is required for the Cufflinks and Cuffdiff steps. This file is not downloaded by the script and must be present in the working directory or its path correctly specified.
* **Reference Genome FASTA File:** A `hg19.fa` file (presumably part of the `hg19.zip` download, or a separate uncompressed version) is required for Bedtools `getfasta` and Bcftools `mpileup`. The script uses `./hg19.fa` and `./hg19/hg19` for the Bowtie2 index, implying the FASTA file might be expected in the root or within the `hg19` directory.

### 3. Pipeline Steps

#### 3.1. Setup and Dependency Installation

This section prepares the environment by downloading the reference genome index and installing the required bioinformatics software.

1.  **Download Human Genome (hg19) Index for Bowtie2:**
    * Command: `wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip`
    * Description: Downloads a pre-built Bowtie2 index for the hg19 human genome from an AWS S3 bucket.
    * Output: `hg19.zip`

2.  **Unzip Genome Index:**
    * Command: `unzip hg19.zip`
    * Description: Extracts the contents of `hg19.zip`, creating a directory (likely named `hg19`) containing the Bowtie2 index files.
    * Input: `hg19.zip`
    * Output: A directory containing hg19 Bowtie2 index files (e.g., `hg19.1.bt2`, `hg19.2.bt2`, etc.).

3.  **Install Bioinformatics Tools:**
    * Commands:
        ```bash
        conda install -c bioconda fastqc
        conda install bioconda::bowtie2
        conda install bioconda::samtools
        pip install MACS2
        apt-get install bedtools
        conda install -c bioconda bcftools
        ```
    * Description: Installs various bioinformatics tools using `conda`, `pip`, and `apt-get`.
        * `fastqc`: For raw read quality control.
        * `bowtie2`: For aligning sequencing reads to a reference genome.
        * `samtools`: For manipulating SAM/BAM alignment files.
        * `MACS2`: For ChIP-seq peak calling.
        * `bedtools`: For genomic feature manipulation (here, extracting FASTA sequences from BED files).
        * `bcftools`: For variant calling and VCF/BCF file manipulation.
    * Note: Mixing `conda`, `pip`, and `apt-get` for package management can sometimes lead to environment conflicts. It's generally recommended to use `conda` for as many packages as possible, especially within a specific environment. `apt-get` installs system-wide packages and may require `sudo` privileges.

#### 3.2. Data Acquisition

This section downloads the raw sequencing read data (FASTQ files) for the ChIP-seq experiments from the ENCODE project.

1.  **Download K562 POLR2A ChIP-seq Reads (Paired-end):**
    * Library: ENCLB464RQG
    * Commands:
        ```bash
        wget https://www.encodeproject.org/files/ENCFF839LPL/@@download/ENCFF839LPL.fastq.gz
        wget https://www.encodeproject.org/files/ENCFF698ICA/@@download/ENCFF698ICA.fastq.gz
        ```
    * Output: `ENCFF839LPL.fastq.gz` (Read 1), `ENCFF698ICA.fastq.gz` (Read 2)

2.  **Download K562 Control Reads (Paired-end):**
    * Library: ENCLB359ZTD
    * Commands:
        ```bash
        wget https://www.encodeproject.org/files/ENCFF002EFD/@@download/ENCFF002EFD.fastq.gz
        wget https://www.encodeproject.org/files/ENCFF002EFA/@@download/ENCFF002EFA.fastq.gz
        ```
    * Output: `ENCFF002EFD.fastq.gz` (Read 1), `ENCFF002EFA.fastq.gz` (Read 2)

3.  **Download MCF7 POLR2A ChIP-seq Reads (Single-end):**
    * Library: ENCLB469RNT
    * Command:
        ```bash
        wget https://www.encodeproject.org/files/ENCFF000SBI/@@download/ENCFF000SBI.fastq.gz
        ```
    * Output: `ENCFF000SBI.fastq.gz`

4.  **Download MCF7 Control ChIP-seq Reads (Single-end):**
    * Library: ENCLB262YBZ
    * Command:
        ```bash
        wget https://www.encodeproject.org/files/ENCFF000SAZ/@@download/ENCFF000SAZ.fastq.gz
        ```
    * Output: `ENCFF000SAZ.fastq.gz`

#### 3.3. Quality Control (FastQC)

This section assesses the quality of the downloaded raw FASTQ files.

* Commands:
    ```bash
    fastqc ENCFF839LPL.fastq.gz
    fastqc ENCFF698ICA.fastq.gz
    fastqc ENCFF002EFD.fastq.gz
    fastqc ENCFF002EFA.fastq.gz
    fastqc ENCFF000SBI.fastq.gz
    fastqc ENCFF000SAZ.fastq.gz
    ```
* Description: `fastqc` generates an HTML report summarizing various quality metrics for each FASTQ file, such as per-base sequence quality, GC content, and adapter content.
* Input: `*.fastq.gz` files.
* Output: For each input file, an HTML report (e.g., `ENCFF839LPL_fastqc.html`) and a ZIP file (e.g., `ENCFF839LPL_fastqc.zip`) containing the report data.

#### 3.4. Timing Initialization

* Command:
    ```bash
    start_time=$(date +%s)
    echo "Start Time: $(date -d @$start_time)"
    ```
* Description: Records the current Unix timestamp (seconds since epoch) as the start time of the main processing steps and prints it in a human-readable format.

#### 3.5. Read Alignment (Bowtie2)

This section aligns the sequencing reads to the hg19 reference genome.

1.  **Alignment for K562 POLR2A (Paired-end):**
    * Command: `bowtie2 -x ./hg19/hg19 -1 ENCFF839LPL.fastq.gz -2 ENCFF698ICA.fastq.gz -S K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam`
    * Description:
        * `-x ./hg19/hg19`: Specifies the basename of the Bowtie2 index for hg19 (located in the `hg19` subdirectory).
        * `-1 ENCFF839LPL.fastq.gz`: Specifies the first read file for paired-end alignment.
        * `-2 ENCFF698ICA.fastq.gz`: Specifies the second read file for paired-end alignment.
        * `-S K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam`: Specifies the output file name in SAM format.
    * Input: `hg19` index files, `ENCFF839LPL.fastq.gz`, `ENCFF698ICA.fastq.gz`.
    * Output: `K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam`

2.  **Alignment for K562 Control (Paired-end):**
    * Command: `bowtie2 -x ./hg19/hg19 -1 ENCFF002EFD.fastq.gz -2 ENCFF002EFA.fastq.gz -S K562_Control_ENCFF002EFD_EFA_Aligned.sam`
    * Input: `hg19` index files, `ENCFF002EFD.fastq.gz`, `ENCFF002EFA.fastq.gz`.
    * Output: `K562_Control_ENCFF002EFD_EFA_Aligned.sam`

3.  **Alignment for MCF7 POLR2A (Single-end):**
    * Command: `bowtie2 -x ./hg19/hg19 -U ENCFF000SBI.fastq.gz -S MCF7_POLR2A_ENCFF000SBI_Aligned.sam`
    * Description:
        * `-U ENCFF000SBI.fastq.gz`: Specifies the input file for single-end alignment.
    * Input: `hg19` index files, `ENCFF000SBI.fastq.gz`.
    * Output: `MCF7_POLR2A_ENCFF000SBI_Aligned.sam`

4.  **Alignment for MCF7 Control (Single-end):**
    * Command: `bowtie2 -x ./hg19/hg19 -U ENCFF000SAZ.fastq.gz -S MCF7_Control_ENCFF000SAZ_Aligned.sam`
    * Input: `hg19` index files, `ENCFF000SAZ.fastq.gz`.
    * Output: `MCF7_Control_ENCFF000SAZ_Aligned.sam`

#### 3.6. Post-Alignment Processing (Samtools)

This section converts SAM files to the more compact BAM format, filters for high-quality alignments, and sorts/indexes the BAM files.

1.  **Convert SAM to BAM:**
    * Commands:
        ```bash
        samtools view -bS K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam > K562_POLR2A_ENCFF839LPL_ICA_Aligned.bam
        samtools view -bS K562_Control_ENCFF002EFD_EFA_Aligned.sam > K562_Control_ENCFF002EFD_EFA_Aligned.bam
        samtools view -bS MCF7_POLR2A_ENCFF000SBI_Aligned.sam > MCF7_POLR2A_ENCFF000SBI_Aligned.bam
        samtools view -bS MCF7_Control_ENCFF000SAZ_Aligned.sam > MCF7_Control_ENCFF000SAZ_Aligned.bam
        ```
    * Description:
        * `samtools view`: Tool for viewing and converting SAM/BAM files.
        * `-bS`: Output in BAM format (`-b`), input is SAM (`-S`).
    * Input: `*.sam` files from the alignment step.
    * Output: `*.bam` files (e.g., `K562_POLR2A_ENCFF839LPL_ICA_Aligned.bam`).

2.  **Filter by Read Quality (MAPQ >= 30):**
    * Commands:
        ```bash
        samtools view -q 30 -b K562_POLR2A_ENCFF839LPL_ICA_Aligned.bam > K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned.bam
        samtools view -q 30 -b K562_Control_ENCFF002EFD_EFA_Aligned.bam > K562_Control_ENCFF002EFD_EFA_UniquelyAligned.bam
        samtools view -q 30 -b MCF7_POLR2A_ENCFF000SBI_Aligned.bam > MCF7_POLR2A_ENCFF000SBI_UniquelyAligned.bam
        samtools view -q 30 -b MCF7_Control_ENCFF000SAZ_Aligned.bam > MCF7_Control_ENCFF000SAZ_UniquelyAligned.bam
        ```
    * Description:
        * `-q 30`: Filters out alignments with a mapping quality (MAPQ) less than 30. This helps to retain uniquely mapped reads with high confidence.
        * `-b`: Output in BAM format.
    * Input: `*_Aligned.bam` files.
    * Output: `*_UniquelyAligned.bam` files (e.g., `K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned.bam`).

3.  **Sort and Index BAM Files:**
    * This is done for each of the four conditions (K562 Control, K562 POLR2A, MCF7 Control, MCF7 POLR2A). Example for K562 Control:
        ```bash
        #K562 Control
        samtools sort K562_Control_ENCFF002EFD_EFA_UniquelyAligned.bam -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam
        samtools index K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam
        ```
    * Description:
        * `samtools sort`: Sorts the BAM file by genomic coordinates.
        * `-o`: Specifies the output file name for the sorted BAM.
        * `samtools index`: Creates an index file (`.bai`) for the sorted BAM file, allowing for fast random access.
    * Input: `*_UniquelyAligned.bam` files.
    * Output: `*_UniquelyAligned_Sorted.bam` files and their corresponding `*.bai` index files (e.g., `K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam` and `K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam.bai`).
    * The script repeats these two commands for:
        * K562 POLR2A: `K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned.bam`
        * MCF7 Control: `MCF7_Control_ENCFF000SAZ_UniquelyAligned.bam`
        * MCF7 POLR2A: `MCF7_POLR2A_ENCFF000SBI_UniquelyAligned.bam`

#### 3.7. Gene Expression Analysis (Cufflinks & Cuffdiff)

This section uses the Cufflinks suite to estimate gene expression levels and identify differentially expressed genes.
* **Prerequisite:** Requires a gene annotation file in GTF format (`hg19_ncbiRefSeq_annotation.gtf`). This file is not downloaded by the script and must be available.

1.  **Analyze Individual Gene Expression (Cufflinks):**
    * Command: `echo "Analyzing individual gene expression"` (informational)
    * Commands:
        ```bash
        cufflinks -G hg19_ncbiRefSeq_annotation.gtf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -o K562_Control_ENCFF002EFD_EFA_GeneExpression
        cufflinks -G hg19_ncbiRefSeq_annotation.gtf K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam -o K562_POLR2A_ENCFF839LPL_ICA_GeneExpression
        cufflinks -G hg19_ncbiRefSeq_annotation.gtf MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam -o MCF7_Control_ENCFF000SAZ_GeneExpression
        cufflinks -G hg19_ncbiRefSeq_annotation.gtf MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam -o MCF7_POLR2A_ENCFF000SBI_GeneExpression
        ```
    * Description:
        * `cufflinks`: Assembles transcripts, estimates their abundances, and tests for differential expression and regulation.
        * `-G hg19_ncbiRefSeq_annotation.gtf`: Provides a reference gene annotation. Cufflinks will use this to estimate expression of known genes.
        * `-o <output_directory>`: Specifies the directory where output files will be written.
    * Input: Sorted BAM files, `hg19_ncbiRefSeq_annotation.gtf`.
    * Output: For each sample, an output directory (e.g., `K562_Control_ENCFF002EFD_EFA_GeneExpression`) containing files like `genes.fpkm_tracking`, `isoforms.fpkm_tracking`, `skipped.gtf`, and `transcripts.gtf`.

2.  **Calculate Differential Gene Expression (Cuffdiff):**
    * Command: `echo "calculating differential gene expression"` (informational)
    * Commands:
        ```bash
        cuffdiff -o ./K562_MCF7_Control hg19_ncbiRefSeq_annotation.gtf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam
        cuffdiff -o ./K562_MCF7_POLR2A hg19_ncbiRefSeq_annotation.gtf K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam
        ```
    * Description:
        * `cuffdiff`: Finds significant changes in transcript expression, splicing, and promoter use between two or more samples.
        * `-o <output_directory>`: Specifies the output directory.
        * `hg19_ncbiRefSeq_annotation.gtf`: The reference annotation file.
        * The remaining arguments are the sorted BAM files for the conditions being compared.
    * Input: Sorted BAM files for the conditions to be compared, `hg19_ncbiRefSeq_annotation.gtf`.
    * Output:
        * `./K562_MCF7_Control` directory: Contains differential expression results comparing K562 Control and MCF7 Control.
        * `./K562_MCF7_POLR2A` directory: Contains differential expression results comparing K562 POLR2A and MCF7 POLR2A.
        * These directories will contain various output files detailing differential expression at gene and isoform levels (e.g., `gene_exp.diff`, `isoform_exp.diff`).

#### 3.8. Peak Calling (MACS2)

This section identifies genomic regions where the POLR2A protein is enriched (peaks) compared to the control samples.

1.  **Call Peaks for K562 POLR2A:**
    * Command: `macs2 callpeak -t K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam -c K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -f BAM --outdir K562_POLR2A -n K562_POLR2A -p 0.00001`
    * Description:
        * `macs2 callpeak`: Main function for calling peaks.
        * `-t ...bam`: Specifies the treatment (ChIP) BAM file.
        * `-c ...bam`: Specifies the control BAM file.
        * `-f BAM`: Indicates the input file format is BAM.
        * `--outdir K562_POLR2A`: Specifies the output directory for K562 POLR2A results.
        * `-n K562_POLR2A`: Sets the prefix for output file names.
        * `-p 0.00001`: Sets the p-value cutoff for peak detection (more stringent for K562).
    * Input: Sorted BAM files for K562 POLR2A (treatment) and K562 Control.
    * Output: A directory `K562_POLR2A` containing peak files, including `K562_POLR2A_peaks.narrowPeak`, `K562_POLR2A_peaks.xls`, `K562_POLR2A_summits.bed`.

2.  **Call Peaks for MCF7 POLR2A:**
    * Command: `macs2 callpeak -t MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam -c MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam -f BAM --outdir MCF_POLR2A -n MCF7_POLR2A -p 0.001`
    * Description: Similar to the K562 peak calling, but with a less stringent p-value cutoff (`-p 0.001`). The output directory is `MCF_POLR2A` and the name prefix is `MCF7_POLR2A`.
    * Input: Sorted BAM files for MCF7 POLR2A (treatment) and MCF7 Control.
    * Output: A directory `MCF_POLR2A` containing peak files, including `MCF7_POLR2A_peaks.narrowPeak`, `MCF7_POLR2A_peaks.xls`, `MCF7_POLR2A_summits.bed`.

#### 3.9. Extract Peak Sequences (Bedtools)

This section extracts the DNA sequences corresponding to the identified peak regions.
* **Prerequisite:** Requires a reference genome FASTA file (`hg19.fa`). The script assumes this file is in the current working directory. It should correspond to the hg19 genome build used for alignment.

1.  **Extract K562 Peak Sequences:**
    * Command: `bedtools getfasta -fi hg19.fa -bed K562_POLR2A_peaks.narrowPeak -fo K562_POLR2A_peaks.narrowPeak_Sequences.fasta`
    * Description:
        * `bedtools getfasta`: Extracts sequences from a FASTA file based on coordinates in a BED file.
        * `-fi hg19.fa`: Specifies the input reference FASTA file.
        * `-bed K562_POLR2A_peaks.narrowPeak`: Specifies the BED file containing peak coordinates (output from MACS2, likely needs path adjustment e.g. `K562_POLR2A/K562_POLR2A_peaks.narrowPeak`).
        * `-fo K562_POLR2A_peaks.narrowPeak_Sequences.fasta`: Specifies the output FASTA file for the peak sequences.
    * Input: `hg19.fa`, `K562_POLR2A/K562_POLR2A_peaks.narrowPeak` (assuming MACS2 output structure).
    * Output: `K562_POLR2A_peaks.narrowPeak_Sequences.fasta`

2.  **Extract MCF7 Peak Sequences:**
    * Command: `bedtools getfasta -fi hg19.fa -bed MCF7_POLR2A_peaks.narrowPeak -fo MCF7_POLR2A_peaks.narrowPeak_Sequences.fasta`
    * Description: Similar to K562 peak sequence extraction.
    * Input: `hg19.fa`, `MCF_POLR2A/MCF7_POLR2A_peaks.narrowPeak` (assuming MACS2 output structure).
    * Output: `MCF7_POLR2A_peaks.narrowPeak_Sequences.fasta`

    *Important Note:* The commands for `bedtools getfasta` assume the `*_peaks.narrowPeak` files are in the current working directory. However, MACS2 creates them within the specified `--outdir`. The paths in the `bedtools` commands might need to be adjusted to `K562_POLR2A/K562_POLR2A_peaks.narrowPeak` and `MCF_POLR2A/MCF7_POLR2A_peaks.narrowPeak` respectively for the script to run correctly.

#### 3.10. Variant Calling (Bcftools)

This section identifies genetic variants (SNPs and indels) in each sample compared to the reference genome.
* **Prerequisite:** Requires the reference genome FASTA file (`./hg19.fa`).

1.  **Generate Variants for K562 Control:**
    * Command: `echo "Generating Variants in VCF format"` (informational)
    * `bcftools mpileup`:
        * Command: `bcftools mpileup -Ov -f ./hg19.fa K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.vcf`
        * Description:
            * `bcftools mpileup`: Generates VCF or BCF files by calculating genotype likelihoods from input BAM files.
            * `-Ov`: Output uncompressed VCF.
            * `-f ./hg19.fa`: Specifies the reference FASTA file.
            * Input BAM: `K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam`
            * `-o ...vcf`: Specifies the output VCF file.
        * Input: Sorted BAM file for K562 Control, `hg19.fa`.
        * Output: `K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.vcf`
    * `bcftools call`:
        * Command: `echo "Now Calling SNPs"` (informational)
        * Command: `bcftools call -mv -Ov -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted_calledSNP.vcf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.vcf -P 32`
        * Description:
            * `bcftools call`: Calls SNPs and indels.
            * `-m`: Alternative model for multiallelic and rare-variant calling (deprecated in newer versions, `-A` might be preferred or it might imply `--multiallelic-caller`).
            * `-v`: Output variant sites only.
            * `-Ov`: Output uncompressed VCF.
            * `-o ..._calledSNP.vcf`: Specifies the output VCF file for called SNPs.
            * `-P 32`: Priors. This option is related to older versions or specific prior settings. In current bcftools, `-P` might refer to `--pval-threshold` for `consensus-caller` or might be a typo/legacy option. The exact behavior might depend on the `bcftools` version. Assuming it's related to some prior probability for SNP calling.
        * Input: VCF file from `mpileup`.
        * Output: `K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted_calledSNP.vcf`

2.  **Generate Variants for K562 POLR2A:**
    * Commands: Similar `bcftools mpileup` and `bcftools call` steps are performed for the `K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam` file.
    * Outputs:
        * `K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.vcf`
        * `K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted_calledSNP.vcf`

3.  **Generate Variants for MCF7 Control:**
    * Commands: Similar `bcftools mpileup` and `bcftools call` steps are performed for the `MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam` file.
    * Outputs:
        * `MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.vcf`
        * `MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted_calledSNP.vcf`

4.  **Generate Variants for MCF7 POLR2A:**
    * Commands: Similar `bcftools mpileup` and `bcftools call` steps are performed for the `MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam` file.
    * Outputs:
        * `MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.vcf`
        * `MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted_calledSNP.vcf`

#### 3.11. Timing Finalization

* Commands:
    ```bash
    end_time=$(date +%s)
    echo "Finish Time: $(date -d @$end_time)"

    time_elapsed=$((end_time - start_time))
    days=$((time_elapsed / 86400))
    hours=$(( (time_elapsed % 86400) / 3600 ))
    minutes=$(( (time_elapsed % 3600) / 60 ))
    seconds=$((time_elapsed % 60))

    echo "Time Elapsed: ${days}d ${hours}h ${minutes}m ${seconds}s"
    ```
* Description: Records the end time, prints it, and then calculates and prints the total time elapsed for the pipeline execution in days, hours, minutes, and seconds.

### 4. Summary of Outputs

* **FastQC Reports:** HTML and ZIP files for each raw FASTQ file.
* **Alignment Files:**
    * `*.sam` (intermediate)
    * `*_Aligned.bam` (intermediate)
    * `*_UniquelyAligned.bam` (filtered, intermediate)
    * `*_UniquelyAligned_Sorted.bam` (final BAM for downstream analysis)
    * `*_UniquelyAligned_Sorted.bam.bai` (BAM index files)
* **Gene Expression Files (Cufflinks/Cuffdiff):**
    * Directories for each sample from `cufflinks` (e.g., `K562_Control_ENCFF002EFD_EFA_GeneExpression`) containing FPKM tracking and transcript files.
    * Directories for differential expression from `cuffdiff` (`K562_MCF7_Control`, `K562_MCF7_POLR2A`) containing tables of differentially expressed genes/isoforms.
* **Peak Calling Files (MACS2):**
    * Directories `K562_POLR2A` and `MCF_POLR2A` containing:
        * `*_peaks.narrowPeak`: Peak locations in BED format.
        * `*_peaks.xls`: Peak information in tabular format.
        * `*_summits.bed`: Peak summit locations.
* **Peak Sequence Files (Bedtools):**
    * `K562_POLR2A_peaks.narrowPeak_Sequences.fasta`
    * `MCF7_POLR2A_peaks.narrowPeak_Sequences.fasta`
* **Variant Calling Files (Bcftools):**
    * `*_UniquelyAligned_Sorted.vcf` (raw VCF from mpileup)
    * `*_UniquelyAligned_Sorted_calledSNP.vcf` (VCF with called variants)
* **Log/Timing Information:** Printed to standard output.

