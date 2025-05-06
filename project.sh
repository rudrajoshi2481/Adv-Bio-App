#BINF05555 PROJECT CODE

#download index file
wget https://genome-idx.s3.amazonaws.com/bt/hg19.zip
unzip hg19.zip

#Download FastQC, Bowtie2, Samtools, MACS2, Bedtools, Bcftools, etc.
conda install -c bioconda fastqc
conda install bioconda::bowtie2
conda install bioconda::samtools
pip install MACS2
apt-get install bedtools
conda install -c bioconda bcftools

#Download read data
#Download the raw reads for library ENCLB464RQG of K562 POLR2A ChIP-seq data (paired-end)
wget https://www.encodeproject.org/files/ENCFF839LPL/@@download/ENCFF839LPL.fastq.gz
wget https://www.encodeproject.org/files/ENCFF698ICA/@@download/ENCFF698ICA.fastq.gz
#Download the raw reads for library ENCLB359ZTD for the K562 control data (paired-end)
wget https://www.encodeproject.org/files/ENCFF002EFD/@@download/ENCFF002EFD.fastq.gz
wget https://www.encodeproject.org/files/ENCFF002EFA/@@download/ENCFF002EFA.fastq.gz
#Download the raw reads for library ENCLB469RNT for the MCF7 POLR2A ChIP-seq data (single-end)
wget https://www.encodeproject.org/files/ENCFF000SBI/@@download/ENCFF000SBI.fastq.gz
#Download the raw reads for library ENCLB262YBZ for the MCF7 control ChIP-seq data (single-end)
wget https://www.encodeproject.org/files/ENCFF000SAZ/@@download/ENCFF000SAZ.fastq.gz

#Use FastQC to analyze
fastqc ENCFF839LPL.fastq.gz
fastqc ENCFF698ICA.fastq.gz
fastqc ENCFF002EFD.fastq.gz
fastqc ENCFF002EFA.fastq.gz
fastqc ENCFF000SBI.fastq.gz
fastqc ENCFF000SAZ.fastq.gz

start_time=$(date +%s)
echo "Start Time: $(date -d @$start_time)"

#Perform alignments
#alignment for LPL_ICA (paired-end ChIP-seq K562)
bowtie2 -x ./hg19/hg19 -1 ENCFF839LPL.fastq.gz -2 ENCFF698ICA.fastq.gz -S K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam
#alignment for EFD_EFA (paired-end control K562)
bowtie2 -x ./hg19/hg19 -1 ENCFF002EFD.fastq.gz -2 ENCFF002EFA.fastq.gz -S K562_Control_ENCFF002EFD_EFA_Aligned.sam
#alignment for SBI (single-end ChIP-seq MCF7)
bowtie2 -x ./hg19/hg19 -U ENCFF000SBI.fastq.gz -S MCF7_POLR2A_ENCFF000SBI_Aligned.sam
#alignment for SAZ (single-end control MCF7)
bowtie2 -x ./hg19/hg19 -U ENCFF000SAZ.fastq.gz -S MCF7_Control_ENCFF000SAZ_Aligned.sam

#Convert .sam to .bam
samtools view -bS K562_POLR2A_ENCFF839LPL_ICA_Aligned.sam >K562_POLR2A_ENCFF839LPL_ICA_Aligned.bam
samtools view -bS K562_Control_ENCFF002EFD_EFA_Aligned.sam >K562_Control_ENCFF002EFD_EFA_Aligned.bam
samtools view -bS MCF7_POLR2A_ENCFF000SBI_Aligned.sam >MCF7_POLR2A_ENCFF000SBI_Aligned.bam
samtools view -bS MCF7_Control_ENCFF000SAZ_Aligned.sam >MCF7_Control_ENCFF000SAZ_Aligned.bam

#Filter by read quality
samtools view -q 30 -b K562_POLR2A_ENCFF839LPL_ICA_Aligned.bam >K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned.bam
samtools view -q 30 -b K562_Control_ENCFF002EFD_EFA_Aligned.bam >K562_Control_ENCFF002EFD_EFA_UniquelyAligned.bam
samtools view -q 30 -b MCF7_POLR2A_ENCFF000SBI_Aligned.bam >MCF7_POLR2A_ENCFF000SBI_UniquelyAligned.bam
samtools view -q 30 -b MCF7_Control_ENCFF000SAZ_Aligned.bam >MCF7_Control_ENCFF000SAZ_UniquelyAligned.bam

#Sort and index
#K562 Control
samtools sort K562_Control_ENCFF002EFD_EFA_UniquelyAligned.bam -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam
samtools index K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam
#K562 POLR2A
samtools sort K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned.bam -o K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam
samtools index K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam
#MCF7 Control
samtools sort MCF7_Control_ENCFF000SAZ_UniquelyAligned.bam -o MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam
samtools index MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam
#MCF7 POLR2A
samtools sort MCF7_POLR2A_ENCFF000SBI_UniquelyAligned.bam -o MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam
samtools index MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam

#Cufflinks
#analyze individual gene expression for each condition
echo "Analyzing individual gene expression"
cufflinks -G hg19_ncbiRefSeq_annotation.gtf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -o K562_Control_ENCFF002EFD_EFA_GeneExpression
cufflinks -G hg19_ncbiRefSeq_annotation.gtf K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam -o K562_POLR2A_ENCFF839LPL_ICA_GeneExpression
cufflinks -G hg19_ncbiRefSeq_annotation.gtf MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam -o MCF7_Control_ENCFF000SAZ_GeneExpression
cufflinks -G hg19_ncbiRefSeq_annotation.gtf MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam -o MCF7_POLR2A_ENCFF000SBI_GeneExpression

#Cuffdiff
#calculate differential gene expression between conditions
echo "calculating differential gene expression"
cuffdiff -o ./K562_MCF7_Control hg19_ncbiRefSeq_annotation.gtf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam
cuffdiff -o ./K562_MCF7_POLR2A hg19_ncbiRefSeq_annotation.gtf K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam

#Call peaks for POLR2A at significance of p<0.00001 for K562 cell lines and p<0.001 for MCF7 cell lines
macs2 callpeak -t K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam -c K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -f BAM --outdir K562_POLR2A -n K562_POLR2A -p 0.00001
macs2 callpeak -t MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam -c MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam -f BAM --outdir MCF_POLR2A -n MCF7_POLR2A -p 0.001

#Extract peak sequences
bedtools getfasta -fi hg19.fa -bed K562_POLR2A_peaks.narrowPeak -fo K562_POLR2A_peaks.narrowPeak_Sequences.fasta
bedtools getfasta -fi hg19.fa -bed MCF7_POLR2A_peaks.narrowPeak -fo MCF7_POLR2A_peaks.narrowPeak_Sequences.fasta

#generate the variants in VCF format
echo "Generating Variants in VCF format"

#K562 Control
bcftools mpileup -Ov -f ./hg19.fa K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.bam -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.vcf
echo "Now Calling SNPs"
bcftools call -mv -Ov -o K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted_calledSNP.vcf K562_Control_ENCFF002EFD_EFA_UniquelyAligned_Sorted.vcf -P 32

#K562 POLR2A
bcftools mpileup -Ov -f ./hg19.fa K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.bam -o K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.vcf
echo "Now Calling SNPs"
bcftools call -mv -Ov -o K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted_calledSNP.vcf K562_POLR2A_ENCFF839LPL_ICA_UniquelyAligned_Sorted.vcf -P 32

#MCF7 Control
bcftools mpileup -Ov -f ./hg19.fa MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.bam -o MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.vcf
echo "Now Calling SNPs"
bcftools call -mv -Ov -o MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted_calledSNP.vcf MCF7_Control_ENCFF000SAZ_UniquelyAligned_Sorted.vcf -P 32

#MCF7 POLR2A 
bcftools mpileup -Ov -f ./hg19.fa MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.bam -o MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.vcf
echo "Now Calling SNPs"
bcftools call -mv -Ov -o MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted_calledSNP.vcf MCF7_POLR2A_ENCFF000SBI_UniquelyAligned_Sorted.vcf -P 32

end_time=$(date +%s)
echo "Finish Time: $(date -d @$end_time)"

time_elapsed=$((end_time - start_time))
days=$((time_elapsed / 86400))
hours=$(( (time_elapsed % 86400) / 3600 ))
minutes=$(( (time_elapsed % 3600) / 60 ))
seconds=$((time_elapsed % 60))

echo "Time Elapsed: ${days}d ${hours}h ${minutes}m ${seconds}s"