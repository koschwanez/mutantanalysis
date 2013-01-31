#! /bin/bash

# Pipeline for paired end fastq to bam file and variant calls
# $1 is prefix for fastq files (note R1 and R2 labels)
# $2 is read group id for bam header. This can be anything but
# should be flow cell name and flow cell. 
# $GATK is location of GenomeAnalysisTK.jar: www.broadinstitute.org/gatk
# $VARSCAN is location of Varscan: varscan.sourceforge.net 
# bwa and samtools need to be in environmental path
# genome needs to be indexed for both bwa and samtools (see man pages)

bwa aln ref_seq/s288c_sgd.fa \
    $1.R1.fastq > $1.R1.sai

bwa aln ref_seq/s288c_sgd.fa \
    $1.R2.fastq > $1.R2.sai
    
bwa sampe \
	-r '@RG\tID:'$2'\tSM:'$1'\tPL:ILLUMINA' \
	ref_seq/s288c_sgd.fa \
	$1.R1.sai $1.R2.sai \
	$1.R1.fastq $1.R2.fastq \
	> $1.unsorted.sam
	
samtools view -Sbt \
	ref_seq/s288c_sgd.fa \
	$1.unsorted.sam > $1.unsorted.bam
	
samtools sort $1.unsorted.bam $1.sorted

samtools index $1.sorted.bam

java -Xmx2g -jar $GATK \
   -I $1.sorted.bam \
   -R ref_seq/s288c_sgd.fa \
   -T RealignerTargetCreator \
   -o $1.intervals



java -Xmx2g -jar $GATK \
	-I $1.sorted.bam \
	-R ref_seq/s288c_sgd.fa \
	-T IndelRealigner \
	-targetIntervals $1.intervals \
	-o $1.bam

samtools index $1.bam

samtools mpileup \
	-f /Users/johnkoschwanez/sequencing/sam_ref/s288c_sgd.fa \
	$1.bam > $1.pileup

java -jar $VARSCAN \
    readcounts $1.pileup --min-coverage 0 --output-file $1.readcounts
