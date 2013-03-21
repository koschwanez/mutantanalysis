#! /bin/bash

# Pipeline for paired end fastq to pileup
# Required options;
#   -b basename of all generated files
#   -r reference fasta sequence. Must be indexed by bwa and samtools
#   -s sample fastq sequence for first pair (R1)
#   -t sample fastq sequence for second pair (R2)
# Optional options;
#   -g location of GenomeAnalysisTK.jar: www.broadinstitute.org/gatk
#       or can be set as $GATK in environemnt
#   -i read group id. Normally flow cell & lane. basename by default.   
# bwa and samtools need to be in environmental path

usage()
{
    echo "Usage: $PROGRAM -b basename [-g GATK_path] [-h help] [-i readgroup_id] -r reference.fasta -s sample_firstreadpair.fastq -t sample_secondreadpair.fastq"
}

usage_and_exit()
{
    usage
    exit $1
}

PROGRAM=`basename $0`
VERSION=0.1

while getopts :b:g:hi:r:s:t: opt
do
    case $opt in
        b) basefilename=$OPTARG
            ;;
        g) gatk=$OPTARG
            ;;
        h) usage_and_exit
            ;;
        i) rgid=$OPTARG
            ;;
        r) refseq=$OPTARG
            refseq_ext=$(echo $refseq |awk -F . '{if (NF>1) {print $NF}}')
            ;;
        s) sampseqR1=$OPTARG
            sampseqR1_ext=$(echo $sampseqR1 |awk -F . '{if (NF>1) {print $NF}}')
            ;;
        t) sampseqR2=$OPTARG
            sampseqR2_ext=$(echo $sampseqR2 |awk -F . '{if (NF>1) {print $NF}}')
            ;;
        '?') echo "$0: invalid option -$optarg" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools must be in environmental path" >&2; exit 1; }
command -v bwa >/dev/null 2>&1 || { echo >&2 "bwa must be in environmental path" >&2; exit 1; }

if ! [ $basefilename ] 
then
    echo "Specify base filename with option -b"
    exit 1
fi

if ! [ ${gatk:=$GATK} ]
then
    echo "gatk needs to be defined with option -g"
    exit 1
fi

echo ${rgid:=$basefilename} > /dev/null

if ! [ -f "$refseq" ] 
then
    echo "reference fasta file required with option -r"
    exit 1
fi

if ! [ "X$refseq_ext" = "Xfa" ] && ! [ "X$refseq_ext" = "Xfasta" ]
then
    echo "reference file must have fa or fasta extension"
    exit 1
fi

if ! [ -f "$sampseqR1" ] || ! [ -f "$sampseqR2" ]
then
    echo "sample fasta files required with options -s and -t"
    exit 1
fi

if ! [ "X$sampseqR1_ext" = "Xfq" ] && ! [ "X$sampseqR1_ext" = "Xfastq" ]
then
    echo "sample files must have fq or fastq extension"
    exit 1
fi

if ! [ "X$sampseqR2_ext" = "Xfq" ] && ! [ "X$sampseqR2_ext" = "Xfastq" ]
then
    echo "sample files must have fq or fastq extension"
    exit 1
fi

echo "now aligning $sampseqR1"
bwa aln $refseq \
    $sampseqR1 > $sampseqR1.sai

echo "now aligning $sampseqR2"
bwa aln $refseq \
    $sampseqR2 > $sampseqR2.sai
    
echo "now creating $basefilename.unsorted.sam"
bwa sampe \
	-r '@RG\tID:'$rgid'\tSM:'$basefilename'\tPL:ILLUMINA' \
	$refseq \
	$sampseqR1.sai $sampseqR2.sai \
	$sampseqR1 $sampseqR2 \
	> $basefilename.unsorted.sam
	
echo "now creating $basefilename.unsorted.bam"
samtools view -Sbt \
	$refseq \
	$basefilename.unsorted.sam > $basefilename.unsorted.bam
	
echo "now sorting $basefilename.unsorted.bam"
samtools sort $basefilename.unsorted.bam $basefilename.sorted

echo "now indexing $basefilename.sorted.bam"
samtools index $basefilename.sorted.bam

echo "now running GATK RealignerTargetCreator"
java -Xmx2g -jar $gatk \
   -I $basefilename.sorted.bam \
   -R $refseq \
   -T RealignerTargetCreator \
   -o $basefilename.intervals

echo "now running GATK IndelRealigner"
java -Xmx2g -jar $gatk \
	-I $basefilename.sorted.bam \
	-R $refseq \
	-T IndelRealigner \
	-targetIntervals $basefilename.intervals \
	-o $basefilename.bam

echo "now indexing realigned file"
samtools index $basefilename.bam

echo "Now creating pileup"
samtools mpileup \
	-f $refseq \
	$basefilename.bam > $basefilename.pileup

