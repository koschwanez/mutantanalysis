#! /bin/bash

# This script prepares the reference sequence for sequence analysis using
# bwa and samtools.
# The only argument is the fasta reference sequence.
# bwa and samtools should be set in the environment.

usage()
{
    echo "Usage: $PROGRAM [-h help] reference.fasta"
}

usage_and_exit()
{
    usage
    exit $1
}

PROGRAM=`basename $0`

command -v samtools >/dev/null 2>&1 || { echo >&2 "samtools must be in environmental path" >&2; exit 1; }
command -v bwa >/dev/null 2>&1 || { echo >&2 "bwa must be in environmental path" >&2; exit 1; }

while getopts :h opt
do
    case $opt in
        h) usage_and_exit
            ;;
        '?') echo "$0: invalid option -$optarg" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if ! [ -f "$1" ] 
then
    echo "reference fasta file required"
    exit 1
fi

refseq_ext=$(echo $1 |awk -F . '{if (NF>1) {print $NF}}')
if ! [ "X$refseq_ext" = "Xfa" ] && ! [ "X$refseq_ext" = "Xfasta" ]
then
    echo "reference file must have fa or fasta extension"
    exit 1
fi

bwa index $1

samtools faidx $1
