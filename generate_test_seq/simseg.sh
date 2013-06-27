#! /bin/bash

# $1 is reference fasta
# $2 is tab-delimited file of mutations as follows:
#   Col 1: original sequence
#   Col 2: substituted sequence
#   Col 3: gene name
#   Col 4: reference mutation only ("ref") or mutation from ancestor ("anc")
#   Col 5: segregation percentage (0, 20, 40, 60, 80, or 100)
#   Col 6: mutation type: "SNP", "deletion", "insertion"
#   Col 7: coding change: "nonsyn_ORF", "promoter", or "syn_ORF"

# Make one line fasta files from reference (makes sed substitution easier)
awk -f fasta2onelineseq.awk $1 > oneline_ref.fasta

# Make temporary sed script to hold next substitution, then run sed
# to make another fasta file. Repeat for all needed sequences.

# First ancestor sequence
awk -f make_mutation_sed.awk \
    r_or_a=^ref$ seg_perc=. $2 > ref.sed
sed -f ref.sed oneline_ref.fasta > ancestor.fasta
art_illumina --noALN --id ancestor_\
    --in ancestor.fasta --len 50 --fcov 20 --out ancestor

# ... seg100 sequence
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^100$ $2 > seg100.sed
sed -f seg100.sed ancestor.fasta > seg100.fasta
art_illumina --noALN --id seg100_ --paired --mflen 300 --sdev 20 \
    --in seg100.fasta --len 100 --fcov 10 --out seg100_r

# ... seg80 sequence
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^80$ $2 > seg80.sed
sed -f seg80.sed seg100.fasta > seg80.fasta
art_illumina --noALN --id seg80_ --paired --mflen 300 --sdev 20\
    --in seg80.fasta --len 100 --fcov 10 --out seg80_r

# ... seg60 sequence
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^60$ $2 > seg60.sed
sed -f seg60.sed seg80.fasta > seg60.fasta
art_illumina --noALN --id seg60_ --paired --mflen 300 --sdev 20\
    --in seg60.fasta --len 100 --fcov 10 --out seg60_r

# ... seg40 sequence
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^40$ $2 > seg40.sed
sed -f seg40.sed seg60.fasta > seg40.fasta
art_illumina --noALN --id seg40_ --paired --mflen 300 --sdev 20\
    --in seg40.fasta --len 100 --fcov 10 --out seg40_r

# ... seg20 sequence
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^20$ $2 > seg20.sed
sed -f seg20.sed seg40.fasta > seg20.fasta
art_illumina --noALN --id seg20_ --paired --mflen 300 --sdev 20\
    --in seg20.fasta --len 100 --fcov 10 --out seg20_r

# ... seg0 sequence (which is also the clone sequence)
awk -f make_mutation_sed.awk \
    r_or_a=anc seg_perc=^0$ $2 > seg0.sed
sed -f seg0.sed seg20.fasta > seg0.fasta
art_illumina --noALN --id seg0_ --paired --mflen 300 --sdev 20\
    --in seg0.fasta --len 100 --fcov 10 --out seg0_r
art_illumina --noALN --id clone_ \
    --in seg0.fasta --len 50 --fcov 20 --out clone

cat seg*_r1.fq > pool_r1.fq
cat seg*_r2.fq > pool_r2.fq
rm seg0_r*.fq seg20_r*.fq seg40_r*.fq seg60_r*.fq seg80_r*.fq seg100_r*.fq
rm seg0.fasta seg20.fasta seg40.fasta seg60.fasta seg80.fasta seg100.fasta
rm ancestor.fasta
rm oneline_ref.fasta
rm *.sed
