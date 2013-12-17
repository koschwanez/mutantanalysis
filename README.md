## mutantanalysis
by John Koschwanez

This is code I wrote to analyze the results of my experimental evolution. The segregant analysis part of the code accepts three sequences: the ancestor, the clone, and the backcrossed pool. The clone and ancestor are compared to find mutations, and then the segregating percentage of the backcrossed pool is found for each mutation to call a segregant.

Please see the `compare.html` file in the `example_output` directory for an example output for read reads. I have included a script to make simulated reads in the `generate_test_seq` directory. The reference sequence and annotation that I used are in the `ref_seq` directory.

Also included in this package is a script for generating simulated sequences from a backcrossed pool, and scripts for generating a samtools pileup from a fastq file. 

#### Required programs

- [bwa](http://bio-bwa.sourceforge.net/ "bwa") Used to align reads to the reference genome.
- [samtools](http://samtools.sourceforge.net/ "samtools") Used to generate sorted alignment files and pileups.
- [GATK](http://www.broadinstitute.org/gatk/ "Genome Analysis Toolkit") Used to realign indels.
- [varscan](http://varscan.sourceforge.net/ "varscan") Used to call variants from the reference and between samples.
- [python](http://www.python.org/ "python") The mutant anlaysis code is written in python.
- [Biopython](http://biopython.org/wiki/Biopython/ "Biopython") Python library used to interface with DNA sequences.
- [pysam](http://code.google.com/p/pysam/ "pysam") Python wrapper for samtools.
- [ClustalW](http://www.clustal.org/clustal2/ "clustalW") Command line ClustalW binary.

- [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) Sequence simulation tools. This is only necessary to generate simulated backcrossed sequences.

### Usage

`mutantanalysis.py` assumes that you have a sorted bam file for the ancestor and clone strains and the backcrossed pool, and Varscan variant files for the clone (compared to the ancestor) and the pool (compared to the reference). The script `segtools acps` will generate the required files and run `mutantanalysis.py`. `segtools` also contains two other functions: `singlef2p` generates a samtools pileup from a single-end fastq file and `pairedf2p` generates a samtools pileup from a paired-end fastq file.

`singlef2p` will generate a sorted and indel-realigned bam file and a pileup file from a single-end read sample.

    segtools singlef2p -a sample_file.fq -c sample_name -g path_to_GATK 
        [-h help] -r reference.fasta -v path_to_Varscan
        OPTIONS
        -a sample_file.fq 
        -c sample_name 
        -g path to GATK (optionally set $GATK instead)
        -h help
        -r reference.fasta If this has not been indexed by bwa, segtools will index it using the bwa index routine.
        -v path to Varscan (optionally set $VARSCAN instead)

`pairedf2p` will generate a sorted and indel-realigned bam file and a pileup file from a paired-end read sample.

    segtools pairedf2p -a sample_file_R1.fq -b sample_file_R2.fq 
        -c sample_name -g path_to_GATK
        [-h help] -r reference.fasta -v path_to_Varscan
        OPTIONS
        -a sample_file_R1.fq 
        -b sample_file_R2.fq
        -c sample_name
        -g path to GATK (optionally set $GATK instead)
        -h help
        -r reference.fasta If this has not been indexed by bwa, segtools will index it using the bwa index routine.
        -v path to Varscan (optionally set $VARSCAN instead)

`acps` will generate sorted and indel-realigned bam files and pileup files from three sets of fastq files: ancestor, clone, and pool. The fastq files can be single-end or paired end. An html output will be generated listing the variants, their classification, and their segregation percentage.

    segtools acps -a ancestor_file_R1.fq [-b ancestor_file_R2.fq] [-c ancestor_sample_name] 
        -d clone_file_R1.fq [-e clone_file_R2.fq] [-f clone_sample_name] 
        -g path_to_GATK [-h help] 
        -i pool_file_R1.fq [-j pool_file_R2.fq] [-k pool_sample_name] 
        [-o output_directory] -r reference.fasta -s reference_features.tab 
        [-t analysis_title] -v path_to_Varscan
        OPTIONS
        -a ancestor_file_R1.fq 
        -b ancestor_file_R2.fq Required if paired-end read.
        -c ancestor_name (optional)
        -d clone_file_R1.fq 
        -e clone_file_R2.fq Required if paired-end read.
        -f clone_name (optional)
        -g path to GATK (optionally set $GATK instead)
        -h help
        -i pool_file_R1.fq 
        -j pool_file_R2.fq Required if paired-end read.
        -k pool_name (optional)
        -o output_directory Set to seg_output by default
        -r reference.fasta If this has not been indexed by bwa,
            segtools will index it using the bwa index routine.
        -s reference_features Tabbed file - see example for yeast.
        -t analysis_title Set to "pool" by default.
        -v path to Varscan (optionally set $VARSCAN instead)

`mutantanalysis.py` takes the three sets of variants and sorted bam files and generates an html output to view mutation classification and segregation.

    python mutantanalysis.py [-h] [--title TITLE] --output_dir OUTPUT_DIR
                         --ref_seq REF_SEQ --ref_feat REF_FEAT
                         [--ref_chrom_num REF_CHROM_NUM] --ancestor_bam
                         ANCESTOR_BAM [--ancestor_name ANCESTOR_NAME]
                         --clone_bam CLONE_BAM --clone_snp CLONE_SNP
                         --clone_indel CLONE_INDEL [--clone_name CLONE_NAME]
                         [--add_clone_bam ADD_CLONE_BAM]
                         [--add_clone_name ADD_CLONE_NAME]
                         [--add_clone_snp ADD_CLONE_SNP]
                         [--add_clone_indel ADD_CLONE_INDEL]
                         [--pool_bam POOL_BAM] [--pool_name POOL_NAME]
                         [--pool_snp POOL_SNP] [--pool_indel POOL_INDEL]
                         [--mut_percent MUT_PERCENT]
                         [--seg_percent SEG_PERCENT]


        OPTIONS:
        -h, --help 
        --title Title of analysis
        --output_dir Output directory
        --ref_seq fasta reference sequence
        --ref_feat Reference features tab file
        --ref_chrom_num Tab-file of chrom_name to be used and chrom_id in ref fasta
        --ancestor_bam Ancestor sorted bam file
        --ancestor_name Name of ancestor
        --clone_bam Clone sorted bam file
        --clone_snp List of ancestor to clone snps generated by Varscan
        --clone_indel List of ancestor to clone indels generated by Varscan
        --clone_name Name of first clone
        --add_clone_bam Additional clone: Sorted bam file
        --add_clone_name Additional clone: Name
        --add_clone_snp Additional clone: snp file
        --add_clone_indel Additional clone: indel file
        --pool_bam Pool sorted bam file. Required for segregation analysis.
        --pool_name Name of pool
        --pool_snp Pool snp file. Required for segregation analysis.
        --pool_indel Pool indel file. Required for segregation analysis.
        --mut_percent Percentage to call mutant from ancestor to clone
        --seg_percent Percentage to call segregant

### Tutorial

The best way to describe the pipeline is to take you through an example. If you have your own set of reads in fastq format (ancestor.fastq, clone.fastq, pool.fastq), then use those. Single reads will be in one file, and paired reads will usually come in a matched set typically labeled with `R1` and `R2`. If your files are gzipped (have a `.tar.gz` at the end), then extract them by typing:

    cd directory\where\your\samples\are
    tar xvfz yoursample.fastq.tar.gz

You may occasionally get multiple files for the same sample from the sequencing facility. If you are sure these are from the same end (i.e. not a matched set of paired end reads), then concatenate them as follows:

    cat yoursample1.fastq yoursample2.fastq yoursample3.fastq > yourtotalsample.fastq

If you don't have a set of samples to use, I have written a script to generate simulated reads for an ancestor, clone, and backcrossed pool: `simseg.sh`. This requires a reference genome in fasta format, and a tab delimited file with substitution descriptions and segregation percentages (0, 20, 40, 60, 80, or 100% segregation.) These files are all included in the `generate_test_seq` directory. Modify `mutation file.txt` as necessary to make your own mutations. To generate sequences, run as follows (assuming you are in the `mutantanalysis` main directory):

    cd generate_test_seq
    ./simseg.sh s288c_sgd.fa mutation_file.txt

This script uses the program [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) to generate 4 files: single-end reads for the ancestor and clone (`ancestor.fq` and `clone.fq`), and paired end reads for the backcrossed pool (`pool_r1.fq` and `pool_r2.fq`). Move these files into a working directory. For example, to make a working directory off the main directory (assuming you are still in the `generate_test_seq` directory):

    mkdir ../working_directory
    mv *.fq ../working_directory
    cd ../working_directory

Now you are ready to run the mutantanalysis. You will use the program segtools (Segregation Tools), and I'm assuming that you have installed all the programs listed above. clustalw, samtools, and bwa should be set in your environmental path. This is a nice explanation of how to do this: [PATH questions](http://superuser.com/questions/284342/what-are-path-and-other-environment-variables-and-how-can-i-set-or-use-them/ "Setting the environmental path"). You will also need a reference fasta, an annotated features file, and (optional) a translation of chromosome numbers that are in the reference fasta. I have included these files for yeast in the directory `ref_seq` as the files `s288c_sgd.fa`, `s288c_ref_annot.txt`, and `s288c_sgd.fa.chrom`. I generated the features file from `SGD_features.tab`, which I downloaded from [SGD](http://www.yeastgenome.org/ "SGD") using `SGD_features_parse.r`. Modify this file as necessary for your strain or organism.

segtools contains three programs:
- `singlef2p` (single-end read fastq to pileup) uses a reference sequence to convert a single end fastq file to a samtools pileup file using bwa for the alignment.
- `pairedf2p` (paired-end read fastq to pileup) does the same except starting with paired-end fastq files.
- `acps` (ancestor clone pool segregation) starts with three sets of fastq reads and generates an html output of the variants and the segregation percentage.

To run acps from the working directory we just created, type the following:

    ../segtools acps \
        -a ancestor.fq \
        -c my_ancestor_name \
        -d clone.fq \
        -f my_clone_name \
        -g /your/path/to/GATK
        -i pool_r1.fq \
        -j pool_r2.fq \
        -k my_pool_name \
        -r ../ref_seq/s288c_sgd.fa \
        -s ../ref_seq/s288c_ref_annot.txt \
        -t My_Analysis_Title
        -v /your/path/to/Varscan
