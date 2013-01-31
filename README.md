# mutantanalysis
by John Koschwanez

This is code I wrote to analyze the results of my experimental evolution. The segregant analysis part of the code accepts three sequences: the ancestor, the clone, and the backcrossed pool. The clone and ancestor are compared to find mutations, and then the segregating percentage of the backcrossed pool is found for each mutation to call a segregant. The pipeline I use for analyzing sequences in in the `paired_fastq_to_bam.sh` shell script and is summarized as follows:

1. Align the sequence to the reference genome using bwa:
    http://bio-bwa.sourceforge.net
2. Generate an indexed bam file using samtools:
    http://samtools.sourceforge.net
3. Realign indels using GATK:
    http://www.broadinstitute.org/gatk
4. Create a pileup using samtools.
5. Call variants using varscan:
    http://varscan.sourceforge.net
    This step is not in the shell script because it differs for each strain. For the ancestor and clone, use the `somatic` command with the ancestor as the "normal" and the clone as the "tumor". This gives a nice, tabbed variant list for both strains. For the segregant pool, use the `snp` and `indel` commands separately. 
6. Run the `segregant_analysis_with_clone` command in `mutantanalysis.py`. I run this from a separate python script - see `example_run.py`. This will generate the `compare.html` output.

The reference sequence and the annotation files I use are in the `ref_seq` directory. Please see an example output in `output_dir`. 

The code can also be run with no segregants and multiple clones using the `compare_multiple_clones` command in `mutantanalysis.py`.

Please help and contribute! Thanks. I am working on a small, sample input for testing and I have a huge TODO list.
