#! usr/bin/Python

from mutantanalysis import *

segregant_analysis_with_clone('Title of analysis', 
		'Ancestor name', 
		'input_dir/ancestor.indexed.bam',

		'ref_seq/s288c_sgd.fa', 
		'ref_seq/s288c_chr_names.txt',
		'ref_seq/s288c_ref_annot.txt',
		'output_dir',

		'Clone name', 
		'input_dir/clone.indexed.bam', 
		'input_dir/clone.snp',  # Varscan output
		'input_dir/clone.indel',  # Varscan output
		75,  # Read percent to call the clone a mutation

		'Segregated pool name',
		'input_dir/pool.indexed.bam', 
		'input_dir/pool.snp',  # Varscan output
		'input_dir/pool.indel',  # Varscan output
		90  # Read percent to call segregant
		)

