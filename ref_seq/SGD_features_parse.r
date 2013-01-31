# SGD types
# "ARS"                                 "ARS consensus sequence"             
# "CDEI"                                "CDEII"                              
# "CDEIII"                              "CDS"                                
# "ORF"                                 "W_region"                           
# "X_element_combinatorial_repeats"     "X_element_core_sequence"            
# "X_region"                            "Y'_element"                         
# "Y_region"                            "Z1_region"                          
# "Z2_region"                           "binding_site"                       
# "centromere"                          "external_transcribed_spacer_region" 
# "five_prime_UTR_intron"               "gene_cassette"                      
# "insertion"                           "internal_transcribed_spacer_region" 
# "intron"                              "long_terminal_repeat"               
# "mating_locus"                        "multigene locus"                    
# "ncRNA"                               "non_transcribed_region"             
# "noncoding_exon"                      "not in systematic sequence of S288C"
# "not physically mapped"               "plus_1_translational_frameshift"    
# "pseudogene"                          "rRNA"                               
# "repeat_region"                       "retrotransposon"                    
# "snRNA"                               "snoRNA"                             
# "tRNA"                                "telomere"                           
# "telomeric_repeat"                    "transposable_element_gene"


library(plyr)
PROMOTER_LENGTH <- 750
TERMINATOR_LENGTH <- 250

convert_to_promoter <- function(input_df) {
	start <- input_df$start_position
	if (input_df$strand[1] == 'W') {
		input_df$start_position <- pmax(1, start - PROMOTER_LENGTH)
		input_df$stop_position <- pmax(1, start - 1)
	}
	if (input_df$strand[1] == 'C') {
	 	input_df$start_position <- start + PROMOTER_LENGTH
		input_df$stop_position <- start + 1
	}
	input_df$type <- "promoter"
	input_df$SGDID_pt <- paste(input_df$SGDID, 'p', sep='_')
	return(input_df)
}

convert_to_terminator <- function(input_df) {
	stop <- input_df$stop_position
	if (input_df$strand[1] == 'W') {
		input_df$start_position <- stop + 1
		input_df$stop_position <- stop + TERMINATOR_LENGTH
	}
	if (input_df$strand[1] == 'C') {
	 	input_df$start_position <- stop - 1
		input_df$stop_position <- stop - TERMINATOR_LENGTH
	}
	input_df$type <- "terminator"
	input_df$SGDID_pt <- paste(input_df$SGDID, 't', sep='_')
	return(input_df)
}

sgd_df <- read.delim('/Users/johnkoschwanez/sequencing/S288C_reference_genome_r64/SGD_features.tab', 
	col.names=c('SGDID', 'type', 'qualifier', 'feature_name',
		'gene_name', 'alias', 'parent_feature_name',
		'secondary_SGDID', 'chromosome', 'start_position',
		'stop_position', 'strand', 'genetic_position',
		'coordinate_version', 'sequence_version',
		'description'))
sgd_df$SGDID_pt <- sgd_df$SGDID
orf_df <- subset(sgd_df, (type == "ORF"))

promoter_df <- ddply(orf_df,
	.(strand),
	convert_to_promoter)

terminator_df <- ddply(orf_df,
	.(strand),
	convert_to_terminator)

sgd_df <- rbind(sgd_df, promoter_df, terminator_df)

sgd_df$min_position <- pmin(sgd_df$start_position, sgd_df$stop_position)
sgd_df$max_position <- pmax(sgd_df$start_position, sgd_df$stop_position)
sgd_df$chromosome <- factor(sgd_df$chromosome)

min_df <- subset(sgd_df, (type == "ORF" |
							type == "promoter" |
							type == "terminator"  |
							type == "ARS" |
							type == "centromere" |
							type == "mating_locus" |
							type == "ncRNA" |
							type == "snRNA" |
							type == "tRNA" | 
							type == "binding site" |
							type == "rRNA" |
							type == "snoRNA" | 
							type == "CDEI" |
							type == "CDEII" |
							type == "CDEIII" |
							type == "intron"
							),
					select=c('chromosome', 'min_position', 'max_position',
						'type', 'gene_name', 'feature_name',
						'strand', 'SGDID', 'SGDID_pt', 'parent_feature_name'))
write.table(min_df, 'SGD_features_min.txt', sep='\t')
write.table(min_df, 'SGD_features_python.txt', quote=F, col.names=F,
 	row.names=F, sep='\t')