args = commandArgs(trailingOnly=TRUE)

input_vcf=args[1]
output_DIR=args[2]
sample_info_file=args[3]

library(FoundHaplo)

Phasing_by_pedigree(input_vcf,output_DIR,sample_info_file)

