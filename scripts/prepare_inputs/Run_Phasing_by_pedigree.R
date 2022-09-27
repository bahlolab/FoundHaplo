args = commandArgs(trailingOnly=TRUE)

input_vcf=args[1]
dir_output=args[2]
sample_info_file=args[3]

library(FoundHaplo)

Phasing_by_pedigree(input_vcf,dir_output,sample_info_file)

