args = commandArgs(trailingOnly=TRUE)

input_vcf=args[1]
dir_output=args[2]
type=args[3]
sample_names=args[4]

library(FoundHaplo)

Phasing_by_pedigree(input_vcf,dir_output,type,sample_names)

