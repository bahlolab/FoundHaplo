args = commandArgs(trailingOnly=TRUE)

input_vcf=args[1]
output_DIR=args[2]
sample_info_file=args[3]
n.cores=args[4]

library(FoundHaplo)
library(foreach)
library(dplyr)
library(vcfR)
library(stringi)
library(stringr)
library(doParallel)
library(parallel)

Phasing_by_pedigree(input_vcf,output_DIR,sample_info_file,n.cores)

