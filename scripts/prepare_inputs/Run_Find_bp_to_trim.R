args = commandArgs(trailingOnly=TRUE)

input_vector=args[1]
dir_geneticMap=args[2]
output_file=args[3]


library(FoundHaplo)

Find_bp_to_trim(input_vector,dir_geneticMap,output_file)

