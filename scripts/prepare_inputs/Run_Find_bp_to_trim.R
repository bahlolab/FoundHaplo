args = commandArgs(trailingOnly=TRUE)

input_vector=args[1]
geneticMap_DIR=args[2]
output_file=args[3]


library(FoundHaplo)

Find_bp_to_trim(input_vector, geneticMap_DIR,output_file)

