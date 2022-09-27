
#' take input from the manifest.txt to run the nextflow pipeline

#' @export
#' @examples


args = commandArgs(trailingOnly=TRUE)


DCV=args[1]
minor_allele_cutoff=args[2]
imputation_quality_score_cutoff_test=args[3]
frequency_type=args[4]
dir_geneticMap=args[5]
dir_disease_files=args[6]
test_file=args[7]
test_name=args[8]
test_list=args[9]
data_type=args[10]
dir_controls_file=args[11]
dir_to_save_report=args[12]


class(imputation_quality_score_cutoff_test)="numeric"
class(minor_allele_cutoff)="numeric"

###check below all the time

library(FoundHaplo)

Generate_FH_score(DCV,minor_allele_cutoff=0,imputation_quality_score_cutoff_test,frequency_type,dir_geneticMap,dir_disease_files,test_file,test_name,test_list,data_type,dir_controls_file,dir_to_save_report)


