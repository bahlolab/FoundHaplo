
#' take input from the manifest.txt to run the nextflow pipeline

#' @export
#' @examples
##
args = commandArgs(trailingOnly=TRUE)


DCV=args[1]
minor_allele_cutoff=args[2]
imputation_quality_score_cutoff_test=args[3]
frequency_type=args[4]
geneticMap_DIR=args[5]
disease_files_DIR=args[6]
test_file=args[7]
test_name=args[8]
test_list=args[9]
data_type=args[10]
controls_file_DIR=args[11]
save_report_DIR=args[12]
TEMP_DIR=args[13]


class(imputation_quality_score_cutoff_test)="numeric"
class(minor_allele_cutoff)="numeric"

###check below all the time

library(FoundHaplo)

Generate_FH_score(DCV,minor_allele_cutoff=0,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,test_list,data_type,controls_file_DIR,save_report_DIR,TEMP_DIR)

