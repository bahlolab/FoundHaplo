
#' take input from the manifest.txt to run the nextflow pipeline

#' @export
#' @examples
##
args = commandArgs(trailingOnly=TRUE)

source_of_disease_haplotypes=args[1]
db_port=args[2]
db_host=args[3]
db_password=args[4]
db_name=args[5]
db_unix_socket=args[6]
DCV=args[7]
minor_allele_cutoff=args[8]
imputation_quality_score_cutoff_test=args[9]
frequency_type=args[10]
geneticMap_DIR=args[11]
disease_files_DIR=args[12]
test_file=args[13]
test_name=args[14]
test_list=args[15]
data_type=args[16]
controls_file_DIR=args[17]
save_report_DIR=args[18]
TEMP_DIR=args[19]


class(imputation_quality_score_cutoff_test)="numeric"
class(minor_allele_cutoff)="numeric"
class(db_port)="numeric
###check below all the time

library(FoundHaplo)

Generate_FH_score(source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff=0,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,test_list,data_type,controls_file_DIR,save_report_DIR,TEMP_DIR)

