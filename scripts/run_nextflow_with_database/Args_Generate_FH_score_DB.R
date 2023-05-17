
#' take input from the manifest.txt to run the nextflow pipeline

#' @export
#' @examples
##
args = commandArgs(trailingOnly=TRUE)

db_port=args[1]
db_host=args[2]
db_password=args[3]
db_name=args[4]
db_unix_socket=args[5]
DCV=args[6]
minor_allele_cutoff=args[7]
imputation_quality_score_cutoff_test=args[8]
frequency_type=args[9]
geneticMap_DIR=args[10]
disease_files_DIR=args[11]
test_file=args[12]
test_name=args[13]
test_list=args[14]
data_type=args[15]
controls_file_DIR=args[16]
save_report_DIR=args[17]
TEMP_DIR=args[18]


class(imputation_quality_score_cutoff_test)="numeric"
class(minor_allele_cutoff)="numeric"
class(db_port)="numeric
###check below all the time

library(FoundHaplo)

Generate_FH_score_DB(db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff=0,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,test_list,data_type,controls_file_DIR,save_report_DIR,TEMP_DIR)

