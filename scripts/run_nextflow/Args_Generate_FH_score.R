
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
gen_allele_mismatch_rate=args[9]
MA_cutoff=args[10]
meiosis=args[11]

imputation_quality_score_cutoff_test=args[12]
frequency_type=args[13]
geneticMap_DIR=args[14]
disease_files_DIR=args[15]
test_file=args[16]
test_name=args[17]
test_list=args[18]
data_type=args[19]
controls_file_DIR=args[20]
save_report_DIR=args[21]
TEMP_DIR=args[22]




class(imputation_quality_score_cutoff_test)="numeric"
class(minor_allele_cutoff)="numeric"
###check below all the time

library(FoundHaplo)

Generate_FH_score(source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff,gen_allele_mismatch_rate,MA_cutoff,meiosis,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_nam,test_list,data_type,controls_file_DIR,save_report_DIR,TEMP_DIR)
