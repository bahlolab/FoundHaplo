#' Create manifest.txt to run batch jobs
#'
#' @description
#' Takes all the input parameters required to run FoundHaplo and create manifest.txt file to parallelise the process
#' @param manifest_FILE Path of the manifest.txt file
#' @param test_sample_chunks_DIR Directory to the .txt files with chunks of test sample IDs
#' @param control_sample_chunks_DIR Directory to the .txt files with chunks of control sample IDs
#' @param source_of_disease_haplotypes Are the disease haplotypes sourced from a "database" or from a "directory"? If from a directory, all the database-related parameters must be set to "invalid". db_port="invalid",db_host="invalid",db_password="invalid",db_name="invalid",db_unix_socket="invalid"
#' @param db_port Network port of the FoundHaplo database, "invalid" if disease haplotypes are sourced from a directory  
#' @param db_host Server to the running FoundHaplo database instance, "invalid" if disease haplotypes are sourced from a directory 
#' @param db_password Password of the remote user, "invalid" if disease haplotypes are sourced from a directory 
#' @param db_name Name of the FoundHaplo database, default is FoundHaploDB, "invalid" if disease haplotypes are sourced from a directory 
#' @param db_unix_socket Path to the unix socket file, default is $FoundHaplo_database_DIR/mysql/run/mysqld/mysqld.sock, "invalid" if disease haplotypes are sourced from a directory 
#' @param DCV Name of the disease-causing variant of interest, i.e. FAME1.chr8.119379052. Use the OMIM abbreviation for the disease. (type \code{"character"})
#' @param minor_allele_cutoff The minimum minor allele frequency of SNPs allowed, we recommend this to be 0 (type \code{"numeric"})
#' @param gen_allele_mismatch_rate Genotype and imputation error rate allowed, default is 0.1
#' @param MA_cutoff Moving average threshold for allowing genotype and imputation errors (derived based on simulation studies), default is -0.4
#' @param meiosis Estimated number of meiosis between disease-test pair, default is 1
#' @param imputation_quality_score_cutoff_test Minimum allowed imputation quality which is R-squared. Recommend to use 0.3 if the cohort has >100 samples ; 0 otherwise (type \code{"numeric"})
#' @param frequency_type Population of the test cohort i.e. one of EUR, AMR, SAS, EAS, AFR etc (type \code{"character"})
#' @param geneticMap_DIR Directory to genetic_map_HapMapII_GRCh37 location (type \code{"character"})
#' @param disease_files_DIR Directory of the disease haplotype VCFs for a single disease variant(type \code{"character"})
#' @param test_file File Path to the test cohort VCF (type \code{"character"})
#' @param test_name A meaningful name for the test cohort  (type \code{"character"})
#' @param controls_file_DIR Directory where the 1000 Genomes control files are stored  (type \code{"character"})
#' @param save_report_DIR Directory to save the required details of the IBD sharing to analyze later  (type \code{"character"})
#' @param TEMP_DIR Directory to save the temporary files  (type \code{"character"})
#' @return Write a manifest.txt file that includes all the parameters required to run FoundHaplo. Each line of the manifest.txt file can be submitted as a seperate job parallely.
#'
#' @examples
#' orig_dir <- getwd()
#' setwd(tempdir())
#' file.remove(list.files())
#' unlink(paste0(temp_dir,"/1/"), recursive = T)
#' unlink(paste0(temp_dir,"/2/"), recursive = T)
#' unlink(paste0(temp_dir,"/3/"), recursive = T)
#' unlink(paste0(temp_dir,"/4/"), recursive = T)
#' if(!dir.exists(paste0(tempdir(), "/1"))){dir.create(paste0(tempdir(), "/1"))} #to save disease haplotypes
#' if(!dir.exists(paste0(tempdir(), "/2"))){dir.create(paste0(tempdir(), "/2.test"))} #to save test_list of sample names in .txt files
#' if(!dir.exists(paste0(tempdir(), "/2"))){dir.create(paste0(tempdir(), "/2.controls"))} #to save test_list of sample names in .txt files
#' if(!dir.exists(paste0(tempdir(), "/3"))){dir.create(paste0(tempdir(), "/3"))}  # controls_file_DIR
#' if(!dir.exists(paste0(tempdir(), "/4"))){dir.create(paste0(tempdir(), "/4"))} # to save manifest.txt
#' temp_dir <- tempdir() # To carryout the main workload
#' library(vcfR)
#' write.vcf(FAME1_disease_cohort,paste0(temp_dir,"/","FAME1_disease_cohort.vcf.gz"))
#' sample_info=data.frame(rbind(c("HG00362_1,HG00362_2","duo"),c("NA11920,Affected_parent_NA11920,Unaffected_parent_NA11920","trio"),c("HG00313_1,HG00313_2","duo")))
#' write.table(sample_info,paste0(temp_dir,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' Phasing_by_pedigree(input_vcf = paste0(temp_dir,"/FAME1_disease_cohort",".vcf.gz"),
#'                   dir_output = paste0(temp_dir,"/1"),
#'                   sample_info_file = paste0(temp_dir,"/","sample_info.txt"))
#' write.vcf(FAME1_test_cohort,paste0(temp_dir,"/","FAME1_test_cohort.vcf.gz"))
#' write.table(file00,paste0(temp_dir,"/2.test/","file00",".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' write.vcf(FAME1_control_cohort,paste0(temp_dir,"/3/","FAME1.chr8.vcf.gz"))
#' control_sample_names=colnames(FAME1_control_cohort@gt)[2:ncol(FAME1_control_cohort@gt)]
#' control_sample_names_chunk=split(control_sample_names, ceiling(seq_along(control_sample_names)/10))
#' for(chunk in 1:length(lengths(control_sample_names_chunk)))
#' {
#' write.table(control_sample_names_chunk[[chunk]],paste0(temp_dir,"/2.controls/","file",chunk,".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' }
#' Create_jobs_DB(manifest_FILE=paste0(temp_dir,"/4/manifest.txt"),test_sample_chunks_DIR=paste0(temp_dir,"/2.test"),control_sample_chunks_DIR=paste0(temp_dir,"/2.controls"),source_of_disease_haplotypes="directory",db_port="invalid",db_host="invalid",db_password="invalid",db_name="invalid",db_unix_socket="invalid",minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR=temp_dir,disease_files_DIR=paste0(temp_dir,"/1"),test_file=paste0(temp_dir,"/","FAME1_test_cohort.vcf.gz"),test_name="FAME1_example_test_cohort",controls_file_DIR=paste0(temp_dir,"/3"),save_report_DIR=paste0(temp_dir,"/4"),TEMP_DIR=temp_dir)
#' print("Example content of a manifest.txt file is below")
#' read.delim(paste0(temp_dir,"/4/manifest.txt"),header=FALSE)
#' 
#' 
#' 
Create_jobs=function(manifest_FILE,test_sample_chunks_DIR,control_sample_chunks_DIR,source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff,gen_allele_mismatch_rate,MA_cutoff,meiosis,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,controls_file_DIR,save_report_DIR,TEMP_DIR)
{ 
  #test
  test.txt=list.files(test_sample_chunks_DIR,full.names = TRUE)
  #specify parameters that should be given to the function Generate_FH_score
  test_cohort_entires=expand.grid(source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff,gen_allele_mismatch_rate,MA_cutoff,meiosis,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,test.txt,"test",controls_file_DIR,save_report_DIR,TEMP_DIR)
  
  #controls
  controls.txt=list.files(control_sample_chunks_DIR,full.names = TRUE)
  #specify parameters that should be given to the function Generate_FH_score
  control_cohort_entires=expand.grid(source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff,gen_allele_mismatch_rate,MA_cutoff,meiosis,imputation_quality_score_cutoff_test,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name,controls.txt,"controls",controls_file_DIR,save_report_DIR,TEMP_DIR)
  
  
  manifest.txt_entries=rbind(test_cohort_entires,control_cohort_entires)
  
  #specify path to save the manifest.txt file
  write.table(manifest.txt_entries,manifest_FILE,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)
  
  print(paste0("writing manifest.txt file to ",manifest_FILE))
  
  
  
}
