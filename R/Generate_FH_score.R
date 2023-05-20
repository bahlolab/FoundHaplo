#' A wrapper function to generate FH score and IBD details for each test/control - disease pair using known disease haplotypes sourced from a directory or from a database of disease haplotypes
#'
#' @description
#' This function save FH score and IBD details for each test/control - disease pair in a .txt file in "save_report_DIR"
#' Run this function on R command line since system() command used in the function to utilise vcfools/bcftools may not work on R studio.
#' Ethnicity of the test cohort should be decided in advance EUR,AMR,SAS,EAS,AFR etc
#' Make sure that imputation quality score R2 or R-squared is the third field of the INFO column of the test VCF file.
#' Genome built should be GRCh37 and genetic map in "geneticMap_DIR" must have chromosomes named with "chr" prefix, add the prefix if needed
#' @param source_of_disease_haplotypes Are the disease haplotypes are sourced from a "database" or from a "directory"?. If from a directory, all the database related parameters must be set to NA. db_port=NA,db_host=NA,db_password=NA,db_name=NA,db_unix_socket=NA
#' @param db_port Network port of the FoundHaplo database 
#' @param db_host Server to the running FoundHaplo database instance
#' @param db_password Password of the remote user
#' @param db_name Name of the FoundHaplo database, default is FoundHaploDB
#' @param db_unix_socket Path to the unix socket file, default is $FoundHaplo_database_DIR/mysql/run/mysqld/mysqld.sock
#' @param DCV Name of the disease causing variant of interest i.e FAME1.chr8.119379052 (type \code{"character"})
#' @param minor_allele_cutoff The minimum minor allele frequency of SNPs allowed, we recommend this to be 0 (type \code{"numeric"})
#' @param gen_allele_mismatch_rate Genotype and imputation error rate allowed, default is 0.1
#' @param MA_cutoff Moving average threshold for allowing genotype and imputation errors (derived based on simulation studies), default is -0.4
#' @param meiosis Estimated number of meiosis between disease-test pair, default is 1
#' @param imputation_quality_score_cutoff_test Minimum allowed imputation quality which is R-squared. Recommend to use 0.3 if the cohort has >100 samples ; 0 otherwise (type \code{"numeric"})
#' @param frequency_type Population of the test cohort i.e one of EUR,AMR,SAS,EAS,AFR etc (type \code{"character"})
#' @param geneticMap_DIR Directory to genetic_map_HapMapII_GRCh37 location (type \code{"character"})
#' @param disease_files_DIR directory of the disease haplotype VCFs for a single disease variant(type \code{"character"})
#' @param test_file File path to the test cohort VCF (type \code{"character"})
#' @param test_name meaningful name for the test cohort  (type \code{"character"})
#' @param test_list .txt File path to the file with chunk of test/control samples names to be analysed from the test/control cohort  (type \code{"character"})
#' @param data_type "test" or "control (type \code{"character"})
#' @param controls_file_DIR Directory where the 1000 Genomes control files are stored  (type \code{"character"})
#' @param save_report_DIR Directory to save the required details of the IBD sharing to analyze later  (type \code{"character"})
#' @param TEMP_DIR Directory to save the temporary files  (type \code{"character"})
#' @return All details of IBD sharing for each test/control sample will be saved in tab delimited text files in save_report_DIR location, with below columns :
#' name of each text file will be data_type.test_name.DCV.disease_individual.test_individual.frequency_type.imputation_quality_score_cutoff_test.txt
#' Each .txt file in save_report_DIR location will correspond to each disease haplotype tested
#'
#' \enumerate{
#' \item data_type, "test" or "control (type \code{"character"})
#' \item test_name, meaningful name for the test cohort (type \code{"character"})
#' \item frequency_type, population of the test cohort i.e one of EUR,AMR,SAS,EAS,AFR etc (type \code{"character"})
#' \item minor_allele_cutoff (type \code{"numeric"})
#' \item imputation_quality_score_cutoff_test (type \code{"numeric"})
#' \item DCV (type \code{"character"})
#' \item disease_haplotype (type \code{"character"})
#' \item test_individual (type \code{"character"})
#' \item FH_score (type \code{"numeric"})
#' \item left_LLR (type \code{"numeric"})
#' \item right_LLR (type \code{"numeric"})
#' \item total_cM_sharing (type \code{"numeric"})
#' \item total_left_cM_sharing (type \code{"numeric"})
#' \item total_right_cM_sharing (type \code{"numeric"})
#' \item number_of_allele_mismatches_in_the_markov_chain (type \code{"numeric"})
#' \item number_of_markers_in_the_markov_chain (type \code{"numeric"})
#' \item numer_of_haplotype_switches_in_the_markov_chain (type \code{"numeric"})
#' \item snp_density_in_data_file (type \code{"numeric"})
#' \item total_number_of_markers_in_data_file (type \code{"numeric"})
#' \item total_cM_span_of_data_file (type \code{"numeric"})
#' }
#' @import pryr
#' @import stringi
#' @import stringr
#' @import tibble
#' @import dplyr
#' @import data.table
#' @import R.utils
#' @import vcfR
#' @import DBI
#' @import RMariaDB
#' @export
#' @examples
#' orig_DIR <- getwd()
#' setwd(tempdir())
#' file.remove(list.files())
#' if(!dir.exists(paste0(tempdir(), "/1"))){dir.create(paste0(tempdir(), "/1"))} #to save disease haplotypes
#' if(!dir.exists(paste0(tempdir(), "/2"))){dir.create(paste0(tempdir(), "/2"))} #to save test_list of sample names in .txt files
#' if(!dir.exists(paste0(tempdir(), "/3"))){dir.create(paste0(tempdir(), "/3"))}  # controls_file_DIR
#' if(!dir.exists(paste0(tempdir(), "/4"))){dir.create(paste0(tempdir(), "/4"))} # to save IBD results
#' temp_DIR <- tempdir() # To carryout the main workload
#' library(vcfR)
#' write.vcf(FAME1_disease_cohort,paste0(temp_DIR,"/","FAME1_disease_cohort.vcf.gz"))
#' sample_info=data.frame(rbind(c("HG00362_1_HG00362_1,HG00362_2_HG00362_2","duo"),c("NA11920_NA11920,Affected_parent_NA11920_Affected_parent_NA11920,Unaffected_parent_NA11920_Unaffected_parent_NA11920","trio"),c("HG00313_1_HG00313_1,HG00313_2_HG00313_2","duo")))
#' write.table(sample_info,paste0(temp_DIR,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' Phasing_by_pedigree(input_vcf = paste0(temp_DIR,"/FAME1_disease_cohort",".vcf.gz"),
#'                   output_DIR = paste0(temp_DIR,"/1"),
#'                   sample_info_file = paste0(temp_DIR,"/","sample_info.txt"))
#' write.vcf(FAME1_test_cohort,paste0(temp_DIR,"/","FAME1_test_cohort.vcf.gz"))
#' write.table(file00,paste0(temp_DIR,"/2/","file00",".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' write.vcf(FAME1_control_cohort,paste0(temp_DIR,"/3/","FAME1.chr8.vcf.gz"))
#' write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' Generate_FH_score(source_of_disease_haplotypes="directory",db_port=NA,db_host=NA,db_password=NA,db_name=NA,db_unix_socket=NA,DCV="FAME1.chr8.119379052.",minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR=temp_DIR,disease_files_DIR=paste0(temp_DIR,"/1"),test_file=paste0(temp_DIR,"/","FAME1_test_cohort.vcf.gz"),test_name="FAME1_example_test_cohort",test_list=paste0(temp_DIR,"/2/","file00.txt"),data_type="test",controls_file_DIR=paste0(temp_DIR,"/3"),save_report_DIR=paste0(temp_DIR,"/4"),temp_DIR)
#' setwd(paste0(temp_DIR,"/4"))
#' read.delim(list.files(paste0(temp_DIR,"/4"))[1],header=FALSE)
#' setwd(orig_DIR)

Generate_FH_score=function(source_of_disease_haplotypes,db_port,db_host,db_password,db_name,db_unix_socket,DCV,minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type,geneticMap_DIR,disease_files_DIR,test_file,test_name="test",test_list,data_type,controls_file_DIR,save_report_DIR,TEMP_DIR)
{
  ###
  if(!dir.exists(geneticMap_DIR)){stop("geneticMap_DIR does not exist")}
  
  g1=gen_allele_mismatch_rate # replacing long variable name (g1=gen_allele_mismatch_rate)
  r2=imputation_quality_score_cutoff_test # # replacing long variable name (r2=imputation_quality_score_cutoff_test)
  rand_string=system(paste0("echo $RANDOM | md5sum | head -c 32"),intern = TRUE)
  rand_string
  
  if(source_of_disease_haplotypes=="database")
  {
  
    if(class(db_port)!="numeric"){stop("db_port must be a number")}
    if(!file.exists(db_unix_socket)){stop("db_unix_socket does not exist")}
    db = dbConnect(RMariaDB::MariaDB(),bigint = 'integer',port=db_port,host=db_host,user ='remote_usr',password=db_password,dbname=db_name,unix.socket=db_unix_socket)
    
    #select mutation_id
    PathogenicMutations=dbSendQuery(db, paste0("SELECT * FROM PathogenicMutations where disease_id=","\"",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),"\"",";"))
    PathogenicMutations <- dbFetch(PathogenicMutations,)
    mutation_num=PathogenicMutations$mutation_id
    
    # select individuals with the given mutation_id
    list_of_disease_individuals=dbSendQuery(db, paste0("SELECT * FROM IndividualsWithKnownMutations where mutation_id=","\"",mutation_num,"\"",";"))
    list_of_disease_individuals <- dbFetch(list_of_disease_individuals,)
    list_of_disease_individuals=list_of_disease_individuals$individual_id
    
    # select sample ids for respective individuals
    
    query_command=capture.output(cat(paste0("SELECT * FROM Samples where individual_id"," in ","("),paste0(list_of_disease_individuals[1:(length(list_of_disease_individuals)-1)],","),list_of_disease_individuals[length(list_of_disease_individuals)],paste0(")",";")))
    sample_id=dbSendQuery(db,query_command)
    sample_id=dbFetch(sample_id,)
    sample_num=sample_id$sample_id
    loop_over_disease_haplotypes=sample_num
  }
  
  #add disease sample names in order
  if(source_of_disease_haplotypes=="directory")
  {
    if(db_port!="NA"){stop("db_port must be NA when source_of_disease_haplotypes is a directory")}
    if(db_host!="NA"){stop("db_host must be NA when source_of_disease_haplotypes is a directory")}
    if(db_password!="NA"){stop("db_password must be NA when source_of_disease_haplotypes is a directory")}
    if(db_name!="NA")){stop("db_name must be NA when source_of_disease_haplotypes is a directory")}
    if(db_unix_socket!="NA")){stop("db_unix_socket must be NA when source_of_disease_haplotypes is a directory")}
    
    list_of_disease_individuals=list.files(disease_files_DIR,full.names = TRUE,pattern = ".vcf")
    loop_over_disease_haplotypes=1:length(list_of_disease_individuals)
  }
  
  
  ###
  
  if(data_type!="controls"  & !grepl(".gz", test_file, fixed = TRUE)) # if we are testing the test individuals of interest
  {
    
    command=paste0("module load bcftools ; bcftools view --samples-file ", test_list," ", test_file, " -Oz -o ",TEMP_DIR,"/",rand_string,".vcf")
    system(command) # invoke a system command that will create a temporary vcf file only with the selected sample ids.
    
    fix_file_to_test <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = c(1:9)) # load the first 9 fixed columns of the vcf file once
  }
  if(data_type!="controls"  & grepl(".gz", test_file, fixed = TRUE)) # if we are testing the test individuals of interest
  {
    
    command=paste0("module load bcftools ; bcftools view --samples-file ", test_list," ", test_file, " -Oz -o ",TEMP_DIR,"/",rand_string,".vcf")
    system(command) # invoke a system command that will create a temporary vcf file only with the selected sample ids.
    
    fix_file_to_test <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = c(1:9)) # load the first 9 fixed columns of the vcf file once
  }
  
  if(data_type=="controls" & !grepl(".gz", test_file, fixed = TRUE))  # if we are testing the controls and the test file is not gzipped
  {
    
    command=paste0("cut -f1-9 ",test_file, " > " ,TEMP_DIR,"/",rand_string,".vcf") # save the first 9 fixed columns of the vcf file in a temporary file
    system(command)
    
    fix_file_to_test <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = c(1:9)) # load the temporary file
    
  }
  
  if(data_type=="controls" & grepl(".gz", test_file, fixed = TRUE)) # if we are testing the controls and the test file is gzipped
  {
    
    command=paste0("zcat ",test_file," | ","cut -f1-9 ", " > " ,TEMP_DIR,"/",rand_string,".vcf") # save the first 9 fixed columns of the vcf file in a temporary file
    system(command)
    
    fix_file_to_test <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = c(1:9)) # load the temporary file
    
  }
  
  # Filtering for imputation_quality_score_cutoff_test
  if(r2>0)
  {
    R2=strsplit(fix_file_to_test$INFO,";",fixed=TRUE) # save all R^2 values
    
    keep_R2=grep("R2",R2)
    
    fix_file_to_test=fix_file_to_test[keep_R2,]
    
    R2=strsplit(fix_file_to_test$INFO,";",fixed=TRUE) # save all R^2 values
    R2_list=vector()
    for(ii in 1:length(lengths(R2))){R2_list[ii]=R2[[ii]]  %>%  str_subset(pattern = "R2")}
    R2=R2_list
    R2=sapply(strsplit(R2,"=",fixed=TRUE),"[[",2)
    
    fix_file_to_test$R2=R2
    fix_file_to_test$R2=as.numeric(as.character(fix_file_to_test$R2))
    fix_file_to_test=subset(fix_file_to_test,fix_file_to_test$R2>=r2) #remove markers with R^2 < r2
    fix_file_to_test=as.data.frame(fix_file_to_test)
    fix_file_to_test=fix_file_to_test[,!(names(fix_file_to_test) %in% "R2")] # remove temporary column R2
  }
  # Filtering for imputation_quality_score_cutoff_test
  
  
  for(j in loop_over_disease_haplotypes) # loop over all disease haplotypes in disease_files_DIR
  {
    
    if(source_of_disease_haplotypes=="database")
    {
      Genotypes=dbSendQuery(db, paste0("SELECT * FROM Genotypes where sample_id=","\"",j,"\"",";"))
      Genotypes <- dbFetch(Genotypes,)
      
      
      GeneticMarkers=dbSendQuery(db,"SELECT * FROM GeneticMarkers;") # will have to change this as the database grow, cant load the entire table to R
      GeneticMarkers <- dbFetch(GeneticMarkers,)
      
      Genotypes_markers=merge(GeneticMarkers, Genotypes, by=c("marker_id"))
      
      # Create MAF of the analysis corresponding to the relevant population
      if(frequency_type=="ALL")
      {
        MAF=Genotypes_markers[,"maf_gnomad_ALL"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      if(frequency_type=="AFR")
      {
        MAF=Genotypes_markers[,"maf_gnomad_AFR"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      
      if(frequency_type=="EUR")
      {
        
        MAF=Genotypes_markers[,"maf_gnomad_NFE"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      
      if(frequency_type=="AMR")
      {
        MAF=Genotypes_markers[,"maf_gnomad_AMR"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      
      if(frequency_type=="EAS")
      {
        
        MAF=Genotypes_markers[,"maf_gnomad_EAS"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      
      if(frequency_type=="SAS")
      {
        MAF=Genotypes_markers[,"maf_gnomad_SAS"]
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      
      database_file=as.data.frame(cbind(Genotypes_markers[,c("chr","position_hg19","ref","alt")],MAF,Genotypes_markers[,"genotype"])) # create a data frame with relevant columns for one disease haplotype
      
      disease_individual=dbSendQuery(db, paste0("SELECT * FROM Samples where sample_id=","\"",j,"\"",";"))
      disease_individual <- dbFetch(disease_individual,)
      disease_individual=disease_individual$external_lab_id
      
      database_file=database_file[,c("#CHROM","POS","REF","ALT","MAF","disease_individual")]
    }
    
    if(source_of_disease_haplotypes=="directory")
    {
      if(frequency_type=="ALL")
      {
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM") # database_file : data frame with information regarding a disease haplotype is referred to as a database_file
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_raw"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_raw")}
        }       # if annotated from 100G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "ALL"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "ALL")}
        }
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      if(frequency_type=="AFR")
      {
        
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM")
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_afr"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_afr")}
        }     # if annotated from 1000G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AFR"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AFR")}
        }
        
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      }
      
      if(frequency_type=="EUR")
      {
        
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM")
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_nfe"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_nfe")}
        }     # if annotated from 1000G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "EUR"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "EUR")}
        }
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      
      if(frequency_type=="AMR")
      {
        
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM")
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_amr"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_amr")}
        }     # if annotated from 1000G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AMR"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AMR")}
        }
        
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      
      if(frequency_type=="EAS")
      {
        
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM")
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_eas"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_eas")}
        }     # if annotated from 1000G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "EAS"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "EAS")}
        }
        
        
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      
      if(frequency_type=="SAS")
      {
        MAF_list=vector()
        database_file <-fread(list_of_disease_individuals[j], skip = "#CHROM")
        MAF=strsplit(database_file$INFO,";",fixed=TRUE)
        
        # if annotated from ANNOVAR
        if(length(MAF[[1]]  %>%  str_subset(pattern = "AF_sas"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "AF_sas")}
        }  # if annotated from 1000G
        if(length(MAF[[1]]  %>%  str_subset(pattern = "SAS"))>0){
          for(ii in 1:length(lengths(MAF))){MAF_list[ii]=MAF[[ii]]  %>%  str_subset(pattern = "SAS")}
        }
        
        
        MAF=MAF_list
        MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[",2)
        
        controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
        
      }
      
      
      database_file=as.data.frame(cbind(database_file[,c("#CHROM","POS","REF","ALT")],MAF,database_file[,"h1"])) # create a data frame with relevant columns for one disease haplotype
      database_file=database_file[,c("#CHROM","POS","REF","ALT","MAF","h1")]
      h1=substr(database_file[,"h1"],0,1) # h1 is the disease haplotype. get the first allele of every marker in the VCF file. both alleles are the same in database_file
      
      database_file[,"h1"]=h1
      
      temp=str_count(list_of_disease_individuals[j], "/")
      
      disease_individual = sapply(strsplit(list_of_disease_individuals[j],"/",fixed=TRUE),"[[", temp+1)
      
      colnames(database_file)[ncol(database_file)]=disease_individual # add the name of the disease individual into the database_file
      
    }
    
    database_file$MAF=as.numeric(database_file$MAF)
    database_file=subset(database_file,database_file$MAF>0)
    
    database_file[,"#CHROM"]=fix_file_to_test[,"#CHROM"][1]
    
    # extract positions common to database and test cohort
    common_markers=merge(database_file, fix_file_to_test, by=c("#CHROM","POS", "REF","ALT"))
    fix=common_markers[,c("#CHROM","POS","REF","ALT")]
    
    # extract positions common to database , test cohort and control cohort
    
    controls_file[,"#CHROM"]=fix_file_to_test[,"#CHROM"][1]
    
    common_markers=merge(fix, controls_file, by=c("#CHROM","POS", "REF","ALT"))
    
    fix=common_markers[,c("#CHROM","POS","REF","ALT")] # positions common to database , test cohort and control cohort
    
    
    if(data_type!="controls")
    {
      
      test_samples=read.delim(test_list,header=FALSE)
      test_samples=as.character(test_samples$V1)
      
      test_fix <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = c(1:9))
      test_main <-fread(paste0(TEMP_DIR,"/",rand_string,".vcf"),skip = "#CHROM",select = test_samples)
      
      test_total=as.data.frame(cbind(test_fix,test_main))
      
      test_total[,"#CHROM"]=fix[,"#CHROM"][1]
      
      test_total_1=merge(test_total, fix, by=c("#CHROM","POS", "REF","ALT")) # Extract positions common to database , test cohort and control cohort
      
      test_total_1=Create_hap_VCF(test_total_1) # Separate two haplotypes
      
      test_total_1=Convert_bp_cM(test_total_1,DCV,geneticMap_DIR) # Add cM column
      
      database_file=merge(database_file, fix, by=c("#CHROM","POS", "REF","ALT")) # Extract positions common to database , test cohort and control cohort
      
      database_file=Convert_bp_cM(database_file,DCV,geneticMap_DIR) # Add cM column
      
      database_file=database_file[,which(!colnames(database_file) %in% c("#CHROM","REF","ALT"))]
      
      save_IBD_report_test=function(y) # function to calculate IBD values for all test samples
      {
        
        final_file=as.data.frame(cbind(database_file,test_total_1[,c(y,y+1)]))
        
        final_file=subset(final_file,final_file$MAF>minor_allele_cutoff)
        is.na(final_file) <- final_file=="."
        
        final_file=na.omit(final_file)
        final_file=distinct(final_file,position_cM,.keep_all= TRUE)
        attach(final_file)
        final_file <- final_file[order(POS),]
        detach(final_file)
        
        ####for X chromosome
        if(sum(final_file[1:nrow(final_file),6]=="")==nrow(final_file))
        {
          final_file[,6]=final_file[,5]
        }
        ####for X chromosome
        
        Final_IBD_score <-Calculate_IBD(final_file,DCV,geneticMap_DIR) # Calculate IBD values for a disease-test pair
        Final_IBD_score=paste(Final_IBD_score,collapse="\t")
        ##
        if(grepl("vcf.gz",colnames(final_file)[4], fixed = TRUE))
        {
          disease_individual=str_remove_all(colnames(final_file)[4], ".vcf.gz")
        }else{
          disease_individual=str_remove_all(colnames(final_file)[4], ".vcf")
        }
        ###
        
        test_individual=unlist(strsplit(colnames(final_file)[5], ":",fixed=TRUE))
        test_individual=test_individual[1]
        
        Final_IBD_score=paste0("test","\t",test_name,"\t",frequency_type,"\t",minor_allele_cutoff,"\t",r2,"\t",DCV,"\t",disease_individual,"\t",test_individual,"\t",Final_IBD_score)
        
        return(Final_IBD_score)
      }
    }
    
    
    if(data_type=="controls") # code is similar to if(data_type!="controls")
    {
      
      controls_samples=read.delim(test_list,header=FALSE)
      controls_samples=as.character(controls_samples$V1)
      
      controls_file_fix <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
      
      controls_file <-fread(paste0(controls_file_DIR,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = controls_samples)
      
      controls_file_total=as.data.frame(cbind(controls_file_fix,controls_file))
      
      controls_file_total[,"#CHROM"]=fix[,"#CHROM"][1]
      
      controls_file_total_1=merge(controls_file_total, fix, by=c("#CHROM","POS", "REF","ALT")) # Extract positions common to database , test cohort and control cohort
      
      controls_file_total_1=Create_hap_VCF(controls_file_total_1) # Separate two haplotypes
      
      controls_file_total_1=Convert_bp_cM(controls_file_total_1,DCV,geneticMap_DIR) # Add cM column
      
      
      database_file=merge(database_file, fix, by=c("#CHROM","POS", "REF","ALT"))  # Extract positions common to database , test cohort and control cohort
      
      
      database_file=Convert_bp_cM(database_file,DCV,geneticMap_DIR) # Add cM column
      
      database_file=database_file[,which(!colnames(database_file) %in% c("#CHROM","REF","ALT"))]
      
      save_IBD_report_controls=function(y) #function to calculate IBD values for all control samples
      {
        
        final_file=as.data.frame(cbind(database_file,controls_file_total_1[,c(y,y+1)]))
        
        final_file=subset(final_file,final_file$MAF>minor_allele_cutoff)
        is.na(final_file) <- final_file=="."
        
        final_file=na.omit(final_file)
        final_file=distinct(final_file,position_cM,.keep_all= TRUE)
        attach(final_file)
        final_file <- final_file[order(POS),]
        detach(final_file)
        
        ####for X chromosome
        if(sum(final_file[1:nrow(final_file),6]=="")==nrow(final_file))
        {
          final_file[,6]=final_file[,5]
        }
        ####for X chromosome
        
        Final_IBD_score <-Calculate_IBD(final_file,DCV,geneticMap_DIR) # Calculate IBD values for a disease-control pair
        Final_IBD_score=paste(Final_IBD_score,collapse="\t")
        ##
        if(grepl("vcf.gz",colnames(final_file)[4], fixed = TRUE))
        {
          disease_individual=str_remove_all(colnames(final_file)[4], ".vcf.gz")
        }else{
          disease_individual=str_remove_all(colnames(final_file)[4], ".vcf")
        }
        ###
        
        test_individual=unlist(strsplit(colnames(final_file)[5], ":",fixed=TRUE))
        test_individual=test_individual[1]
        
        Final_IBD_score=paste0("controls","\t",test_name,"\t",frequency_type,"\t",minor_allele_cutoff,"\t",r2,"\t",DCV,"\t",disease_individual,"\t",test_individual,"\t",Final_IBD_score)
        return(Final_IBD_score)
      }
      
      
    }
    
    
    if(data_type=="test") # loop over save_IBD_report_test for all the test samples
    {
      run_test=seq(11,ncol(test_total_1),2)
      
      results=lapply(run_test,save_IBD_report_test)
      
      final_file=data.frame(do.call("rbind",results))
      dim(final_file)
      
      Final_IBD_score_name=paste("test",test_name,frequency_type,minor_allele_cutoff,r2,DCV,disease_individual,sapply(strsplit(test_list,"/"),"[[",str_count(test_list,"/")+1),sep=".",collapse=NULL)
      
      path1<-paste0(save_report_DIR,"/",Final_IBD_score_name)
      
      write.table(final_file,path1,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
      
      
    }
    if(data_type=="controls")  # loop over save_IBD_report_test for all the control samples
    {
      run_controls=seq(11,ncol(controls_file_total_1),2)
      
      results=lapply(run_controls,save_IBD_report_controls)
      
      final_file=data.frame(do.call("rbind",results))
      dim(final_file)
      
      Final_IBD_score_name=paste("controls",test_name,frequency_type,minor_allele_cutoff,r2,DCV,disease_individual,sapply(strsplit(test_list,"/"),"[[",str_count(test_list,"/")+1),sep=".",collapse=NULL)
      
      path1<-paste0(save_report_DIR,"/",Final_IBD_score_name)
      
      write.table(final_file,path1,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
    }
    
    
  }
  
  system(paste0("rm -rf ",TEMP_DIR,"/",rand_string,".vcf")) # delete the temporary vcf file created
 
  if(source_of_disease_haplotypes=="database")
  {
   dbDisconnect(db)
  }
  
}
