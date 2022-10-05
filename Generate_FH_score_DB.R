#' A wrapper function to generate FH score and IBD details for each test/control - disease pair
#'
#' @description
#' This function save FH score and IBD details for each test/control - disease pair in a .txt file in "dir_to_save_report"
#' Run this function on R command line since system() command used in the function doesn't work on R studio.
#' Ethnicity of the test cohort should be decided in advance EUR,AMR,SAS,EAS,AFR etc
#' Make sure that imputation quality score R2 or R-squared is the third field of the INFO column of the test VCF file.
#' Genome built should be GRCh37 and genetic map in "dir_geneticMap" must have chromosomes named with "chr" prefix, add the prefix if needed
#' @param DCV Name of the disease causing variant of interest i.e FAME1.chr8.119379052 (type \code{"character"})
#' @param minor_allele_cutoff The minimum minor allele frequncy of SNPs allowed, we recommend this to be 0 (type \code{"numeric"})
#' @param imputation_quality_score_cutoff_test Minimum allowed imputation quality which is R-squared. Recommend to use 0.3 if the cohort has >100 samples ; 0 otherwise (type \code{"numeric"})
#' @param frequency_type Population of the test cohort i.e one of EUR,AMR,SAS,EAS,AFR etc (type \code{"character"})
#' @param dir_geneticMap Directory path to genetic_map_HapMapII_GRCh37 location (type \code{"character"})
#' @param dir_disease_files directory of the disease haplotype VCFs for a single disease variant(type \code{"character"})
#' @param test_file path of the test cohort file (type \code{"character"})
#' @param test_name meaningful name for the test cohort  (type \code{"character"})
#' @param test_list path to a set of 100 test samples from the test cohort  (type \code{"character"})
#' @param data_type "test" or "control (type \code{"character"})
#' @param dir_controls_file Directory where the 1000genome control files are stored  (type \code{"character"})
#' @param dir_to_save_report Directory path to save the required details of the IBD sharing to analyze later  (type \code{"character"})
#' @return All details of IBD sharing for each test/control sample will be saved in a seperate tab delimitted text file in dir_to_save_report location, with below columns :
#' name of each text file will be data_type.test_name.DCV.disease_individual.test_individual.frequency_type.imputation_quality_score_cutoff_test.txt
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
#' @import RMariaDB
#' @import DBI
#' @export
#' @examples


Generate_FH_score_DB=function(db_port,db_host,db_pwd,db_name,unix_socket,DCV,minor_allele_cutoff=0,imputation_quality_score_cutoff_test=0,frequency_type,dir_geneticMap,dir_disease_files,test_file,test_name="test",test_list,data_type,dir_controls_file,dir_to_save_report)
{
  db = dbConnect(RMariaDB::MariaDB(),bigint = 'integer',port=db_port,host=db_host,user ='remote_usr',password=db_pwd,dbname=db_name,unix.socket=unix_socket)

  #select mutation_id
  q1=dbSendQuery(db, paste0("SELECT * FROM PathogenicMutations where disease_id=","\"",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),"\"",";"))
  PathogenicMutations <- dbFetch(q1,)
  mutation_id=PathogenicMutations$mutation_id

  # select individuals with the given mutation_id
  q2=dbSendQuery(db, paste0("SELECT * FROM IndividualsWithKnownMutations where mutation_id=","\"",mutation_id,"\"",";"))
  list_of_disease_individuals <- dbFetch(q2,)
  list_of_disease_individuals=list_of_disease_individuals$individual_id

  # select sample ids for respective individuals

  c1=capture.output(cat(paste0("SELECT * FROM Samples where individual_id"," in ","("),paste0(list_of_disease_individuals[1:(length(list_of_disease_individuals)-1)],","),list_of_disease_individuals[length(list_of_disease_individuals)],paste0(")",";")))
  q3=dbSendQuery(db,c1)
  sample_id=dbFetch(q3,)
  sample_id=sample_id$sample_id

  gen_allele_mismatch_rate = 0.01 # genotype/imputation allele_mismatch rate allowed
  g1=gen_allele_mismatch_rate

  r2=imputation_quality_score_cutoff_test # minimum R^2 value allowed in imputed genotypes

  #add disease sample names in order




  tmpfile <- tempfile('1', tempdir())


  if(data_type!="controls"  & !grepl(".gz", test_file, fixed = TRUE)) # if we are testing the test individuals of interest
  {

    command=paste0("module load bcftools ; bcftools view --samples-file ", test_list," ", test_file, " -Oz -o ",tmpfile,".vcf")
    system(command) # invoke a system command that will create a temporary vcf file only with the selected sample ids.

    fix_file_to_test <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = c(1:9)) # load the first 9 fixed columns of the vcf file once
  }
  if(data_type!="controls"  & grepl(".gz", test_file, fixed = TRUE)) # if we are testing the test individuals of interest
  {

    command=paste0("module load bcftools ; bcftools view --samples-file ", test_list," ", test_file, " -Oz -o ",tmpfile,".vcf")
    system(command) # invoke a system command that will create a temporary vcf file only with the selected sample ids.

    fix_file_to_test <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = c(1:9)) # load the first 9 fixed columns of the vcf file once
  }

  if(data_type=="controls" & !grepl(".gz", test_file, fixed = TRUE)) # if the test file is not gzipped
  {

    command=paste0("cut -f1-9 ",test_file, " > " ,tmpfile,".vcf") # save the first 9 fixed columns of the vcf file in a temporary file
    system(command)

    fix_file_to_test <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = c(1:9)) # load the temporary file

  }

  if(data_type=="controls" & grepl(".gz", test_file, fixed = TRUE)) # if the test file is gzipped
  {

    command=paste0("zcat ",test_file," | ","cut -f1-9 ", " > " ,tmpfile,".vcf") # save the first 9 fixed columns of the vcf file in a temporary file
    system(command)

    fix_file_to_test <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = c(1:9)) # load the temporary file

  }


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
    fix_file_to_test=fix_file_to_test[,!(names(fix_file_to_test) %in% "R2")] # remove temorary column R2



  }



  for(j in sample_id)
  {


    q4=dbSendQuery(db, paste0("SELECT * FROM Genotypes where sample_id=","\"",sample_id,"\"",";"))
    Genotypes <- dbFetch(q4,)


    q5=dbSendQuery(db,"SELECT * FROM GeneticMarkers;") # will have to change this as the database grow, cant load the entire table to R
    GeneticMarkers <- dbFetch(q5,)

    Genotypes_markers=merge(GeneticMarkers, Genotypes, by=c("marker_id"))

#   c1=capture.output(cat(paste0("SELECT * FROM GeneticMarkers where marker_id"," in ","("),paste0(Genotypes$marker_id[1:(length(Genotypes$marker_id)-1)],","),Genotypes$marker_id[length(Genotypes$marker_id)],paste0(")",";")))


    #load the relevant controls file and the relevant population frequency
    if(frequency_type=="ALL")
    {

      MAF=Genotypes_markers[,"maf_gnomad_ALL"]

      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
    }
    if(frequency_type=="AFR")
    {

      MAF=Genotypes_markers[,"maf_gnomad_AFR"]
      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))
    }

    if(frequency_type=="EUR")
    {

      MAF=Genotypes_markers[,"maf_gnomad_NFE"]
      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))

    }

    if(frequency_type=="AMR")
    {

      MAF=Genotypes_markers[,"maf_gnomad_AMR"]
      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))

    }

    if(frequency_type=="EAS")
    {

      MAF=Genotypes_markers[,"maf_gnomad_EAS"]
      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))

    }

    if(frequency_type=="SAS")
    {
      MAF=Genotypes_markers[,"maf_gnomad_SAS"]
      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))

    }

    database_file=as.data.frame(cbind(Genotypes_markers[,c("chr","position_hg19","ref","alt")],MAF,Genotypes_markers[,"genotype"]))

    q6=dbSendQuery(db, paste0("SELECT * FROM Samples where sample_id=","\"",j,"\"",";"))
    disease_individual <- dbFetch(q6,)
    disease_individual=disease_individual$external_lab_id

     colnames(database_file)=c("#CHROM","POS", "REF","ALT","MAF",disease_individual)

    database_file$MAF=as.numeric(database_file$MAF)
    database_file=subset(database_file,database_file$MAF>0)

    database_file[,"#CHROM"]=fix_file_to_test[,"#CHROM"][1]


    # extract positions common to database and test cohort
    common_markers=merge(database_file, fix_file_to_test, by=c("#CHROM","POS", "REF","ALT"))
    fix=common_markers[,c("#CHROM","POS","REF","ALT")]


    # extract positions common to database , test cohort and control cohort
    controls_file[,"#CHROM"]=fix_file_to_test[,"#CHROM"][1]

    common_markers=merge(fix, controls_file, by=c("#CHROM","POS", "REF","ALT"))

    fix=common_markers[,c("#CHROM","POS","REF","ALT")]




    if(data_type!="controls")
    {

      test_samples=read.delim(test_list,header=FALSE)
      test_samples=as.character(test_samples$V1)

      test_fix <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = c(1:9))
      test_main <-fread(paste0(tmpfile,".vcf"),skip = "#CHROM",select = test_samples)

      test_total=as.data.frame(cbind(test_fix,test_main))

      test_total[,"#CHROM"]=fix[,"#CHROM"][1]

      test_total_1=merge(test_total, fix, by=c("#CHROM","POS", "REF","ALT"))

      test_total_1=Create_hap_VCF(test_total_1)

      test_total_1=Convert_bp_cM(test_total_1,DCV,dir_geneticMap)


      database_file=merge(database_file, fix, by=c("#CHROM","POS", "REF","ALT"))


      database_file=Convert_bp_cM(database_file,DCV,dir_geneticMap)

      database_file=database_file[,which(!colnames(database_file) %in% c("#CHROM","REF","ALT"))]



      save_IBD_report_test=function(y) #function to calculate IBD values for all test samples
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


        Final_IBD_score <-Calculate_IBD(final_file,DCV,dir_geneticMap)
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

        Final_IBD_score_name=paste("test",test_name,frequency_type,minor_allele_cutoff,r2,DCV,disease_individual,test_individual,sep=".",collapse=NULL)

        Final_IBD_score=paste0("test","\t",test_name,"\t",frequency_type,"\t",minor_allele_cutoff,"\t",r2,"\t",DCV,"\t",disease_individual,"\t",test_individual,"\t",Final_IBD_score)

        dir_to_save_report=paste0(dir_to_save_report,"/",Final_IBD_score_name,".txt")


        write.table(Final_IBD_score,dir_to_save_report,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
        rm(final_file)


      }
    }


    if(data_type=="controls")
    {

      controls_samples=read.delim(test_list,header=FALSE)
      controls_samples=as.character(controls_samples$V1)


      controls_file_fix <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = c(1:9))

      controls_file <-fread(paste0(dir_controls_file,"/",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 1),".",sapply(strsplit(DCV,".",fixed=TRUE),"[[", 2),".vcf.gz"),skip = "#CHROM", select = controls_samples)


      controls_file_total=as.data.frame(cbind(controls_file_fix,controls_file))

      controls_file_total[,"#CHROM"]=fix[,"#CHROM"][1]

      controls_file_total_1=merge(controls_file_total, fix, by=c("#CHROM","POS", "REF","ALT"))

      controls_file_total_1=Create_hap_VCF(controls_file_total_1)

      controls_file_total_1=Convert_bp_cM(controls_file_total_1,DCV,dir_geneticMap)


      database_file=merge(database_file, fix, by=c("#CHROM","POS", "REF","ALT"))


      database_file=Convert_bp_cM(database_file,DCV,dir_geneticMap)

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


        Final_IBD_score <-Calculate_IBD(final_file,DCV,dir_geneticMap)
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


        Final_IBD_score_name=paste("controls",test_name,frequency_type,minor_allele_cutoff,r2,DCV,disease_individual,test_individual,sep=".",collapse=NULL)

        Final_IBD_score=paste0("controls","\t",test_name,"\t",frequency_type,"\t",minor_allele_cutoff,"\t",r2,"\t",DCV,"\t",disease_individual,"\t",test_individual,"\t",Final_IBD_score)

        dir_to_save_report=paste0(dir_to_save_report,"/",Final_IBD_score_name,".txt")
        write.table(Final_IBD_score,dir_to_save_report,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)

        rm(final_file)

      }




    }


    if(data_type=="test")
    {
      run_test=seq(11,ncol(test_total_1),2)

      sapply(run_test,save_IBD_report_test)
    }
    if(data_type=="controls")
    {
      run_controls=seq(11,ncol(controls_file_total_1),2)

      sapply(run_controls,save_IBD_report_controls)
    }



  }

  system(paste0("rm -rf ",tmpfile,".vcf")) # delete the temporary vcf file created
  dbHasCompleted(q1,q2,q3,q4,q5,q6)
  dbClearResult(q1,q2,q3,q4,q5,q6)
  #close the connection once everythings done
  dbDisconnect(db)


}












