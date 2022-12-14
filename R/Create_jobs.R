#' Create manifest.txt to run batch jobs
#'
#' @description
#' Takes all the input parameters required to run FoundHaplo and create manifest.txt file to parallelise the process
#' @param path_manifest Path of the manifest.txt file
#' @param path_test_sample_chunks Path to the .txt files with chunks of test sample IDs
#' @param path_control_sample_chunks Path to the .txt files with chunks of control sample IDs
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
#' @param dir_TEMP Directory path to save the temporary files  (type \code{"character"})
#' @return Write a manifest.txt file that includes all the parameters required to run FoundHaplo. Each line of the manifest.txt file can be submitted as a seperate job parallely.
#'
#' @export
#' @examples
#' add example

Create_jobs=function(path_manifest,path_test_sample_chunks,path_control_sample_chunks,DCV,minor_allele_cutoff,imputation_quality_score_cutoff_test,frequency_type,dir_geneticMap,dir_disease_files,test_file,test_name,dir_controls_file,dir_to_save_report,dir_TEMP)
{

    #test
    test_list=list.files(path_test_sample_chunks,full.names = TRUE)
    #specify parameters that should be given to the function Generate_FH_score
    test_cohort_entires=expand.grid(DCV,minor_allele_cutoff,imputation_quality_score_cutoff_test,frequency_type,dir_geneticMap,dir_disease_files,test_file,test_name,test_list,"test",dir_controls_file,dir_to_save_report,dir_TEMP)


    #controls
    controls_list=list.files(path_control_sample_chunks,full.names = TRUE)
    #specify parameters that should be given to the function Generate_FH_score
    control_cohort_entires=expand.grid(DCV,minor_allele_cutoff,imputation_quality_score_cutoff_test,frequency_type,dir_geneticMap,dir_disease_files,test_file,test_name,controls_list,"controls",dir_controls_file,dir_to_save_report,dir_TEMP)

    manifest.txt_entries=rbind(test_cohort_entires,control_cohort_entires)

    #specify path to save the manifest.txt file
    write.table(manifest.txt_entries,path_manifest,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)



}



