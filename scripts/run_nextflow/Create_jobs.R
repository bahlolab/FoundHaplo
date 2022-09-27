#specify path to all the test cohort .txt files with sample IDs and path to all the control cohort .txt files with sample IDs
sample_files=c(PATH_test_list_test,PATH_test_list_controls)
#specify path to all the test cohort .txt files with sample IDs and path to all the control cohort .txt files with sample IDs


#test
test_list=list.files(sample_files[1],full.names = TRUE)
#specify parameters that should be given to the function Generate_FH_score
test_cohort_entires=expand.grid(DCV,minor_allele_cutoff,imputation_quality_score_cutoff_test,frequency_type,root_geneticMap,path_to_save_report,disease_directory,test_file,test_name,test_list,"test",path_to_controls_file,TEMP_dir)


#controls
test_list=list.files(sample_files[2],full.names = TRUE)
#specify parameters that should be given to the function Generate_FH_score
control_cohort_entires=expand.grid(DCV,minor_allele_cutoff,imputation_quality_score_cutoff_test,frequency_type,root_geneticMap,path_to_save_report,disease_directory,test_file,test_name,test_list,"controls",path_to_controls_file,TEMP_dir)

manifest.txt_entries=rbind(test_cohort_entires,control_cohort_entires)

#specify path to save the manifest.txt file
write.table(file,PATH_manifest,sep="\t",quote=FALSE,col.names = FALSE,row.names = FALSE)





