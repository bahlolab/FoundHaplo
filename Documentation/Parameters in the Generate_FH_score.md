[Generate_FH_score](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) is the main wrapper function that generates FH score values for each test/control - disease pair using known disease haplotypes (either sourced from a database or a directory) by taking user input shown below.

* Note: Generate_FH_score function easily works on the command line as it uses the system() function to query from VCFtools and BCFtools. Hence recommended to use a nextflow pipeline as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing%20with%20Nextflow.md).

```R
Generate_FH_score(source_of_disease_haplotypes="directory",db_port="invalid",db_host="invalid",db_password="invalid",db_name="invalid",db_unix_socket="invalid",DCV="FAME1.chr8.119379052",minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR="/mypath/FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37",disease_files_DIR="/mypath/FoundHaplo/input_files/input_vcf_data/disease_haplotypes",test_file="/mypath/FoundHaplo/input_files/input_vcf_data/test_cohort/imputed_phased_FAME1_test_cohort.snp.0.98.sample.0.98.chr8.imputed.trimmed.vcf.gz",test_name="example_test",test_list="/mypath/FoundHaplo/input_files/input_vcf_data/test_cohort/samples.txt",data_type="test",controls_file_DIR="/mypath/FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR",save_report_DIR="/mypath/FoundHaplo/results/output",TEMP_DIR="/mypath/FoundHaplo/temp")
```

# All the parameters that the user has to specify are described below

1. **source_of_disease_haplotypes** Are the disease haplotypes sourced from a "database" or from a "directory"? If from a directory, all the database-related parameters must be set to "invalid". db_port="invalid",db_host="invalid",db_password="invalid",db_name="invalid",db_unix_socket="invalid". "database" works ONLY if you have installed the FoundHaplo database using instructions [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20database%20with%20known%20disease%20haplotypes.md)
2. **db_port** Network port of the FoundHaplo database, "invalid" if disease haplotypes are sourced from a directory 
3. **db_host** Server to the running FoundHaplo database instance, "invalid" if disease haplotypes are sourced from a directory
4. **db_password** Password of the remote user, "invalid" if disease haplotypes are sourced from a directory
5. **db_name Name** of the FoundHaplo database, default is FoundHaploDB, "invalid" if disease haplotypes are sourced from a directory
6. **db_unix_socket** Path to the unix socket file, default is $FoundHaplo_database_DIR/mysql/run/mysqld/mysqld.sock, "invalid" if disease haplotypes are sourced from a directory
7. **DCV** Name of the disease-causing variant of interest, i.e. FAME1.chr8.119379052. Use the OMIM abbreviation for the disease.
8. **minor_allele_cutoff** The minimum minor allele frequency of SNPs allowed; we recommend this to be 0 
9. **gen_allele_mismatch_rate** Genotype and imputation error rate allowed; default is 0.1
10. **MA_cutoff** Moving average threshold for allowing genotype and imputation errors (derived based on simulation studies), default is -0.4
11. **meiosis** Estimated number of meiosis between disease-test pair; default is 1
12. **imputation_quality_score_cutoff_test** Minimum allowed imputation quality cut-off, which is R-squared for the test cohort. Recommend to use 0.3 if the cohort has >100 samples ; 0 otherwise 
13. **frequency_type** population of the test cohort, i.e. one of EUR, AMR, SAS, EAS, AFR etc 
14. **geneticMap_DIR** directory path to genetic_map_HapMapII_GRCh37 files, which are in FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37/
15. **disease_files_DIR** directory of the disease haplotype VCFs for a single disease variant, which is FoundHaplo/input_files/input_vcf_data/disease_haplotypes by default. "invalid" if disease haplotypes are sourced from a database (type \code{"character"})
16. **test_file** path of the test cohort file gzipped, which is in FoundHaplo/input_files/input_vcf_data/test_cohort.
17. **test_name** meaningful name for the test cohort 
18. **test_list** path to a .txt file that includes sample names from the test/control cohort. If running in parallel, We recommend a list of 1000 sample names. List of sample names in a test/control cohort can be split into chunks of samples as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing%20with%20Nextflow.md). Sample names for the example dataset of the test cohort are in FoundHaplo/input_files/input_vcf_data/test_cohort and sample names for the example dataset of the control cohort are in FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population.  
19. **data_type** "test" or "controls"
20. **controls_file_DIR** directory, where the 1000 Genomes gzipped VCF control files are stored. Select the control population cohort the same as the test cohort i.e FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR.
21. **save_report_DIR** directory path to save the output of FoundHaplo IBD sharing for further analysis
22. **TEMP_DIR** directory path to save the temporary files

The function returns all the details of IBD sharing for each test/control sample and will be saved in a separate tab-delimited text file in the save_report_DIR location, with the below columns:

1. **data_type** "test" or "control" 
2. **test_name** meaningful name for the test cohort 
3. **frequency_type** population of the test cohort i.e. one of EUR, AMR, SAS, EAS, AFR etc 
4. **minor_allele_cutoff**
5. **imputation_quality_score_cutoff_test**
6. **DCV** 
7. **disease_haplotype**
8. **test_control_individual**
9. **FH_score** 
10. **left_LLR** log-likelihood ratio to the left of the DCV in the Markov chain
11. **right_LLR** log-likelihood ratio to the right of the DCV in the Markov chain
12. **total_cM_sharing**
13. **total_left_cM_sharing**
14. **total_right_cM_sharing**
15. **number_of_allele_mismatches_in_the_markov_chain** 
16. **number_of_markers_in_the_markov_chain** 
17. **numer_of_haplotype_switches_in_the_markov_chain** 
18. **snp_density_in_data_file** i.e number of SNPs in 1cM in the input files
19. **total_number_of_markers_in_data_file**
20. **total_cM_span_of_data_file**

The name of each text file corresponds to a single run of FoundHaplo (or a job if the pipeline explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing%20with%20Nextflow.md) is used)
i.e. data_type.test_name.DCV.disease_individual.test_individual.frequency_type.imputation_quality_score_cutoff_test.txt.

Concatenate all the .txt files in save_report_DIR and generate a single text file for further analysis as below.

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
cat "$FoundHaplo_DIR/results/output/"*.txt > $FoundHaplo_DIR/results/FH_IBD_scores/results.txt 
```

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).


