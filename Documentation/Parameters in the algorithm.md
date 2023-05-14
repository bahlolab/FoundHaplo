Generate_FH_score is the main wrapper function that generates FH score values for each test/control - disease pair by taking user input shown below.

* Note: Generate_FH_score function easily works on the command line as it uses system() function to query from VCFtools and BCFtools. Hence recomended to use a nextflow pipeline as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing.md).

```R
Generate_FH_score(DCV="FAME1.chr8.119379052",minor_allele_cutoff=0,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR="FoundHaplo_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37",disease_files_DIR="FoundHaplo_PATH/input_files/input_vcf_data/disease_haplotypes",test_file="FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz.imputed.trimmed.vcf.gz",test_name="example_test",test_list="FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/samples/samples.txt",data_type="test",controls_file_DIR="FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR
",save_report_DIR="FoundHaplo_PATH/results/output",TEMP_DIR="FoundHaplo_PATH/temp")
```

# All the parameters that user has to specify are described below

1. **DCV** Name of the disease causing variant of interest i.e FAME1.chr8.119379052 
2. **minor_allele_cutoff** The minimum minor allele frequncy of SNPs allowed, we recommend this to be 0 
3. **imputation_quality_score_cutoff_test** Minimum allowed imputation quality cut off, which is R-squared for the test cohort. Recommend to use 0.3 if the cohort has >100 samples ; 0 otherwise 
4. **frequency_type** population of the test cohort i.e one of EUR,AMR,SAS,EAS,AFR etc 
5. **geneticMap_DIR** directory path to genetic_map_HapMapII_GRCh37 files, which are in FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37/
6. **disease_files_DIR** directory path of the disease haplotype VCFs gzipped, which is FoundHaplo/input_files/input_vcf_data/disease_haplotypes.
7. **test_file** path of the test cohort file gzipped, which is in FoundHaplo/input_files/input_vcf_data/test_cohort.
8. **test_name** meaningful name for the test cohort 
9. **test_list** path to a .txt file that includes sample names from the test/control cohort. If running in parallel, We recommend a list of 1000 sample names if the test cohort has more than 50,000 samples and a list of 100 sample names otherwise. List of sample names in a test/control cohort can be split into chunks of samples as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing.md). Sample names for the test cohort are in /wehisan/bioinf/lab_bahlo/users/robertson.e/FH_test_1/FoundHaplo/input_files/input_vcf_data/test_cohort and sample names for the control cohort are in FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population.  
10. **data_type** "test" or "controls"
11. **controls_file_DIR** directory where the 1000 Genomes gzipped VCF control files are stored. Select the control population cohort the same as the test cohort i.e FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR.
12. **save_report_DIR** directory path to save the output of FoundHaplo IBD sharing for further analysis
13. **TEMP_DIR** directory path to save the temporary files

The function returns all the details of IBD sharing for each test/control sample and will be saved in a seperate tab delimitted text file in dir_to_save_report location, with below columns:

1. **data_type** "test" or "control" 
2. **test_name** meaningful name for the test cohort 
3. **frequency_type** population of the test cohort i.e one of EUR,AMR,SAS,EAS,AFR etc 
4. **minor_allele_cutoff**
5. **imputation_quality_score_cutoff_test**
6. **DCV** 
7. **disease_haplotype**
8. **test_control_individual**
9. **FH_score** 
10. **left_LLR** log likelihood ratio to the left of the DCV in the markov chain
11. **right_LLR** log likelihood ratio to the right of the DCV in the markov chain
12. **total_cM_sharing**
13. **total_left_cM_sharing**
14. **total_right_cM_sharing**
15. **number_of_allele_mismatches_in_the_markov_chain** 
16. **number_of_markers_in_the_markov_chain** 
17. **numer_of_haplotype_switches_in_the_markov_chain** 
18. **snp_density_in_data_file** i.e number of SNPs in 1cM in the inut files
19. **total_number_of_markers_in_data_file**
20. **total_cM_span_of_data_file**

Name of each text file will correspond to a single job sumbitted by the pipeline explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parallel%20processing.md) i.e. data_type.test_name.DCV.disease_individual.test_individual.frequency_type.imputation_quality_score_cutoff_test.txt.

Concatenate all the .txt files in dir_to_save_report and generate a single text file for further analysis as below.

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
cat "$FoundHaplo_DIR/results/output/"*.txt > $FoundHaplo_DIR/results/FH_IBD_scores/results.txt 
```

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).


