## Run FoundHaplo parallelly using Nextflow

1. Make sure you have installed the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```

We recommend using a nextflow pipeline to run the main function [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) in FoundHaplo. Refer https://www.nextflow.io/docs/latest/process.html for more details

2. Use [create_sample_chunks.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/create_sample_chunks.sh) script to split test and control samples into chunks of sample IDs to parallely run chunks of samples in one nextflow job. The script will create a folder named "samples" in the directory where the sample ID files already are. 

```bash
FoundHaplo_PATH=
TEST_SAMPLES_FILE=
CONTROL_SAMPLES_FILE=
CHUNK_SIZE=
$FoundHaplo_PATH/scripts/run_nextflow/create_sample_chunks.sh "$FoundHaplo_PATH" "$FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/samples.txt" "$FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt" "100"
```
Set the variables as below,

* FoundHaplo_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* TEST_SAMPLES_FILE : Path to .txt file with test sample IDs 
* CONTROL_SAMPLES_FILE :  Path to .txt file with control sample IDs 
* CHUNK_SIZE: Number of samples in one chunk. We recommend a CHUNK_SIZE of 1000 sample names if the test cohort has more than 50,000 samples and CHUNK_SIZE of 100 otherwise.

3. Nextflow pipeline requires below scripts and files which are in /FoundHaplo/scripts/run_nextflow/.

* Tab delimitted [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) to read the parameters from.

  manifest.txt file can be easily generated using the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R), which requires all the parameters (except test_list and data_type) in the main R script [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20algorithm.md), and three additional parameters which are,

(i) path_manifest : Path to save the manifest.txt file (/FoundHaplo/scripts/run_nextflow/manifest.txt)

(ii) path_test_sample_chunks : Path of the .txt files with chunks of test sample IDs (FoundHaplo/input_files/input_vcf_data/test_cohort/samples)

(iii) path_control_sample_chunks : Path of the .txt files with chunks of control sample IDs (FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples)

```R
Create_jobs(path_manifest="FoundHaplo_PATH/scripts/run_nextflow/manifest.txt",path_test_sample_chunks="FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/samples",path_control_sample_chunks="FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples",DCV="FAME1.chr8.119379052",minor_allele_cutoff=0,imputation_quality_score_cutoff_test=0,frequency_type="EUR",dir_geneticMap="FoundHaplo_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37",dir_disease_files="FoundHaplo_PATH/input_files/input_vcf_data/disease_haplotypes",test_file="FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz.imputed.trimmed.vcf.gz",test_name="example_test",dir_controls_file="FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR",dir_to_save_report="FoundHaplo_PATH/results",dir_TEMP="FoundHaplo_PATH/temp")
```

* [Args_Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Args_Generate_FH_score.R) to submit parameters as arguments to the main R script, [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R). 

* The nextflow configuration file, [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/nextflow.config).
* The main script to run FoundHaplo jobs in parallel, [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf).

4. Re-write the [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) as necessary.

5. Run all nextflow jobs in the [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) file as below. Each line in the manifest.txt file will be a single independant job that is run parallely.
```bash
module load nextflow
nohup ./run_nextflow.nf
```

Results will be saved in the path given in "dir_to_save_report" in the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R).

6. Concatenate all the .txt files in dir_to_save_report and generate a single text file for further analysis as below.

```bash
FoundHaplo_PATH=
cat "$FoundHaplo_PATH/results/"*.txt > $FoundHaplo_PATH/results/FH_IBD_scores/results.txt 
```

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
