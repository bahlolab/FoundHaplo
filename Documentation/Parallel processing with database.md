## Run FoundHaplo parallelly using Nextflow, when you have FoundHaplo database installed

1. Make sure you have installed the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```
Make sure have installed the FoundHaplo database using instructions [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20database%20with%20known%20disease%20haplotypes.md)

We recommend using a nextflow pipeline to run the main function [Generate_FH_score_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score_DB.R) in FoundHaplo. Refer https://www.nextflow.io/docs/latest/process.html for more details

2. Use [create_sample_chunks.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/create_sample_chunks.sh) script to split test and control samples into chunks of sample IDs to parallely run chunks of samples in one nextflow job. The script will create a folder named "samples" in the directory where the sample ID files already are. 

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
TEST_SAMPLES_FILE=/mypath/FoundHaplo/samples.txt
CONTROL_SAMPLES_FILE=/mypath/FoundHaplo/ethnicity.txt
CHUNK_SIZE=CHUNK_SIZE

$FoundHaplo_DIR/scripts/run_nextflow_with_database/create_sample_chunks.sh $FoundHaplo_DIR $TEST_SAMPLES_FILE $CONTROL_SAMPLES_FILE $CHUNK_SIZE
```
Run example below,

```bash
$FoundHaplo_DIR/scripts/run_nextflow_with_database/create_sample_chunks.sh $FoundHaplo_DIR $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples.txt $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt 100
```

Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e path/FoundHaplo
* TEST_SAMPLES_FILE : Path to .txt file with test sample IDs 
* CONTROL_SAMPLES_FILE :  Path to .txt file with control sample IDs 
* CHUNK_SIZE: Number of samples in one chunk. We recommend a CHUNK_SIZE of 1000 sample names.

3. Nextflow pipeline requires below scripts and files which are in /FoundHaplo/scripts/run_nextflow_with_database/.

* Tab delimitted [manifest_DB.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/manifest_DB.txt) to read the parameters from.

Generate the manifest_DB.txt file using the script [Create_jobs_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/Create_jobs_DB.R), which requires all the parameters in the main R script [Generate_FH_score_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score_DB.R) as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20Generate_FH_score_DB.md), and three additional parameters which are,

(i) manifest_FILE : Path to save the manifest_DB.txt file (/FoundHaplo/scripts/run_nextflow_with_database/manifest_DB.txt)

(ii) test_sample_chunks_DIR : Directory of the .txt files with chunks of test sample IDs (FoundHaplo/input_files/input_vcf_data/test_cohort/samples)

(iii) control_sample_chunks_DIR : Directory of the .txt files with chunks of control sample IDs (FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples)

```R
Create_jobs(manifest_FILE="FoundHaplo_DIR/scripts/run_nextflow_with_database/manifest_DB.txt",test_sample_chunks_DIR="FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples",control_sample_chunks_DIR="FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples",DCV="FAME1.chr8.119379052",minor_allele_cutoff=0,imputation_quality_score_cutoff_test=0,frequency_type="EUR",FoundHaplo_DIR="FoundHaplo_DIR/input_files/public_data/genetic_map_HapMapII_GRCh37",disease_files_DIR="FoundHaplo_DIR/input_files/input_vcf_data/disease_haplotypes",test_file="FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/imputed_phased_FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz.imputed.trimmed.vcf.gz",test_name="example_test",controls_file_DIR="FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR",save_report_DIR="FoundHaplo_DIR/results/output",TEMP_DIR="FoundHaplo_DIR/temp")
```

* [Args_Generate_FH_score_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/Args_Generate_FH_score_DB.R) to submit parameters as arguments to the main R script, [Generate_FH_score_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score_DB.R). 

* The nextflow configuration file, [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/nextflow.config).
* The main script to run FoundHaplo jobs in parallel, [run_nextflow_DB.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/run_nextflow_DB.nf). Re-write the [run_nextflow_DB.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/run_nextflow_DB.nf) as necessary.

4. Run all nextflow jobs in the [manifest_DB.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/manifest_DB.txt) file as below. Each line in the manifest_DB.txt file will be a single independant job that is run parallely.
```bash
module load nextflow
FoundHaplo_DIR=/mypath/FoundHaplo
cd $FoundHaplo_DIR/scripts/run_nextflow_with_database
nohup ./run_nextflow_DB.nf
```

Results will be saved in the path given in "save_report_DIR" (default: /mypath/FoundHaplo_DIR/results/output) in the script [Create_jobs_DB.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow_with_database/Create_jobs_DB.R).

5. Concatenate all the .txt files in save_report_DIR and generate a single text file for further analysis as below.

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
cat $FoundHaplo_DIR/results/output/*.txt > $FoundHaplo_DIR/results/FH_IBD_scores/results.txt 
```

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
