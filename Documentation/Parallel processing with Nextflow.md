## Run FoundHaplo parallelly using a Nextflow pipeline

1. Make sure you have installed the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```

We recommend using a nextflow pipeline to run the main function [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) in FoundHaplo. Refer https://www.nextflow.io/docs/latest/process.html for more details

2. Use [create_sample_chunks.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/create_sample_chunks.sh) script to split test and control samples into chunks of sample IDs to parallely run chunks of samples in one nextflow job. The script will create a folder named "samples" in the directory where the sample ID files already are. 

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
TEST_SAMPLES_FILE=/mypath/FoundHaplo/samples.txt
CONTROL_SAMPLES_FILE=/mypath/FoundHaplo/ethnicity.txt
CHUNK_SIZE=CHUNK_SIZE

$FoundHaplo_DIR/scripts/run_nextflow/create_sample_chunks.sh $FoundHaplo_DIR $TEST_SAMPLES_FILE $CONTROL_SAMPLES_FILE $CHUNK_SIZE
```
Run example below,

```bash
$FoundHaplo_DIR/scripts/run_nextflow/create_sample_chunks.sh $FoundHaplo_DIR $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples.txt $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt 100
```

Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e path/FoundHaplo
* TEST_SAMPLES_FILE : Path to .txt file with test sample IDs, Saved in the same directory as the test file in samples.txt
* CONTROL_SAMPLES_FILE :  Path to .txt file with control sample IDs, Saved in $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population
* CHUNK_SIZE: Number of samples of the test cohort that should be included in one chunk. We recommend a CHUNK_SIZE of 1000 sample names if the test cohort is large (>=10,000) or CHUNK_SIZE of 100 otherwise. Control cohort will be in chunks of 100 samples by default.

3. Nextflow pipeline requires below scripts and files which are in /FoundHaplo/scripts/run_nextflow/.

* Tab delimitted [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) to read the parameters from.

Generate the manifest.txt file in Foundhaplo/scripts/run_nextflow/ using the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R) (not included in the FoundHaplo R package), which requires all the parameters (except for data_type and test_list) in the main R script [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) as explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20Generate_FH_score.md), and three additional parameters which are,

(i) manifest_FILE : Path to save the manifest.txt file (FoundHaplo/scripts/run_nextflow/manifest.txt)

(ii) test_sample_chunks_DIR : Directory of the .txt files with chunks of test sample IDs (FoundHaplo/input_files/input_vcf_data/test_cohort/samples)

(iii) control_sample_chunks_DIR : Directory of the .txt files with chunks of control sample IDs (FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples)

```R
Create_jobs(manifest_FILE="/mypath/FoundHaplo/scripts/run_nextflow/manifest.txt",test_sample_chunks_DIR="/mypath/FoundHaplo/input_files/input_vcf_data/test_cohort/samples",control_sample_chunks_DIR="/mypath/FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/samples",source_of_disease_haplotypes="directory",db_port="invalid",db_host="invalid",db_password="invalid",db_name="invalid",db_unix_socket="invalid",DCV="FAME1.chr8.119379052",minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR="/mypath/FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37",disease_files_DIR="/mypath/FoundHaplo/input_files/input_vcf_data/disease_haplotypes",test_file="/mypath/FoundHaplo/input_files/input_vcf_data/test_cohort/imputed_phased_FAME1_test_cohort.snp.0.98.sample.0.98.chr8.imputed.trimmed.vcf.gz",test_name="example_test",controls_file_DIR="/mypath/FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR",save_report_DIR="/mypath/FoundHaplo/results/output",TEMP_DIR="/mypath/FoundHaplo/temp")
```

* [Args_Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Args_Generate_FH_score.R) to submit parameters as arguments to the main R script, [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R). 

* The nextflow configuration file, [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/nextflow.config).
* The main script to run FoundHaplo jobs in parallel, [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf). Re-write the [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) in Foundhaplo/scripts/run_nextflow/ as necessary.

4. Run all nextflow jobs in the [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) file as below. Each line in the manifest.txt file will be a single independant job that is run parallely.
```bash
module load nextflow
FoundHaplo_DIR=/mypath/FoundHaplo
cd $FoundHaplo_DIR/scripts/run_nextflow
nohup ./run_nextflow.nf
```

Results will be saved in the path given in "save_report_DIR" (default: /mypath/FoundHaplo_DIR/results/output) in the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R).

5. Concatenate all the .txt files in save_report_DIR and generate a single text file for further analysis as below.

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
cat $FoundHaplo_DIR/results/output/*.txt > $FoundHaplo_DIR/results/FH_IBD_scores/results.txt 
```

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
