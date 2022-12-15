## Run FoundHaplo parallelly using Nextflow

1. Make sure you have installed the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```

We recommend using a nextflow pipeline to run FoundHaplo. Refer https://www.nextflow.io/docs/latest/process.html for more details

2. Use [create_sample_chunks.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/create_sample_chunks.sh) script to split test and control samples into chunks of sample IDs to parallely run chunks of samples in one nextflow job. The script will create a folder named "samples" in the directory where the sample ID files already are.

```bash
FoundHaplo_PATH=
TEST_SAMPLES_FILE=
CONTROL_SAMPLES_FILE=
CHUNK_SIZE=$4 
$FoundHaplo_PATH/scripts/run_nextflow/create_sample_chunks.sh "$FoundHaplo_PATH" "$FoundHaplo_PATH/input_files/input_vcf_data/test_cohort/samples.txt" "$FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt" "100"
```
Set the variables as below,

* FoundHaplo_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* TEST_SAMPLES_FILE : Path to .txt file with test sample IDs 
* CONTROL_SAMPLES_FILE :  Path to .txt file with control sample IDs 
* CHUNK_SIZE: Number of samples in one chunk, recomended a max of 1000

3. Nextflow pipeline requires two nextflow scripts [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/nextflow.config) and [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) , a tab delimitted [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) to read the parameters from and [Args_Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Args_Generate_FH_score.R) to submit parameters as arguments to the [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R). 

manifest.txt file can be easily generated using the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R), which has all the parameters in [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R) as explained [here](https://github.com/bahlolab/FoundHaplo/edit/main/Documentation/Parameters%20in%20the%20algorithm.md), and three additional parameters which are,

* **path_manifest** Path to save the manifest.txt file
* **path_test_sample_chunks** Path to save the .txt files with chunks of test sample IDs
* **path_control_sample_chunks** Path to save the .txt files with chunks of control sample IDs

4. Re-write the [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) as necessary.

5. Run all nextflow jobs in the [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) file as below. Each line in the manifest.txt file will be a single independant job that is run parallely.
```bash
./run_nextflow.nf
```


Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
