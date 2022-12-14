## Run FoundHaplo parallelly using Nextflow

1. Install the R package using devtools packages

```R
devtools::install_github("bahlolab/FoundHaplo")
```

We recommend using a nextflow pipeline to run FoundHaplo. Refer https://www.nextflow.io/docs/latest/process.html for more details

2. Use [create_sample_chunks.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/create_sample_chunks.sh) script to split test and control samples into chunks of sample IDs to parallely run chunks of samples in one nextflow job.

Nextflow pipeline requires two nextflow scripts [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/nextflow.config) and [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) , a tab delimitted [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) to read the parameters from and [Args_Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Args_Generate_FH_score.R) to submit parameters as arguments to the [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R). 

manifest.txt file can be easily generated using the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/R/Create_jobs.R)

3. Re-write the [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) and the [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R) as necessary.

4. Run all then nextflow jobs in the [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) file as below.
```bash
./run_nextflow.nf
```


Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
