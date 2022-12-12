## Run FoundHaplo parallelly using Nextflow

Refer https://www.nextflow.io/docs/latest/process.html for more details

We recommend using a nextflow pipeline to run test/control cohort of 100 samples together in one job. i.e if there are 10000 samples in the test cohort nextflow will run 100 jobs.

This requires two nextflow scripts [Nextflow.config](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/nextflow.config) and [run_nextflow.nf](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/run_nextflow.nf) , a tab delimitted [manifest.txt](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/manifest.txt) to read the parameters from and [Args_Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Args_Generate_FH_score.R) to submit parameters as arguments to the [Generate_FH_score.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R)

manifest.txt file can be easily generated using the script [Create_jobs.R](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/run_nextflow/Create_jobs.R)


Go back to the [main page](https://github.com/bahlolab/FoundHaplo).
