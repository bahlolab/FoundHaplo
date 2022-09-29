# Example

This folder contains an example of running FoundHaplo for two disease variants i) FAME1 variant on chromosome 8, SAMD12 that causes Familial Adult Myclonic Epilepsy ii) GEFS+ variant on chromosome 19, SCN1B that causes General Epilepsy with Febrile Seizures Plus type 1. 

Example data are simulated using 1000 Genomes WGS haplotypes

*Test VCF files are predominantely European and are in /Example/Cohort_A/FAME1/ and /Example/Cohort_A/SCN1B_GEFS+/

*Disease sample VCF files are in /Example/Derived_disease_haplotypes/FAME1/Input_disease_VCF and /Example/Derived_disease_haplotypes/SCN1B_GEFS+/Input_disease_VCF
    Phase the disease sample VCF files using [R/Phasing_by_pedigree.R](https://github.com/bahlolab/FoundHaplo/blob/main/R/Phasing_by_pedigree.R), and save the output in /Example/Derived_disease_haplotypes/FAME1/Disease_haplotypes and /Example/Derived_disease_haplotypes/SCN1B_GEFS+/Disease_haplotypes 

*Control cohort VCF files are in /Example/Control_cohort/EUR/ for both disease variants. 

*Install FoundHaplo with devtools::install_github("bahlolab/FoundHaplo") and run nextflow pipeline in Example/Cohort_A/FAME1/Run_FH_Nextflow/ and Example/Cohort_A/SCN1B_GEFS+/Run_FH_Nextflow/

nextflow - 
nextflow.config - Configuration file
manifest.txt - Tab delimitted .txt file to read the input parameters from listing [Parameters that must be specified by the user](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Parameters%20in%20the%20algorithm.md)
create_jobs.R - R script to easily generate the manifest.txt file
Args_Generate_FH_score.R - R script to run arguments in the manifest.txt files on the [Generate_FH_score](https://github.com/bahlolab/FoundHaplo/blob/main/R/Generate_FH_score.R)
run_nextflow.nf - Main nextflow file that parallely runs manifest.txt

```bash
module load nextflow
nohup ./run_nextflow.nf 
```bash

