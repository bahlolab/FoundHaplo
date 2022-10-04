# Example

This folder contains an example nextflow scripts of running FoundHaplo for FAME1 variant on chromosome 8, SAMD12 that causes Familial Adult Myclonic Epilepsy 

Example data are simulated using 1000 Genomes WGS haplotypes

*Test VCF files are predominantely European and must be created in  /FoundHaplo/input_files/input_vcf_data/test_cohort/  using /FoundHaplo/scripts/prepare_inputs/

*Disease sample VCF files must be created in /FoundHaplo/input_files/input_vcf_data/disease_haplotypes/ using /FoundHaplo/scripts/prepare_inputs/

*Control cohort VCF files must be created in /FoundHaplo/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant

*Install FoundHaplo with devtools::install_github("bahlolab/FoundHaplo") and run nextflow pipeline using /FoundHaplo/scripts/prepare_inputs/run_nextflow/run_nextflow.nf
Edit all teh scripts as necesssary. 

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

