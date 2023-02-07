Follow below steps everytime you either have a new test cohort to run or you need to test the same test cohort for a new disease-causing variant. You can skip step 1 if the new disease-causing variant is on a chromsome that is already being imputed.

1. Run [FoundHaplo/scripts/prepare_inputs/test_samples_pre_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/test_samples_pre_imputation.sh) to prepare the test cohort of SNP genotyping data.

```bash
FoundHaplo_PATH= 
INPUT_PLINK_PATH= 
INPUT_PLINK_BASE_NAME=
CHROMOSOME=
GENOTYPEHARMONIZER_PATH=
PLINK_TOOL_PATH= 
$FoundHaplo_PATH/scripts/prepare_inputs/test_samples_pre_imputation.sh $FoundHaplo_PATH $FoundHaplo_PATH/example FAME1_test_cohort 8 $GENOTYPEHARMONIZER_PATH $PLINK_TOOL_PATH
```
Set the variables as below,

* FoundHaplo_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* INPUT_PLINK_PATH :  Path to VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_PLINK_BASE_NAME : File name of the INPUT_VCF 
* CHROMOSOME: Chromosome relevant to the interested disease variant without "chr" prefix
* GENOTYPEHARMONIZER_PATH : Path to GenotypeHarmonizer.jar
* PLINK_TOOL_PATH : Path to Plink executable 

2. Impute the generated VCF file in FoundHaplo/temp/ using Michigan imputation server with 1000G phase 3 as reference panel and hg19 built. Do not filter for imutation quality yet. 

3. Create the test cohort with test individual haplotypes using the imputed VCF file generated in step 2 and [FoundHaplo/scripts/prepare_inputs/test_samples_post_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/test_samples_post_imputation.sh)

```bash
FoundHaplo_PATH= 
INPUT_VCF_PATH= 
INPUT_VCF_BASE_NAME=
DCV= 
$FoundHaplo_PATH/scripts/prepare_inputs/test_samples_post_imputation.sh $FoundHaplo_PATH $FoundHaplo_PATH/temp/FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz FAME1.chr8.119379052.
```
Set the variables as below,

* FoundHaplo_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* INPUT_VCF_PATH :  Path to imputed VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_VCF_BASE_NAME : File name of the imputed INPUT_VCF 
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052.

The test_samples_post_imputation.sh will create a VCF file with test individuals in FoundHaplo/input_files/input_vcf_data/

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
