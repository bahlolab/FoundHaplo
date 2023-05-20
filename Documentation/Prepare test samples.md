Follow below steps everytime you either have a new test cohort to run or you need to test the same test cohort for a new disease-causing variant. You can skip steps 1 and 2 if the new disease-causing variant is on a chromsome that is already being imputed.

1. Run [FoundHaplo/scripts/prepare_inputs/test_samples_pre_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/test_samples_pre_imputation.sh) to prepare the test cohort of SNP genotyping data.

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
INPUT_PLINK_DIR=/mypath/PLINK_file_DIR
INPUT_PLINK_BASE_NAME=PLINK_file
CHROMOSOME=CHROMOSOME
GENOTYPEHARMONIZER_JAR=/mypath/GenotypeHarmonizer.jar
PLINK_TOOL_EXECUTABLE=/mypath/PLINK_executable 

$FoundHaplo_DIR/scripts/prepare_inputs/test_samples_pre_imputation.sh $FoundHaplo_DIR $INPUT_PLINK_DIR INPUT_PLINK_BASE_NAME $CHROMOSOME $GENOTYPEHARMONIZER_JAR $PLINK_TOOL_EXECUTABLE
```
Run example below,

```bash
$FoundHaplo_DIR/scripts/prepare_inputs/test_samples_pre_imputation.sh $FoundHaplo_DIR $FoundHaplo_DIR/example FAME1_test_cohort 8 $GENOTYPEHARMONIZER_JAR $PLINK_TOOL_EXECUTABLE
```
Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e /mypath/FoundHaplo
* INPUT_PLINK_DIR : Directory to the PLINK file that has individuals to be tested for known disease variants
* INPUT_PLINK_BASE_NAME : File name of the PLINK file 
* CHROMOSOME: Chromosome relevant to the interested disease variant without "chr" prefix
* GENOTYPEHARMONIZER_JAR : Path to GenotypeHarmonizer.jar
* PLINK_TOOL_EXECUTABLE : Path to Plink executable 

2. Impute the generated VCF file in FoundHaplo/temp/ using Michigan imputation server (Genotype Imputaion Minimac 4) with 1000G phase 3 V5 (GRCh37/hg19) as reference panel with the hg19 human genome build. Use below parameters when using [Michigan imputation server](https://imputationserver.sph.umich.edu/). 

Array build GRCh37/hg19.

Do not filter for imputation quality yet (keep rsq Filter off). 

Phase using Eagle v2.4.

Select the relevant ancestral population for the cohort (FoundHaplo example is of EUR ancestry).

Mode "Quality Control and Phasing Only".

Rename and save the resulting imputed and phased file with the "imputed_phased_" prefix, and its original file name in the same location.

3. Create the test cohort with test individual haplotypes using the imputed VCF file generated in step 2 and [FoundHaplo/scripts/prepare_inputs/test_samples_post_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/test_samples_post_imputation.sh)

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
INPUT_VCF_FILE=/mypath/VCF_file  
INPUT_VCF_BASE_NAME=VCF_file
DCV=DCV

$FoundHaplo_DIR/scripts/prepare_inputs/test_samples_post_imputation.sh $FoundHaplo_DIR $INPUT_VCF_FILE $INPUT_VCF_BASE_NAME $DCV
```
Run example below,

```bash
$FoundHaplo_DIR/scripts/prepare_inputs/test_samples_post_imputation.sh $FoundHaplo_DIR $FoundHaplo_DIR/temp/imputed_phased_FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz imputed_phased_FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz FAME1.chr8.119379052.
```
Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e path/FoundHaplo
* INPUT_VCF_FILE : Imputed VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_VCF_BASE_NAME : File name of the imputed INPUT_VCF 
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052.

The test_samples_post_imputation.sh will create a VCF file with test individuals in FoundHaplo/input_files/input_vcf_data/test_cohort

Go back to the [documentation](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).
