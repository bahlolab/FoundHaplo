Follow below steps eveytime when new disease haplotypes are available to create disease haplotypes from known affected samples using SNP genotyping data.

FH model is sensitive to phasing errors on the disease samples, hence it requires accurately formed disease haplotypes. This can be achieved by phasing the affected individuals by trio , duo or any type of pedigree phasing to derive the disease haplotypes shared among affected individuals instead of using traditional phasing tools. 

**You can proceed if you have any of below**.

* Genotype data of offsprings known to have a disease varaint and both parents, where you know which parent is affected.
* Genotype data of any two related individuals who share a disease variant.

Can merge as many samples from different families into a single VCF file as long as all are generated from the same SNPchip design.

1. Run [FoundHaplo/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh).

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
INPUT_PLINK_DIR=/mypath/PLINK_file_DIR
INPUT_PLINK_BASE_NAME=PLINK_file
CHROMOSOME=CHROMOSOME
GENOTYPEHARMONIZER_JAR=/mypath/GenotypeHarmonizer.jar
PLINK_TOOL_EXECUTABLE=/mypath/PLINK_executable

$FoundHaplo_DIR/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh $FoundHaplo_DIR $INPUT_PLINK_DIR $INPUT_PLINK_BASE_NAME $CHROMOSOME $GENOTYPEHARMONIZER_JAR $PLINK_TOOL_EXECUTABLE
```

Run example below,
```bash
$FoundHaplo_DIR/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh $FoundHaplo_DIR $FoundHaplo_DIR/example FAME1_disease_cohort 8 $GENOTYPEHARMONIZER_JAR $PLINK_TOOL_EXECUTABLE
```

Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e /mypath/FoundHaplo
* INPUT_PLINK_DIR : Directory to the PLINK file that has individuals with known diease variants that can be pedigree phased.
* INPUT_PLINK_BASE_NAME : File name of the PLINK file 
* CHROMOSOME: Chromosome relevant to the interested disease variant without "chr" prefix
* GENOTYPEHARMONIZER_JAR : Path to GenotypeHarmonizer.jar
* PLINK_TOOL_EXECUTABLE : Path to Plink executable 

2. Impute the generated VCF file in FoundHaplo/temp/ using Michigan imputation server with 1000G phase 3 V5 (GRCh37/hg19) as reference panel with the hg19 human genome build. Use below parameters when using [Michigan imputation server](https://imputationserver.sph.umich.edu/). 

Do not filter for imputation quality yet (keep rsq Filter off). 

Array build GRCh37/hg19.

Phase using Eagle v2.4

Mode "Quality Control and Phasing Only"

Rename and save the resulting imputed and phased file with the "imputed_phased_" prefix, and its original file name in the same location.

3. Create disease haplotypes using the imputed VCF file generated in step 2 and [FoundHaplo/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh](https://github.com/bahlolab/FoundHaplo/blob/main/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh) 

```bash
FoundHaplo_DIR=/mypath/FoundHaplo
INPUT_VCF_FILE=/mypath/VCF_file 
INPUT_VCF_BASE_NAME=VCF_file
DCV=DCV
ANNOVAR_DIR=/mypath/annovar 
ANNOVAR_HUMANDB_DIR=/mypath/annovar/humandb 
SAMPLE_INFO_FILE=/mypath/sample_info.txt 

$FoundHaplo_DIR/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh $FoundHaplo_DIR $INPUT_VCF_FILE $INPUT_VCF_BASE_NAME $DCV $ANNOVAR_DIR $ANNOVAR_HUMANDB_DIR $SAMPLE_INFO_FILE
```
Run example below,
```bash
$FoundHaplo_DIR/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh $FoundHaplo_DIR $FoundHaplo_DIR/temp/imputed_phased_FAME1_disease_cohort.snp.0.98.sample.0.98.chr8.vcf.gz imputed_phased_FAME1_disease_cohort.snp.0.98.sample.0.98.chr8.vcf.gz FAME1.chr8.119379052. $ANNOVAR_DIR $ANNOVAR_HUMANDB_DIR $FoundHaplo_DIR/example/sample_info.txt
```

Set the variables as below,

* FoundHaplo_DIR : FoundHaplo directory i.e /mypath/FoundHaplo
* INPUT_VCF_FILE : Imputed VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_VCF_BASE_NAME : File name of the imputed INPUT_VCF 
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052.
* ANNOVAR_DIR : ANNOVAR directory
* ANNOVAR_HUMANDB_DIR : ANNOVAR database directory
* SAMPLE_INFO_FILE : Path to a tab delimitted .txt file with sample names and type of phasing to be used included in a new line, include sample names as in the VCF file in mentioned order.   
For the type "trio", affected-offspring,affected-parent,unaffected-parent trio  
For the type "duo", affected-offspring,affected-parent duo
For the type "related", affected-sample-1,affected-sample2 related

The post imputation script will create a seperate file with VCF columns for each derived disease haplotype in FoundHaplo/input_files/input_vcf_data/disease_haplotypes/, additionaly it will remove multiallelic SNPs using gnomAD frequency files downloaded with ANNOVAR. 

The created disease haplotypes in FoundHaplo/input_files/input_vcf_data/disease_haplotypes/ can be efficiently managed in a MySQL database system. How to create and manage the FoundHaplo database is explained [here](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Prepare%20database%20with%20known%20disease%20haplotypes.md). 

Go back to the [documentaton](https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md).



