Follow below steps to create disease haplotypes from known affected samples using SNP chip data.

FH model is sensitive to phasing errors on the disease samples, hence it requires accurately formed disease haplotypes. This can be achieved by phasing the affected individuals by trio , duo or any type of pedigree phasing to derive the disease haplotypes shared among affected individuals instead of using traditional phasing tools. 

**You can proceed if you have any of below**.

* Genotype data of offsprings known to have a disease varaint and both parents, where you know which parent is affected.
* Genotype data of any two related individuals who share a disease variant.

Can merge as many samples from differnt families into a single VCF file as long as all are generated from the same SNPchip design.

1. Run FoundHaplo/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh.

```bash
FoundHaplo/scripts/prepare_inputs/disease_haplotypes_pre_imputation.sh MAIN_PATH INPUT_VCF_PATH INPUT_VCF_BASE_NAME CHROMOSOME GENOTYPEHARMONIZER_PATH PLINK_PATH VCFTOOLS_PATH
```
Arguments are explained below,

* MAIN_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* INPUT_VCF_PATH :  Path to VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_VCF_BASE_NAME : File name of the INPUT_VCF 
* CHROMOSOME: Chromosome relevant to the interested disease variant without "chr" prefix
* GENOTYPEHARMONIZER_PATH : Path to GenotypeHarmonizer.jar
* PLINK_PATH : Path to Plink executable 
* VCFTOOLS_PATH : Path to vcftoools executable 

2. Impute the generated VCF file in MAIN_PATH/temp/ using Michigan imputation server with 1000G phase 3 as refereance panel and hg19 built. Do not filter for imutation quality yet. 

3. Create disease haplotypes using the imputed VCF file generated in step 2 and FoundHaplo/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh. 
```bash
FoundHaplo/scripts/prepare_inputs/disease_haplotypes_post_imputation.sh MAIN_PATH INPUT_VCF_PATH INPUT_VCF_BASE_NAME CHROMOSOME DCV VCFTOOLS_PATH ANNOVAR_PATH ANNOVAR_HUMANDB_DIR_PATH BCFTOOLS_PATH TYPE SAMPLE_NAMES
```

Arguments are explained below,

* MAIN_PATH : Path to FoundHaplo directory i.e path/FoundHaplo
* INPUT_VCF_PATH :  Path to imputed VCF file with individuals with known diease variants that can be pedigree phased.
* INPUT_VCF_BASE_NAME : File name of the imputed INPUT_VCF 
* CHROMOSOME : Chromosome relevant to the interested disease variant without "chr" prefix
* DCV : Name the disease variant of interest in the format of disease.chr.position. i.e FAME1.chr8.119379052.
* VCFTOOLS_PATH : Path to vcftoools executable 
* ANNOVAR_PATH : Path to ANNOVAR directory
* ANNOVAR_HUMANDB_DIR_PATH : Path to ANNOVAR database
* BCFTOOLS_PATH : Path to bcftoools executable
* TYPE : Type of phasing , "trio" or "duo"
* SAMPLE_NAMES : Give sample names to phase in specific order,
               If TYPE is "trio", SAMPLE_NAMES="affected offspring" "affected parent" "unaffected parent"
               If TYPE is "duo", SAMPLE_NAMES="affected individual 1" "affected individual 2"

The disease_haplotypes_post_imputation.sh will create a seperate file with VCF columns for each derived disease haplotype in FoundHaplo/input_files/input_vcf_data/disease_haplotypes/, additionaly it will remove multiallelic SNPs using gnomAD frequency files downlaoded with ANNOVAR. 




