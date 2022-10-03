#!/bin/bash
set -eu

MAIN_PATH=$1 # path of FoundHaplo directory
INPUT_PLINK_PATH=$2 # example: FoundHaplo/example
INPUT_PLINK_BASE_NAME=$3 # example: FAME1_disease_cohort
CHROMOSOME=$4 # no "chr" prefix # example: 8

GENOTYPEHARMONIZER_PATH=$5
PLINK_TOOL_PATH=$6


module unload vcftools
module unload htslib
module load vcftools
module load htslib # to bgzip ## check if you can give a path to vcftools

#Harmonisation to 1000G

java -jar $GENOTYPEHARMONIZER_PATH/GenotypeHarmonizer.jar \
--input $INPUT_PLINK_PATH/$INPUT_PLINK_BASE_NAME \
--inputType PLINK_BED \
--ref $MAIN_PATH/input_files/public_data/1000G_plink/1000G_phase3_common_norel \
--refType PLINK_BED \
--update-id \
--update-reference-allele \
--output $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME \
--outputType PLINK_BED


# remove SNPs with overall callrate<98

$PLINK_TOOL_PATH \
--bfile $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME \
--allow-no-sex \
--geno 0.02 \
--make-bed \
--out $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98


# remove samples with callrate<98%
$PLINK_TOOL_PATH \
--bfile $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98 \
--allow-no-sex \
--mind 0.02 \
--make-bed \
--out $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

$PLINK_TOOL_PATH --bfile $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98  --recode vcf --out $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

vcftools --vcf $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.vcf --chr $CHROMOSOME --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz


echo "impute" $INPUT_PLINK_BASE_NAME_snp_0.98_sample_0.98.chr$CHROMOSOME.vcf.gz "using Michigan server"
