#!/bin/bash
set -eu

FoundHaplo_PATH=$1 # path of FoundHaplo directory
INPUT_PLINK_PATH=$2 # example: FoundHaplo/example
INPUT_PLINK_BASE_NAME=$3 # example: FAME1_disease_cohort
CHROMOSOME=$4 # no "chr" prefix # example: 8

GENOTYPEHARMONIZER_PATH=$5
PLINK_TOOL_PATH=$6


module unload vcftools
module unload htslib
module load vcftools
module load htslib # to bgzip ## check if you can give a path to vcftools

echo "Harmonisation to 1000 Genomes data."

java -jar $GENOTYPEHARMONIZER_PATH \
--input $INPUT_PLINK_PATH/$INPUT_PLINK_BASE_NAME \
--inputType PLINK_BED \
--ref $FoundHaplo_PATH/input_files/public_data/1000G_plink/1000G_phase3_common_norel \
--refType PLINK_BED \
--update-id \
--update-reference-allele \
--output $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME \
--outputType PLINK_BED


echo "Removing SNPs with overall callrate<98%."

$PLINK_TOOL_PATH \
--bfile $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME \
--allow-no-sex \
--geno 0.02 \
--make-bed \
--out $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98

echo "Removing samples with callrate<98%."

$PLINK_TOOL_PATH \
--bfile $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98 \
--allow-no-sex \
--mind 0.02 \
--make-bed \
--out $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

echo "Converting to a VCF file in $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98"

$PLINK_TOOL_PATH --bfile $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98  --recode vcf --out $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

echo "Extracting chromsome $CHROMOSOME." 
echo "Saving to " "$FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz"
vcftools --vcf $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.vcf --chr $CHROMOSOME --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz


echo "impute" $FoundHaplo_PATH/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz "using Michigan server"
