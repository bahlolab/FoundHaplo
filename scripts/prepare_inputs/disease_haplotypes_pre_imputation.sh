#!/bin/bash
set -euxo pipefail

FoundHaplo_DIR=$1 # FoundHaplo directory
INPUT_PLINK_FILE=$2 # example: FoundHaplo/example
INPUT_PLINK_BASE_NAME=$3 # example: FAME1_disease_cohort
CHROMOSOME=$4 # no "chr" prefix # example: 8

GENOTYPEHARMONIZER_JAR=$5 # GENOTYPEHARMONIZER_JAR file
PLINK_TOOL_EXECUTABLE=$6 # /mypath/plink2


module unload vcftools
module unload htslib
module load vcftools
module load htslib # to bgzip ## check if you can give a path to vcftools
  
echo "Harmonisation to 1000 Genomes data."

java -jar $GENOTYPEHARMONIZER_JAR \
--input $INPUT_PLINK_FILE/$INPUT_PLINK_BASE_NAME \
--inputType PLINK_BED \
--ref $FoundHaplo_DIR/input_files/public_data/1000G_plink/1000G_phase3_common_norel \
--refType PLINK_BED \
--update-id \
--update-reference-allele \
--output $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME \
--outputType PLINK_BED


echo "Removing SNPs with overall callrate<98%."

$PLINK_TOOL_EXECUTABLE \
--bfile $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME \
--allow-no-sex \
--geno 0.02 \
--make-bed \
--out $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98

echo "Removing samples with callrate<98%."

$PLINK_TOOL_EXECUTABLE \
--bfile $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98 \
--allow-no-sex \
--mind 0.02 \
--make-bed \
--out $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

echo "Converting to a VCF file in $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98"

$PLINK_TOOL_EXECUTABLE --bfile $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98  --recode vcf --out $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98

echo "Extracting chromsome $CHROMOSOME." 
echo "Saving to " "$FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz"
vcftools --vcf $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.vcf --chr $CHROMOSOME --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz


echo "impute" $FoundHaplo_DIR/temp/$INPUT_PLINK_BASE_NAME.snp.0.98.sample.0.98.chr$CHROMOSOME.vcf.gz "using Michigan server"
