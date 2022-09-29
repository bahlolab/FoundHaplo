#!/usr/bin/env bash

MAIN_PATH=$1 # path of FoundHaplo directory
INPUT_VCF_PATH=$2 # imputed vcf file path # example: FoundHaplo/temp/FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz
INPUT_VCF_BASE_NAME=$3 # example : FAME1_test_cohort.snp.0.98.sample.0.98.chr8
CHROMOSOME=$4 # no "chr" prefix # example : 8
DCV=$5 # example : "FAME1.chr8.119379052."

module unload vcftools
module unload bcftools
module load vcftools
module load bcftools


module unload R
module load R/4.2.0 # edit this line accordingly. load the R version with FoundHaplo


Rscript $MAIN_PATH/scripts/prepare_inputs/Find_bp_to_trim.R $DCV $MAIN_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37 $MAIN_PATH/temp/DCV_bp.txt
START_BP=$(cut -f2 $MAIN_PATH/temp/DCV_bp.txt)
END_BP=$(cut -f3 $MAIN_PATH/temp/DCV_bp.txt)


mkdir $MAIN_PATH/input_files/input_vcf_data
mkdir $MAIN_PATH/input_files/input_vcf_data/test_cohort

vcftools --gzvcf $INPUT_VCF_PATH --chr $CHROMOSOME --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/input_vcf_data/test_cohort/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz

mkdir $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples

bcftools query -l $MAIN_PATH/input_files/input_vcf_data/test_cohort/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz > $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples/samples.txt
cd $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples
split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples/samples.txt  $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples/file

rm -rf $MAIN_PATH/input_files/input_vcf_data/test_cohort/samples/samples.txt
rm -rf $MAIN_PATH/temp/*

