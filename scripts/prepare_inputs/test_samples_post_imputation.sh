#!/bin/bash
set -euxo pipefail

FoundHaplo_DIR=$1 # path of FoundHaplo directory
INPUT_VCF_FILE=$2 # imputed vcf file path # example: FoundHaplo/temp/FAME1_test_cohort.snp.0.98.sample.0.98.chr8.vcf.gz
INPUT_VCF_BASE_NAME=$3 # example : FAME1_test_cohort.snp.0.98.sample.0.98.chr8
DCV=$4 # example : "FAME1.chr8.119379052."

CHROMOSOME=$(echo "$DCV" | cut -d'.' -f2)
prefix="chr"
CHROMOSOME=${CHROMOSOME#"$prefix"}

echo "Finding start and end base pair positions to trim the VCF file."

Rscript $FoundHaplo_DIR/scripts/prepare_inputs/Run_Find_bp_to_trim.R $DCV $FoundHaplo_DIR/input_files/public_data/genetic_map_HapMapII_GRCh37 $FoundHaplo_DIR/temp/DCV_bp.txt
START_BP=$(cut -f2 $FoundHaplo_DIR/temp/DCV_bp.txt)
END_BP=$(cut -f3 $FoundHaplo_DIR/temp/DCV_bp.txt)


mkdir -p $FoundHaplo_DIR/input_files/input_vcf_data
mkdir -p $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort

echo "Saving to $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz." 

vcftools --gzvcf $INPUT_VCF_FILE --chr $CHROMOSOME --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz

mkdir -p $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples

echo "Saving sample IDs to $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples/samples.txt"

bcftools query -l $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz > $FoundHaplo_DIR/input_files/input_vcf_data/test_cohort/samples/samples.txt

rm -rf $FoundHaplo_DIR/temp/*

