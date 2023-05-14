#!/bin/bash
# * Run this once for one DCV* 
set -euxo pipefail
FoundHaplo_DIR=$1
DCV=$2

CHROMOSOME=$(echo "$DCV" | cut -d'.' -f2)
prefix="chr"
CHROMOSOME=${CHROMOSOME#"$prefix"}

echo "Finding start and end base pair positions to trim the VCF files."

Rscript $FoundHaplo_DIR/scripts/prepare_inputs/Run_Find_bp_to_trim.R $DCV $FoundHaplo_DIR/input_files/public_data/genetic_map_HapMapII_GRCh37 $FoundHaplo_DIR/temp/DCV_bp.txt
START_BP=$(cut -f2 $FoundHaplo_DIR/temp/DCV_bp.txt)
END_BP=$(cut -f3 $FoundHaplo_DIR/temp/DCV_bp.txt)

OUTPUT_NAME="$(echo $DCV | cut -d'.' -f1-2)"

echo "creating control cohorts for the disease variant" $DCV "for all five super populations in $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant."

vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AMR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/SAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AFR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR/$OUTPUT_NAME.vcf.gz
#####
