#!/bin/bash
# * Run this once for one DCV* 
set -eu
FoundHaplo_PATH=$1
DCV=$2

CHROMOSOME=$(echo "$DCV" | cut -d'.' -f2)
prefix="chr"
CHROMOSOME=${CHROMOSOME#"$prefix"}

module unload R
module load R/4.2.0 # edit this line accordingly. load the R version with FoundHaplo [--> this doesn't belong in your script because this is only certain specific HPC setups, people should manage modules themselves outside running any scripts you provide (and people who don't use module system or have different version numbers will get errors)]

Rscript $FoundHaplo_PATH/scripts/prepare_inputs/Find_bp_to_trim.R $DCV $FoundHaplo_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37 $FoundHaplo_PATH/temp/DCV_bp.txt
START_BP=$(cut -f2 $FoundHaplo_PATH/temp/DCV_bp.txt)
END_BP=$(cut -f3 $FoundHaplo_PATH/temp/DCV_bp.txt)

OUTPUT_NAME="$(echo $DCV | cut -d'.' -f1-2)"

module load vcftools
module load htslib
# create control cohorts for the disease variant for all five super populations.

vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AMR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/SAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AFR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR/$OUTPUT_NAME.vcf.gz
#####
