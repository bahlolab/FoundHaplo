# * Run this once for one DCV* 

MAIN_PATH=$1
DCV=$2
CHROMOSOME=$3
set -eu
module unload R
module load R/4.2.0 # edit this line accordingly. load the R version with FoundHaplo

Rscript $MAIN_PATH/scripts/prepare_inputs/Find_bp_to_trim.R $DCV $MAIN_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37 $MAIN_PATH/temp/DCV_bp.txt
START_BP=$(cut -f2 $MAIN_PATH/temp/DCV_bp.txt)
END_BP=$(cut -f3 $MAIN_PATH/temp/DCV_bp.txt)

OUTPUT_NAME="$(echo $DCV | cut -d'.' -f1-2)"

module load vcftools
module load htslib
# create control cohorts for the disease variant for all five super populations.

vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AMR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/SAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS/$OUTPUT_NAME.vcf.gz

vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EAS.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS/$OUTPUT_NAME.vcf.gz


vcftools --gzvcf $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/ALL.chr$CHROMOSOME.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --chr $CHROMOSOME --keep $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AFR.txt --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR/$OUTPUT_NAME.vcf.gz
#Add as many disease variants as you want
#####
