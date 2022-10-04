#!/bin/bash
set -eu

MAIN_PATH=$1 # path of FoundHaplo directory
INPUT_VCF_FILE=$2 # imputed vcf file path # example: FoundHaplo/temp/FAME1_disease_cohort.snp.0.98.sample.0.98.chr8.vcf.gz
INPUT_VCF_BASE_NAME=$3 # example : FAME1_disease_cohort.snp.0.98.sample.0.98.chr8
CHROMOSOME=$4 # no "chr" prefix # example : 8
DCV=$5 # example : "FAME1.chr8.119379052."
ANNOVAR_PATH=$6 # path o ANNOVAR directory
ANNOVAR_HUMANDB_DIR_PATH=$7 # path to ANNOVAR databases

SAMPLE_INFO_FILE=$8 # Path to a tab delimitted .txt file with sample names and type of phasing to be used included in a new line, include sample names as in the VCF file in mentioned order. # example : FoundHaplo/sample_info.txt

#' For the type "trio", affected-offspring,affected-parent,unaffected-parent trio
#' For the type "duo" or "related", affected-offspring,affected-parent,unaffected-parent duo


module unload vcftools
module unload bcftools
module unload htslib
module load vcftools
module load bcftools
module load htslib

module unload R
module load R/4.2.0 # edit this line accordingly. load the R version with FoundHaplo

Rscript $MAIN_PATH/scripts/prepare_inputs/Run_Find_bp_to_trim.R $DCV $MAIN_PATH/input_files/public_data/genetic_map_HapMapII_GRCh37 $MAIN_PATH/temp/DCV_bp.txt
START_BP=$(cut -f2 $MAIN_PATH/temp/DCV_bp.txt)
END_BP=$(cut -f3 $MAIN_PATH/temp/DCV_bp.txt)

vcftools --gzvcf $INPUT_VCF_FILE --chr $CHROMOSOME --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $MAIN_PATH/temp/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz

ANNOVAR_SCRIPT=$ANNOVAR_PATH/table_annovar.pl

ANNOVAR_OUTPUT_FILENAME_BASE=$MAIN_PATH/temp/$INPUT_VCF_BASE_NAME.imputed.trimmed # name of the annotated file

# annotate population frequencies form gnomAD
perl $ANNOVAR_SCRIPT $ANNOVAR_OUTPUT_FILENAME_BASE.vcf.gz \
$ANNOVAR_HUMANDB_DIR_PATH -buildver hg19 \
-vcfinput -out $ANNOVAR_OUTPUT_FILENAME_BASE -remove \
-protocol gnomad211_genome -operation f -nastring .


# $OUTPUT_FILENAME_BASE.hg19_multianno.vcf has annotated INFO based on gnomAD
bcftools annotate -x ^INFO/R2,^INFO/AF_raw,^INFO/AF_afr,^INFO/AF_sas,^INFO/AF_amr,^INFO/AF_eas,^INFO/AF_nfe,^INFO/AF_fin "$ANNOVAR_OUTPUT_FILENAME_BASE".hg19_multianno.vcf > $MAIN_PATH/temp/ready.to.phase.vcf

# phase by pedigrees


# Rscript will create seperate VCF files with known disease haplotypes

mkdir -p $MAIN_PATH/input_files/input_vcf_data/disease_haplotypes

Rscript $MAIN_PATH/scripts/prepare_inputs/Run_Phasing_by_pedigree.R $MAIN_PATH/temp/ready.to.phase.vcf $MAIN_PATH/input_files/input_vcf_data/disease_haplotypes $SAMPLE_INFO_FILE


ls $MAIN_PATH/input_files/input_vcf_data/disease_haplotypes/*.vcf | xargs -n1 bgzip
