#!/bin/bash
set -euxo pipefail

FoundHaplo_DIR=$1 # FoundHaplo directory
INPUT_VCF_FILE=$2 # Imputed vcf file path # example: FoundHaplo/temp/FAME1_disease_cohort.snp.0.98.sample.0.98.chr8.vcf.gz
INPUT_VCF_BASE_NAME=$3 # example : FAME1_disease_cohort.snp.0.98.sample.0.98.chr8
DCV=$4 # example : "FAME1.chr8.119379052."
ANNOVAR_DIR=$5 # ANNOVAR directory
ANNOVAR_HUMANDB_DIR=$6 # ANNOVAR database directory

SAMPLE_INFO_FILE=$7 # Tab delimitted .txt file with sample names and type of phasing to be used included in a new line, include sample names as in the VCF file in mentioned order. # example : FoundHaplo/example/sample_info.txt

#' For the type "trio", affected-offspring,affected-parent,unaffected-parent trio
#' For the type "duo" or "related", affected-offspring,affected-parent,unaffected-parent duo

CHROMOSOME=$(echo "$DCV" | cut -d'.' -f2)
prefix="chr"
CHROMOSOME=${CHROMOSOME#"$prefix"}

module unload vcftools
module unload bcftools
module unload htslib
module load vcftools
module load bcftools
module load htslib

module unload R
module load R/4.2.0 # edit this line accordingly. load the R version with FoundHaplo [MB: --> this doesn't belong in your script because this is only certain specific HPC setups, people should manage modules themselves outside running any scripts you provide (and people who don't use module system or have different version numbers will get errors)]

echo "Finding start and end base pair positions to trim the VCF file."

Rscript $FoundHaplo_DIR/scripts/prepare_inputs/Run_Find_bp_to_trim.R $DCV $FoundHaplo_DIR/input_files/public_data/genetic_map_HapMapII_GRCh37 $FoundHaplo_DIR/temp/DCV_bp.txt
START_BP=$(cut -f2 $FoundHaplo_DIR/temp/DCV_bp.txt)
END_BP=$(cut -f3 $FoundHaplo_DIR/temp/DCV_bp.txt)

vcftools --gzvcf $INPUT_VCF_FILE --chr $CHROMOSOME --remove-indels --min-alleles 2 --max-alleles 2 --from-bp $START_BP --to-bp $END_BP --recode --recode-INFO-all --stdout | bgzip -c > $FoundHaplo_DIR/temp/$INPUT_VCF_BASE_NAME.imputed.trimmed.vcf.gz

ANNOVAR_SCRIPT=$ANNOVAR_DIR/table_annovar.pl

ANNOVAR_OUTPUT_FILENAME_BASE=$FoundHaplo_DIR/temp/$INPUT_VCF_BASE_NAME.imputed.trimmed # name of the annotated file

echo "Annotating gnomAD population allele frequencies."

# annotate population frequencies form gnomAD
perl $ANNOVAR_SCRIPT $ANNOVAR_OUTPUT_FILENAME_BASE.vcf.gz \
$ANNOVAR_HUMANDB_DIR -buildver hg19 \
-vcfinput -out $ANNOVAR_OUTPUT_FILENAME_BASE -remove \
-protocol gnomad211_genome -operation f -nastring .


# $OUTPUT_FILENAME_BASE.hg19_multianno.vcf has annotated INFO based on gnomAD

# Note : Example 1000 Genomes WGS data does not have the R2 (Imputation quality) tag 

n_lines_R2=$(cat "$ANNOVAR_OUTPUT_FILENAME_BASE".hg19_multianno.vcf | tail -1 | grep "R2=" |wc -l)

if [ $n_lines_R2 -eq 0 ] 

then  
    bcftools annotate -x ^INFO/AF_raw,^INFO/AF_afr,^INFO/AF_sas,^INFO/AF_amr,^INFO/AF_eas,^INFO/AF_nfe,^INFO/AF_fin "$ANNOVAR_OUTPUT_FILENAME_BASE".hg19_multianno.vcf > $FoundHaplo_DIR/temp/ready.to.phase.vcf
else
    bcftools annotate -x ^INFO/R2,^INFO/AF_raw,^INFO/AF_afr,^INFO/AF_sas,^INFO/AF_amr,^INFO/AF_eas,^INFO/AF_nfe,^INFO/AF_fin "$ANNOVAR_OUTPUT_FILENAME_BASE".hg19_multianno.vcf > $FoundHaplo_DIR/temp/ready.to.phase.vcf
fi 

# phase by pedigrees
# Rscript will create seperate VCF files with known disease haplotypes

echo "Phasing disease haplotypes by pedigree information."

mkdir -p $FoundHaplo_DIR/input_files/input_vcf_data/disease_haplotypes

Rscript $FoundHaplo_DIR/scripts/prepare_inputs/Run_Phasing_by_pedigree.R $FoundHaplo_DIR/temp/ready.to.phase.vcf $FoundHapFoundHaplo_DIRlo_PATH/input_files/input_vcf_data/disease_haplotypes $SAMPLE_INFO_FILE

echo "bgzipping all the VCF files with disease haplotypes and saving to $FoundHaplo_DIR/input_files/input_vcf_data/disease_haplotypes."
ls $FoundHaplo_DIR/input_files/input_vcf_data/disease_haplotypes/*.vcf | xargs -n1 bgzip
