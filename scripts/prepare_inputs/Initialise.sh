#!/usr/bin/env bash

# Make sure you have installed FoundHaplo in R
FoundHaplo_PATH=$1

#Publicly available reference files required to run FoundHaplo

# 1) 1000 Genomes phased 3 bed/bim/fam files to harmonize input VCF data.
# 2) 1000 Genomes phase 3 haplotypes are used as the control cohort when running FoundHaplo.
# 3) 1000G sample names by super populations
# 4) Hapmap recombination files in hg19 to calculate the recombination rates (in /input_files).

# download 1000 Genomes phased 3 bed/bim/fam files

mkdir -p $FoundHaplo_PATH/results # to save final IBD report in .txt file

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_plink

cd $FoundHaplo_PATH/input_files/public_data/1000G_plink
wget https://ndownloader.figshare.com/files/17838962 --output-document "1000G_plink.zip"
unzip 1000G_plink.zip

# download 1000 Genomes phase 3 haplotypes


mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original
cd $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original

wget --recursive --no-parent http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/

mv $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/* $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original

rm -r $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/hgdownload.cse.ucsc.edu

# 1000 Genomes sample names are in /FoundHaplo/input_files/public_data/1000G_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt, EUR.txt, AMR.txt, EAS.txt, SAS.txt and AFR.txt.

# create control cohorts for the disease variant for all five super populations.

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR

# 1000Genome sample IDs are in $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population

# Downlaod all software tools listed

# 1) GenotypeHarmonizer tool is used to harmomize input VCF data to 1000Genomes.
# 2) Plink is used to perform quality control steps on input VCF data.
# 3) ANNOVAR tool is used to annotate population frequencies form gnomAD data.
# 4) vcftools and bcftools to perform queries on VCF files

