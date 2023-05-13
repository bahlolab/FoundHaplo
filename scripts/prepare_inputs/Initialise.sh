#!/usr/bin/env bash
set -euxo pipefail

echo "Make sure you have installed FoundHaplo in R"
FoundHaplo_PATH=$1

#Publicly available reference files required to run FoundHaplo

echo "This script downloads, 1) 1000 Genomes phased 3 bed/bim/fam files to harmonize input VCF data and 2) 1000 Genomes phase 3 haplotypes are used as the control cohort when running FoundHaplo."
echo "1000 Genomes sample names by super populations are already in $FoundHaplo_PATH/input_files/public_data/1000G_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt, EUR.txt, AMR.txt, EAS.txt, SAS.txt and AFR.txt."
echo "Hapmap recombination files in hg19 to calculate the recombination rates are in $FoundHaplo_PATH/input_files."

echo "Creating results and temp directories"

mkdir -p $FoundHaplo_PATH/results/FH_IBD_scores # to save the final IBD report in .txt file
mkdir -p $FoundHaplo_PATH/results/FH_Analysis # to save results after analysing FH scores in the final IBD report
mkdir -p $FoundHaplo_PATH/temp # to save temporary files

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_plink


echo "Downloading 1000 Genomes phased 3 bed/bim/fam files to $FoundHaplo_PATH/input_files/public_data/1000G_plink."

cd $FoundHaplo_PATH/input_files/public_data/1000G_plink
wget https://ndownloader.figshare.com/files/17838962 --output-document "1000G_plink.zip"
unzip 1000G_plink.zip

echo "Downloading 1000 Genomes phase 3 haplotypes to $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original."

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original
cd $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original

wget --recursive --no-parent http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/

echo "main files are downloaded"

mv $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/* $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original

rm -r $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_original/hgdownload.cse.ucsc.edu

echo "creating control cohorts for the disease variant for all five super populations in $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant."

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant

mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS
mkdir -p $FoundHaplo_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR

echo "Download complete for all required input data to run FoundHaplo"

echo "Make sure you have all software tools listed in modules or please download them : List of software tools required are GenotypeHarmonizer,plink,ANNOVAR,vcftools and bcftools."

echo "Download all software tools listed below"

echo "1) GenotypeHarmonizer tool is used to harmomize input VCF data to 1000Genomes."
echo "2) Plink 1.9 is used to perform quality control steps on input VCF data."
echo "3) ANNOVAR tool is used to annotate population frequencies form gnomAD data."
echo "4) VCFtools version v0.1.13 and BCFtools version 1.16 to perform queries on VCF files."
echo "5) Install Nextflow to parallely run FoundHaplo (Optional)."

