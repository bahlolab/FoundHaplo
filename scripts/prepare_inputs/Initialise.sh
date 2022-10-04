#!/usr/bin/env bash

# Make sure you have installed FoundHaplo in R
#### Do not change the code in between
MAIN_PATH=$1
#Publicly available reference files required to run FoundHaplo

# 1) 1000 Genomes phased 3 bed/bim/fam files to harmonize input VCF data.
# 2) 1000 Genomes phase 3 haplotypes are used as the control cohort when running FoundHaplo.
# 3) 1000G sample names by super populations
# 4) Hapmap recombination files in hg19 to calculate the recombination rates (in /input_files).

# download 1000 Genomes phased 3 bed/bim/fam files

mkdir $MAIN_PATH/results # to save final IBD report in .txt file

mkdir $MAIN_PATH/public_data
mkdir $MAIN_PATH/public_data/1000G_plink

cd $MAIN_PATH/public_data/1000G_plink
wget https://ndownloader.figshare.com/files/17838962 --output-document "1000G_plink.zip"
unzip 1000G_plink.zip

# download 1000 Genomes phase 3 haplotypes


mkdir $MAIN_PATH/input_files/public_data/1000G_haplotypes/1000G_haplotypes_original
cd $MAIN_PATH/input_files/public_data/1000G_haplotypes/1000G_haplotypes_original

wget --recursive --no-parent http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/

# 1000 Genomes sample names are in /FoundHaplo/input_files/public_data/1000G_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt, EUR.txt, AMR.txt, EAS.txt, SAS.txt and AFR.txt.

# create control cohorts for the disease variant for all five super populations.

mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant

mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR

mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL/samples
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR/samples
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AMR/samples
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS/samples
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS/samples
mkdir $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR/samples

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/ALL.txt  $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/ALL/samples/file

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EUR.txt  $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EUR/samples/file

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypesG_haplotypes/1000G_haplotypes_samples_by_population/AMR.txt  $MAIN_PATH/input_files/public_data/1000G_haplotypes/1000G_control_haplotypes/AMR/samples/file

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/SAS.txt  $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/SAS/samples/file

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/EAS.txt  $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/EAS/samples/file

split -l 100 -d --additional-suffix=.txt $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_samples_by_population/AFR.txt  $MAIN_PATH/input_files/public_data/1000G_control_haplotypes/1000G_haplotypes_by_variant/AFR/samples/file

#### Do not change the code in between

# change paths of software tools as necessary

# 1) GenotypeHarmonizer tool is used to harmomize input VCF data to 1000Genomes.
# 2) Plink is used to perform quality control steps on input VCF data.
# 3) ANNOVAR tool is used to annotate population frequencies form gnomAD data.
# 4) vcftools and bcftools to perform queries on VCF files


#Publicly available software tools required to run FoundHaplo

#download GenotypeHarmonizer
mkdir $MAIN_PATH/software
mkdir $MAIN_PATH/software/GenotypeHarmonizer
cd $MAIN_PATH/software/GenotypeHarmonizer

wget https://github.com/molgenis/systemsgenetics/releases/download/1.4.0_20-8.1/GenotypeHarmonizer-1.4.23-dist.tar.gz
tar -xvf GenotypeHarmonizer-1.4.23-dist.tar.gz

#download Plink based on the OS using https://zzz.bwh.harvard.edu/plink/plink2.shtml

#download ANNOVAR using https://annovar.openbioinformatics.org/en/latest/user-guide/download/ by getting registerered


# VCFtools  https://vcftools.github.io/examples.html

# BCFtools  http://www.htslib.org/download/

# htslib http://www.htslib.org/download/

# nextflow https://www.nextflow.io/
