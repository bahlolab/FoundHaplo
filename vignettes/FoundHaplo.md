---
title: "FoundHaplo Introduction"
author: "Erandee Robertson"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
%\VignetteIndexEntry{FoundHaplo Introduction}
%\VignetteEngine{knitr::rmarkdown}
%\VignetteEncoding{UTF-8}
---

Individuals who inherit the same genetic mutation from a common ancestor also share genomic regions either side of the shared disease-causing variant. This suggests that the presence of a disease-causing genetic variant can be inferred by assessing identity by descent sharing of the variant-associated haplotype between an individual known to have the disease-causing mutation and a patient with unknown aetiology.

Using a Hidden Markov Model, Foundhaplo calculates the FH score, which is the likelihood of identity by descent between the samples of interest and a known disease causing haplotype and evaluate the resulted scores against a control cohort to predict individuals that carry disease causing haplotypes using SNP chip data.


## Data Format

There are three types of input data required by FoundHaplo R package

1. bgzipped imputed VCF file with individuals with known diseases and atleast one of their family memebers to create the disease haplotype
2. bgzipped imputed VCF file with test samples
3. bgzipped VCF file with 1000 Genomes control samples for a specific population from EUR,AMR,SAS,EAS and AFR.

Disease and test VCF files are generated from PLINK BED, BIM and FAM files, which are the binary files of .ped/map (http://www.cog-genomics.org/plink2).

### VCF files for disease and test cohort 

The input for FoundHaplo is imputed  bgzipped VCF files with haplotype data with below columns.

1. CHROM  : Chromosome 
2. POS : Base-pair position
3. ID : rs# or variant identifier
4. REF :  Reference allele
5. ALT  : Alternate allele
6. QUAL : A quality score associated with the inference of the given alleles.
7. FILTER : A flag indicating which of a given set of filters the variation has failed or PASS if all the filters were passed successfully.
8. INFO : An extensible list of key-value pairs (fields) describing the variation. See below for some common fields. Multiple fields are separated by semicolons with optional values in the format
9. FORMAT :  An (optional) extensible list of fields for describing the samples
+ SAMPLES

INFO column in the disease cohort must have the gnomAD population frequencies annotated from ANNOVAR as explained in https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md

explain freq tags

### Phased VCF files for the 1000 Genomes control cohort

Refer https://github.com/bahlolab/FoundHaplo/blob/main/Documentation/Guide%20to%20run%20FoundHaplo.md to download and prepare 1000 Genomes haplotype VCFs.

### Reference Data

HapMap genetic maps are needed to annotate genetic map distances and are in https://github.com/bahlolab/Foundhaplo/tree/main/input_files/public_data/genetic_map_HapMapII_GRCh37
gnomAD population frequencies to disease cort VCF files are annotated from ANNOVR, https://annovar.openbioinformatics.org/en/latest/user-guide/download/

## Key Functions

Function | Description
------------------ | --------------------------------------------------
Find_bp_to_trim | Outputs the base pair positions corresponding to +-20cM from a given locus
Create_hap_VCF | Takes a phased input data file with columns in VCF format and separate haplotypes into two columns
Convert_bp_cM | Adds a column with centiMorgan values to the input file
Phasing_by_pedigree | Phase samples with a disease variant by pedigrees
Calculate_IBD |  Calculate the FH score and IBD details between a known disease haplotype and a test individual for a disease variant of interest
Generate_FH_score | A wrapper function to generate FH scores and IBD details for each test/control - disease pair

## Example

Below is an IBD analysis performed on simulated data to demonstrate the steps involved in FoundHaplo. Change the parameters for each function to see how it affects the results.

**NOTE: Generate_FH_score must not be run on R studio and only be run on command line as it uses system() function to invoke vcftools and bcftools queries whcih can only be accessed from commanline.**

### Calculate the corresponding base pair positions

First we calculate the corresponding base pair positions of 10cM both sides of the disease causing variant.

```{r}
# load FoundHaplo library
library(FoundHaplo)
#lets look at an example Hapmap genetic map file
str(genetic_map_GRCh37_chr8)
orig_dir <- getwd()
temp_dir <- tempdir()
setwd( temp_dir )
#save the Hapmap genetic map file
write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE) # save genetic_map_GRCh37_chr8 as a text file to read from
Find_bp_to_trim(input_vector=c("FAME1.chr8.119379052."),dir_geneticMap=temp_dir,output_file=paste0(temp_dir,"/DCV_bp.txt"))
setwd( orig_dir )
```

### Create haplotype columns from the test cohort VCF file

```{r}
# FAME1_test_cohort is a VCF file
haplotype_info=Create_hap_VCF(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)))
```
### Adds a column with centiMorgan values to the input file

```{r}
orig_dir <- getwd()
temp_dir <- tempdir()
setwd( temp_dir )
write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE) # save genetic_map_GRCh37_chr8 as a text file to read from
haplotype_file_cM=Convert_bp_cM(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)),DCV="FAME1.chr8.119379052.",dir_geneticMap=temp_dir)
setwd( orig_dir )
```
### Phase samples with a disease variant by pedigrees

```{r}
orig_dir <- getwd()
temp_dir <- tempdir()
setwd( temp_dir )
library(vcfR)
write.vcf(FAME1_disease_cohort,paste0(temp_dir,"/","FAME1_disease_cohort.vcf.gz")) # save FAME1_disease_cohort as a VCF file to read from
sample_info # with three lines to generate three disease haplotypes in three seperate VCF files
sample_info=data.frame(rbind(c("HG00362_1,HG00362_2","duo"),c("NA11920,Affected_parent_NA11920,Unaffected_parent_NA11920","trio"),c("HG00313_1,HG00313_2","duo")))
write.table(sample_info,paste0(temp_dir,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save sample_info as a tab delimitted text file to read from
Phasing_by_pedigree(input_vcf = paste0(temp_dir,"/FAME1_disease_cohort",".vcf.gz"),
dir_output = temp_dir,
sample_info_file = paste0(temp_dir,"/","sample_info.txt"))
# temp_dir will have disease haplotypes in VCF format. disease haplotype is named as "h1" with homozygous genotypes in every VCF file.
setwd( orig_dir )
```
### Calculate the FH score and IBD details between a known disease haplotype and a test individual for a disease variant of interest

```{r}
library(vcfR)
library(data.table)
library(dplyr)
orig_dir <- getwd()
temp_dir <- tempdir()
setwd( temp_dir )
write.vcf(FAME1_disease_cohort,paste0(temp_dir,"/","FAME1_disease_cohort.vcf.gz")) # save FAME1_disease_cohort as a VCF file to read from
sample_info=data.frame(rbind(c("HG00362_1,HG00362_2","duo"),c("NA11920,Affected_parent_NA11920,Unaffected_parent_NA11920","trio"),c("HG00313_1,HG00313_2","duo")))
write.table(sample_info,paste0(temp_dir,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save sample_info as a tab delimitted text file to read from
Phasing_by_pedigree(input_vcf = paste0(temp_dir,"/FAME1_disease_cohort",".vcf.gz"),
dir_output = temp_dir,
sample_info_file = paste0(temp_dir,"/","sample_info.txt"))
disease_file <-fread(paste0(temp_dir,"/",paste(c("NA11920","Affected_parent_NA11920","Unaffected_parent_NA11920"), collapse = ','),".vcf"), skip = "#CHROM")
write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE) # save genetic_map_GRCh37_chr8 as a text file to read from
test_haplotype_info=Create_hap_VCF(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)))
test_haplotype_info_cM=Convert_bp_cM(haplotype_file=test_haplotype_info,DCV="FAME1.chr8.119379052.",dir_geneticMap=temp_dir)
MAF=sapply(strsplit(disease_file$INFO,";",fixed=TRUE),"[[", 9)
MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[", 2)
disease_file=as.data.frame(cbind(disease_file[,c("#CHROM","POS","REF","ALT")],MAF,disease_file[,"h1"]))
disease_file=disease_file[,c("#CHROM","POS","REF","ALT","MAF","h1")]
h1=substr(disease_file[,"h1"],0,1) # h1 is the disease haplotype. get the first allele of every marker in the VCF file. both alleles are the same in database_file
disease_file[,"h1"]=h1
colnames(disease_file)[ncol(disease_file)]=paste(c("NA11920","Affected_parent_NA11920","Unaffected_parent_NA11920"), collapse = ',') # add the name of the disease individual into the database_file
disease_file$MAF=as.numeric(disease_file$MAF)
disease_file=subset(disease_file,disease_file$MAF>0)
colnames(test_haplotype_info_cM)[1]="#CHROM"
common_markers=merge(disease_file, test_haplotype_info_cM, by=c("#CHROM","POS", "REF","ALT"))
fix=common_markers[,c("#CHROM","POS","REF","ALT")]
test_haplotype_info_cM=merge(test_haplotype_info_cM, fix, by=c("#CHROM","POS", "REF","ALT"))
disease_file=merge(disease_file, fix, by=c("#CHROM","POS", "REF","ALT"))
disease_file=Convert_bp_cM(haplotype_file=disease_file,DCV="FAME1.chr8.119379052.",dir_geneticMap=temp_dir)
disease_file=disease_file[,which(!colnames(disease_file) %in% c("#CHROM","REF","ALT"))]
final_file=as.data.frame(cbind(disease_file,test_haplotype_info_cM[,c(11,11+1)]))
final_file=subset(final_file,final_file$MAF>0)
is.na(final_file) <- final_file=="."
final_file=na.omit(final_file)
final_file=distinct(final_file,position_cM,.keep_all= TRUE)
attach(final_file)
final_file <- final_file[order(POS),]
detach(final_file) # final_file is ready to compare disease and the two test haplotypes for IBD relateness
Final_IBD_score <-Calculate_IBD(final_file,"FAME1.chr8.119379052.",dir_geneticMap=temp_dir)
setwd( orig_dir )
```

### A wrapper function to generate FH score and IBD details for each test/control - disease pair

```{r}
# This function will call Create_hap_VCF , haplotype_file_cM and Calculate_IBD
orig_dir <- getwd()
setwd(tempdir())
file.remove(list.files())
if(!dir.exists(paste0(tempdir(), "/1"))){dir.create(paste0(tempdir(), "/1"))} #to save disease haplotypes
if(!dir.exists(paste0(tempdir(), "/2"))){dir.create(paste0(tempdir(), "/2"))} #to save test_list of sample names in .txt files
if(!dir.exists(paste0(tempdir(), "/3"))){dir.create(paste0(tempdir(), "/3"))}  # dir_controls_file
if(!dir.exists(paste0(tempdir(), "/4"))){dir.create(paste0(tempdir(), "/4"))} # to save IBD results
temp_dir <- tempdir() # To carryout the main workload
library(vcfR)
write.vcf(FAME1_disease_cohort,paste0(temp_dir,"/","FAME1_disease_cohort.vcf.gz")) # save FAME1_disease_cohort as a VCF file to read from
sample_info=data.frame(rbind(c("HG00362_1,HG00362_2","duo"),c("NA11920,Affected_parent_NA11920,Unaffected_parent_NA11920","trio"),c("HG00313_1,HG00313_2","duo")))
write.table(sample_info,paste0(temp_dir,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save sample_info as a tab delimitted text file to read from
Phasing_by_pedigree(input_vcf = paste0(temp_dir,"/FAME1_disease_cohort",".vcf.gz"),
dir_output = paste0(temp_dir,"/1"),
sample_info_file = paste0(temp_dir,"/","sample_info.txt"))
write.vcf(FAME1_test_cohort,paste0(temp_dir,"/","FAME1_test_cohort.vcf.gz")) # save FAME1_test_cohort as a VCF file to read from
write.table(file00,paste0(temp_dir,"/2/","file00",".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save file00 with maximum of 100 sample names of the FAME1_test_cohort as a text file to read from
write.vcf(FAME1_control_cohort,paste0(temp_dir,"/3/","FAME1.chr8.vcf.gz")) # save FAME1_control_cohort as a VCF file to read from
Generate_FH_score(DCV="FAME1.chr8.119379052.",minor_allele_cutoff=0,imputation_quality_score_cutoff_test=0,frequency_type="EUR",dir_geneticMap=temp_dir,dir_disease_files=paste0(temp_dir,"/1"),test_file=paste0(temp_dir,"/","FAME1_test_cohort.vcf.gz"),test_name="FAME1_example_test_cohort",test_list=paste0(temp_dir,"/2/","file00.txt"),data_type="test",dir_controls_file=paste0(temp_dir,"/3"),dir_to_save_report=paste0(temp_dir,"/4"))
setwd(paste0(temp_dir,"/4"))
read.delim(list.files(paste0(temp_dir,"/4"))[1],header=FALSE)
setwd(orig_dir)
```

### Plotting Results


