## All the functions are explained as a workflow below
##-----------------------------------------------------------------------------
# load FoundHaplo library
library(FoundHaplo)
library(vcfR)
library(data.table)
library(dplyr)

file.remove(list.files())
# lets look at an example Hapmap genetic map file
str(genetic_map_GRCh37_chr8)
orig_DIR <- getwd()
temp_DIR <- tempdir()
setwd( temp_DIR )
# save the Hapmap genetic map file to read from
write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE) # save genetic_map_GRCh37_chr8 as a text file to read from
# Find_bp_to_trim gives the corresponding base pair positions of +-10cM from a given disease locus
Find_bp_to_trim(input_vector=c("FAME1.chr8.119379052."),geneticMap_DIR=temp_DIR,output_file=paste0(temp_DIR,"/DCV_bp.txt"))

## -----------------------------------------------------------------------------
# FAME1_test_cohort is a VCF file with samples to test using FoundHaplo
#Create_hap_VCF seperates haplotype data into two columns of a dataframe.
test_haplotype_info=Create_hap_VCF(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)))
head(test_haplotype_info)
## -----------------------------------------------------------------------------
# Convert_bp_cM adds a new column with centiMorgan positions corresponding to the base pairs in haplotype_file data
test_haplotype_info_cM=Convert_bp_cM(haplotype_file=test_haplotype_info,DCV="FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
head(test_haplotype_info_cM)
## -----------------------------------------------------------------------------
if(!dir.exists(paste0(tempdir(), "/1"))){dir.create(paste0(tempdir(), "/1"))} # create a temporafy directory to save disease haplotypes

# FAME1_disease_cohort is a VCF file with genotype data of disease individuals
write.vcf(FAME1_disease_cohort,paste0(temp_DIR,"/","FAME1_disease_cohort.vcf.gz")) # save FAME1_disease_cohort as a VCF file to read from
sample_info=data.frame(rbind(c("HG00362_1_HG00362_1,HG00362_2_HG00362_2","duo"),c("NA11920_NA11920,Affected_parent_NA11920_Affected_parent_NA11920,Unaffected_parent_NA11920_Unaffected_parent_NA11920","trio"),c("HG00313_1_HG00313_1,HG00313_2_HG00313_2","duo")))
sample_info # sample_info has three lines with instructions to pedigree phase three disease haplotypes in three seperate VCF files, using the FAME1_disease_cohort
write.table(sample_info,paste0(temp_DIR,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE) # save sample_info as a tab delimitted text file to read from
# Phasing_by_pedigree will create disease haplotypes in a given directory
Phasing_by_pedigree(input_vcf = paste0(temp_DIR,"/FAME1_disease_cohort",".vcf.gz"),
                    output_DIR = paste0(tempdir(), "/1"),
                    sample_info_file = paste0(temp_DIR,"/","sample_info.txt"))
# temp_DIR/1 will have disease haplotypes in VCF format. disease haplotype column in each file is named as "h1" with homozygous genotypes in every VCF file.

## -----------------------------------------------------------------------------
# Read one of the created disease haplotype files to test using FoundHaplo
disease_file <-fread(paste0(tempdir(), "/1","/",paste(c("NA11920_NA11920","Affected_parent_NA11920_Affected_parent_NA11920","Unaffected_parent_NA11920_Unaffected_parent_NA11920"), collapse = ','),".vcf"), skip = "#CHROM")
MAF=sapply(strsplit(disease_file$INFO,";",fixed=TRUE),"[[", 6) # Taking "AF_nfe" as MAF
MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[", 2)
disease_file=as.data.frame(cbind(disease_file[,c("#CHROM","POS","REF","ALT")],MAF,disease_file[,"h1"]))
disease_file=disease_file[,c("#CHROM","POS","REF","ALT","MAF","h1")]
h1=substr(disease_file[,"h1"],0,1) # h1 is the disease haplotype. get the first allele of every marker in the VCF file. both alleles are the same in database_file
disease_file[,"h1"]=h1
colnames(disease_file)[ncol(disease_file)]=paste(c("NA11920","Affected_parent_NA11920","Unaffected_parent_NA11920"), collapse = ',') # add the name of the disease individual into the database_file
disease_file$MAF=as.numeric(disease_file$MAF)
disease_file=subset(disease_file,disease_file$MAF>0)
# test_haplotype_info_cM is the file with test samples
colnames(test_haplotype_info_cM)[1]="#CHROM"
common_markers=merge(disease_file, test_haplotype_info_cM, by=c("#CHROM","POS", "REF","ALT"))
fix=common_markers[,c("#CHROM","POS","REF","ALT")]
test_haplotype_info_cM=merge(test_haplotype_info_cM, fix, by=c("#CHROM","POS", "REF","ALT"))
disease_file=merge(disease_file, fix, by=c("#CHROM","POS", "REF","ALT"))
disease_file=Convert_bp_cM(haplotype_file=disease_file,DCV="FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
disease_file=disease_file[,which(!colnames(disease_file) %in% c("#CHROM","REF","ALT"))]
final_file=as.data.frame(cbind(disease_file,test_haplotype_info_cM[,c(11,11+1)]))
final_file=subset(final_file,final_file$MAF>0)
is.na(final_file) <- final_file=="."
final_file=na.omit(final_file)
final_file=distinct(final_file,position_cM,.keep_all= TRUE)
attach(final_file)
final_file <- final_file[order(POS),]
detach(final_file) # final_file is ready to compare disease and the two test haplotypes for IBD relateness
#Calculate_IBD generates the IBD report and saves in a given directory as a text file
Final_IBD_score <-Calculate_IBD(final_file,"FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
Final_IBD_score
setwd( orig_dir )

## -----------------------------------------------------------------------------
## Run the wrapper function Generate_FH_score that uses all the above functions
file.remove(list.files())
orig_DIR <- getwd()
setwd(tempdir())
file.remove(list.files())
if(!dir.exists(paste0(tempdir(), "/1"))){dir.create(paste0(tempdir(), "/1"))} #to save disease haplotypes
if(!dir.exists(paste0(tempdir(), "/2"))){dir.create(paste0(tempdir(), "/2"))} #to save test_list of sample names in .txt files
if(!dir.exists(paste0(tempdir(), "/3"))){dir.create(paste0(tempdir(), "/3"))}  # controls_file_DIR
if(!dir.exists(paste0(tempdir(), "/4"))){dir.create(paste0(tempdir(), "/4"))} # to save IBD results
temp_DIR <- tempdir() # To carryout the main workload
library(vcfR)
# FAME1_disease_cohort is the disease file
write.vcf(FAME1_disease_cohort,paste0(temp_DIR,"/","FAME1_disease_cohort.vcf.gz"))
sample_info=data.frame(rbind(c("HG00362_1_HG00362_1,HG00362_2_HG00362_2","duo"),c("NA11920_NA11920,Affected_parent_NA11920_Affected_parent_NA11920,Unaffected_parent_NA11920_Unaffected_parent_NA11920","trio"),c("HG00313_1_HG00313_1,HG00313_2_HG00313_2","duo")))
write.table(sample_info,paste0(temp_DIR,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
# Phase FAME1_disease_cohort by pedigree
Phasing_by_pedigree(input_vcf = paste0(temp_DIR,"/FAME1_disease_cohort",".vcf.gz"),
                   output_DIR = paste0(temp_DIR,"/1"),
                   sample_info_file = paste0(temp_DIR,"/","sample_info.txt"))
# FAME1_test_cohort is the test cohort
write.vcf(FAME1_test_cohort,paste0(temp_DIR,"/","FAME1_test_cohort.vcf.gz"))
write.table(file00,paste0(temp_DIR,"/2/","file00",".txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
write.vcf(FAME1_control_cohort,paste0(temp_DIR,"/3/","FAME1.chr8.vcf.gz"))
write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
# Generate_FH_score calculates FH score values for every disease-test pair and save results in .txt files for every disease haplotype
Generate_FH_score(DCV="FAME1.chr8.119379052.",minor_allele_cutoff=0,gen_allele_mismatch_rate=0.01,MA_cutoff=-0.4,meiosis=1,imputation_quality_score_cutoff_test=0,frequency_type="EUR",geneticMap_DIR=temp_DIR,disease_files_DIR=paste0(temp_DIR,"/1"),test_file=paste0(temp_DIR,"/","FAME1_test_cohort.vcf.gz"),test_name="FAME1_example_test_cohort",test_list=paste0(temp_DIR,"/2/","file00.txt"),data_type="test",controls_file_DIR=paste0(temp_DIR,"/3"),save_report_DIR=paste0(temp_DIR,"/4"),temp_DIR)
setwd(paste0(temp_DIR,"/4"))
FH_IBD_scores=read.delim(list.files(paste0(temp_DIR,"/4"))[1],header=FALSE) # read the first .txt file with results to review
# Concatanate all the results in list.files(paste0(temp_DIR,"/4"))  into one text file for analysis
setwd(orig_DIR)

## -----------------------------------------------------------------------------
# Analyse_FH lists the samples predicted by the FoundHaplo algorithm and also produce graphical results in a given directory
orig_DIR <- getwd()
temp_DIR <- tempdir()
setwd( temp_DIR )
# FH_IBD_scores represents all the concatanated FH scores from multiple .txt files
write.table(FH_IBD_scores,"FH_IBD_scores.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
Analyse_FH("FH_IBD_scores.txt",temp_DIR,0.99)
setwd(orig_DIR)
