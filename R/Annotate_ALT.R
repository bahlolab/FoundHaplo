#' Annotate the ALT column of bi-allelic SNPs in a GATK VCF file
#'
#' @description
#' Annotate the ALT column of bi-allelic SNPs in a GATK VCF file
#' @details
#' Make sure the input_vcf chromosome column doesn't have the "chr" prefix.
#' Takes a GATK generated VCF file with homozygous calls and annotate the ALT column using ANNOVAR reference files.
#' Download ANNOVAR reference .txt files and specify the path to the main folder in "dir_input_ref"
#' @param input_vcf Path to a VCF file  (type \code{"character"})
#' @param dir_input_ref Directory path to ANNOVAR reference .txt files with general MAF  (type \code{"character"})
#' @param dir_output Path to a directory to save the output VCF files (type \code{"character"})
#' @param output_file_name Name of the annotated vcf file (type \code{"character"})
#' @return Write a valid new VCF file with ALT column annotated
#' @import vcfR
#' @importFrom stringr str_remove
#' @importFrom plyr join
#' @importFrom dplyr mutate_all
#' @export
#' @examples
#' orig_dir <- getwd()
#' temp_dir <- tempdir()
#' setwd( temp_dir )
#' write.vcf(FAME1_test_cohort_GATK,paste0(temp_dir,"/","FAME1_test_cohort_GATK.vcf.gz"))
#' write.table(trimmed_chr8_hg19_ALL.sites.2015_08,paste0(temp_dir,"/","hg19_ALL.sites.2015_08.txt"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' Annotate_ALT(input_vcf = paste0(temp_dir,"/","FAME1_test_cohort_GATK.vcf.gz"), # read input from tempdir()
#'                          dir_input_ref = temp_dir,
#'                          dir_output = temp_dir,
#'                          output_file_name = "example_Annotate_ALT")
#' setwd( orig_dir )


Annotate_ALT=function(input_vcf,dir_input_ref,dir_output,output_file_name)
{
input_file_vcf=read.vcfR(file=input_vcf,verbose = FALSE)

data_file=as.data.frame(cbind(input_file_vcf@fix,input_file_vcf@gt))
chr=unique(data_file$CHROM)
chr=str_remove(chr, "chr")

ALL=read.delim(paste0(dir_input_ref,"/hg19_ALL.sites.2015_08.txt"),header=FALSE,stringsAsFactors=FALSE)


# REF and ALT alleles
ALL.f=subset(ALL,ALL$V1 %in% chr)

ALL.f=ALL.f[!(duplicated(ALL.f$V2) | duplicated(ALL.f$V2, fromLast = TRUE)), ]

ALL.f=subset(ALL.f,nchar(ALL.f$V3)<2 & nchar(ALL.f$V4)<2)

colnames(ALL.f)[2:4]=c("POS","REF","ALT")
ALL.alleles=ALL.f[,c(2,3,4)]

ALL.alleles.REF=ALL.alleles[,c(1,2)]
colnames(ALL.alleles.REF)[2]="REF.1"

ALL.alleles.ALT=ALL.alleles[,c(1,3)]
colnames(ALL.alleles.ALT)[2]="ALT.1"


data_file=data_file[!duplicated(data_file[c(2)]),]

data_file <- join(data_file, ALL.alleles.REF, by = "POS")
data_file <- join(data_file, ALL.alleles.ALT, by = "POS")

data_file=subset(data_file,!is.na(data_file$ALT.1))
data_file=subset(data_file,!is.na(data_file$REF.1))
data_file$POS=as.numeric(paste(data_file$POS))


data_file=data_file %>% mutate_all(as.character)


data_file=subset(data_file,data_file$REF==data_file$REF.1)


data_file=subset(data_file,!data_file$ALT.1 %in% 1)

attach(data_file)
data_file=data_file[order(POS),]
data_file$ALT=data_file$ALT.1
data_file=data_file[,!(names(data_file) %in% c("REF.1", "ALT.1"))]

colnames(data_file)[1]="#CHROM"
data_file[["QUAL"]][is.na(data_file[["QUAL"]])] <- "."
data_file[["ID"]][is.na(data_file[["ID"]])] <- "."
data_file[["FILTER"]][is.na(data_file[["FILTER"]])] <- "."
data_file[["INFO"]][is.na(data_file[["INFO"]])] <- "."
#####
write.table(input_file_vcf@meta,paste0(dir_output,"/",output_file_name,".vcf"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
write.table(data_file,paste0(dir_output,"/",output_file_name,".vcf"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = TRUE,append=TRUE)


}
