#' Adds a column with centiMorgan values to the input file
#'
#' @description
#' "position_cM" column will be added after column "POS" in to a data frame file, where first nine columns are in VCF format.
#' Genome built should be GRCh37 and genetic map in "geneticMap_DIR" must have chromosomes named with "chr" prefix, add the prefix if needed.
#' @param haplotype_file A dataframe file including genotype data, with first nine columns in VCF format after removing meta data, if any
#' @param DCV Name of the disease causing variant of interest i.e FAME1.chr8.119379052
#' @param geneticMap_DIR Directory to genetic_map_HapMapII_GRCh37 location
#' @return A dataframe with cM distance annotated
#' @import tibble
#' @export
#' @examples
#' orig_DIR <- getwd()
#' temp_DIR <- tempdir()
#' setwd( temp_DIR )
#' write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' haplotype_file_cM=Convert_bp_cM(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)),DCV="FAME1.chr8.119379052",geneticMap_DIR=temp_DIR)
#' setwd(orig_DIR)

Convert_bp_cM=function(haplotype_file,DCV,geneticMap_DIR)
{
  if(!dir.exists(geneticMap_DIR)){stop("geneticMap_DIR does not exist")}
  if(!is.data.frame(haplotype_file)){stop("haplotype_file is not a data frame")}

  #Extract the relevant chromosome of the DCV
  chr=strsplit(as.character(DCV),".",fixed=TRUE)
  chr=sapply(chr, "[[", 2)

  recombination_map=read.delim(paste0(geneticMap_DIR,"/genetic_map_GRCh37_",chr,".txt"))
  colnames(recombination_map)=c("Chromosome","position_bp","Rate.cM.Mb.","position_cM")

  position_bp=haplotype_file[,"POS"]
  position_cM=as.data.frame(matrix(0,ncol=1,nrow=nrow(haplotype_file))) # to store cM values
  position_bp_cM=cbind(position_bp,position_cM)
  colnames(position_bp_cM)=c("position_bp","position_cM")

  position_bp_cM$position_cM <- recombination_map$position_cM[match(position_bp_cM$position_bp, recombination_map$position_bp)]
  position_bp_cM$position_bp=as.numeric(paste(position_bp_cM$position_bp))

  #function to calculate cM distance using linear interpolation
  fun_to_cM <- approxfun(recombination_map$position_bp,recombination_map$position_cM,ties=mean)
  position_cM=fun_to_cM(position_bp_cM$position_bp)

  haplotype_file=add_column(haplotype_file,position_cM, .after = 2) # Make cM column the third column
  haplotype_file=subset(haplotype_file,!is.na(haplotype_file$position_cM))

  return(haplotype_file)


}


