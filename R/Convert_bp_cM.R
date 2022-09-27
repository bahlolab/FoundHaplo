#' Adds a column with centiMorgan values to the input file
#'
#' @description
#' "position_cM" column will be added after column "POS" in to file with VCF format
#' Genome built should be GRCh37 and genetic map in "dir_geneticMap" must have chromosomes named with "chr" prefix, add the prefix if needed
#' @param haplotype_file File in VCF format after removing meta data
#' @param DCV Name of the disease causing variant of interest i.e FAME1.chr8.119379052
#' @param dir_geneticMap Directory to genetic_map_HapMapII_GRCh37 location
#' @return A file with cM distance annotated
#' @export
#' @examples
#' orig_dir <- getwd()
#' temp_dir <- tempdir()
#' setwd( temp_dir )
#' write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' haplotype_file_cM=Convert_bp_cM(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)),DCV="FAME1.chr8.119379052.",dir_geneticMap=temp_dir)
#' setwd( orig_dir )

Convert_bp_cM=function(haplotype_file,DCV,dir_geneticMap)
{
    #Extract the relevant chromosome of the DCV
    chr=strsplit(as.character(DCV),".",fixed=TRUE)
    chr=sapply(chr, "[[", 2)


    path=paste0(dir_geneticMap,"/genetic_map_GRCh37_",chr,".txt")
    recombination_map=read.delim(path)
    colnames(recombination_map)=c("Chromosome","position_bp","Rate.cM.Mb.","position_cM")


    position_bp=haplotype_file[,"POS"]
    position_cM=as.data.frame(matrix(0,ncol=1,nrow=nrow(haplotype_file))) # to store cM values
    position_bp_cM=cbind(position_bp,position_cM)
    colnames(position_bp_cM)=c("position_bp","position_cM")



    position_bp_cM$position_cM <- recombination_map$position_cM[match(position_bp_cM$position_bp, recombination_map$position_bp)]
    position_bp_cM$position_bp=as.numeric(paste(position_bp_cM$position_bp))

    #function to calculate cM distance using linear interpolation
    fun <- approxfun(recombination_map$position_bp,recombination_map$position_cM,ties=mean)
    position_cM=fun(position_bp_cM$position_bp)

    haplotype_file=add_column(haplotype_file,position_cM, .after = 2) # Make cM column the third column
    haplotype_file=subset(haplotype_file,!is.na(haplotype_file$position_cM))

    return(haplotype_file)


}


