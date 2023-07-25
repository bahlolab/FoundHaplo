#' Takes a data frame with phased genotype data with columns in VCF format and separates haplotypes into two columns
#'
#' @description
#' sample.A with alleles 0|1 will be split to two columns with sample.A:a=0 and sample.A:b=1
#' @param haplotype_file A dataframe file containing phased SNP genotyped data in VCF format (including columns in a typical VCF file) after removing meta data
#' @return A data frame with haplotypes of test individuals in two columns
#' @import vcfR
#' @importFrom dplyr mutate_all
#' @export
#' @examples
#' haplotype_info=Create_hap_VCF(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)))

Create_hap_VCF=function(haplotype_file)
{
  if(!is.data.frame(haplotype_file)){stop("haplotype_file is not a data frame")}

  genotype_info=haplotype_file[,-c(1:9)]
  sample_column_names=colnames(haplotype_file)[-c(1:9)]
  haplotype_info=data.frame(matrix(NA,nrow=nrow(haplotype_file),ncol=length(sample_column_names)*2)) #initiate data frame to save haplotype data

  n1=sample_column_names # initiate vector to save the first haplotype of the samples
  n2=sample_column_names # initiate vector to save the second haplotype of the samples
  n1=paste(n1,sep=":","a")
  n2=paste(n2,sep=":","b")
  n3=data.frame(n1,n2)

  n3[] <- lapply(n3, as.character)
  haplotype_names_in_order=vector("character",(2*nrow(n3)))
  k=1
  for (i in 1:nrow(n3))
  {
    for(j in 1:2)
    {
      haplotype_names_in_order[k]=n3[i,j]
      k=k+1
    }
  }

  # haplotype_names_in_order has both haplotypes of the same sample together
  colnames(haplotype_info)=haplotype_names_in_order

  # Seperating the alleles to two haplotypes
  for(i in 1:ncol(genotype_info)){
    haplotype_info[,(2*i-1)]=substr(genotype_info[,i],0,1)
    haplotype_info[,(2*i)]=substr(genotype_info[,i],3,3)
  }
  # Seperating the alleles to two haplotypes

  #merge first 9 columns of the original haplotype_file and the haplotype_info into one dataframe
  haplotype_info=as.data.frame(cbind(haplotype_file[,c(1:9)],haplotype_info))

  haplotype_info=haplotype_info %>% mutate_all(as.character)

  #Further makes sure that we have only biallelic markers in the file
  haplotype_info=subset(haplotype_info,nchar(haplotype_info$REF)<2 & nchar(haplotype_info$ALT)<2)

  return(haplotype_info)

}



