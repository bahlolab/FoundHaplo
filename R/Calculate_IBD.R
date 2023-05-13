#' Calculate the FH score and IBD details between a known disease haplotype and a test individual for a disease variant of interest
#'
#' @description
#' Calculate the FH score and IBD details between a known disease haplotype and a test individual for a disease variant of interest.
#' Genome built should be GRCh37 and genetic map in "geneticMap_DIR" must have chromosomes named with "chr" prefix, add the prefix if needed
#' @param data_file A data frame with 6 columns containing genotype data of disease haplotype and test individual with
#' \enumerate{
#' \item position_bp, position in base pairs (type \code{"numeric"})
#' \item position_cM, position in centi Morgans (type \code{"numeric"})
#' \item Minor_allele_frequency,Minor allele frequency (type \code{"numeric"})
#' \item disease_haplotype, known disease haplotype (type \code{"numeric"})
#' \item test_id:a, haplotype 1 of the test individual (type \code{"numeric"})
#' \item test_id:b, haplotype 2 of the test individual  (type \code{"numeric"})
#' }
#' Make sure that data_file sample names do not have underscore
#' @param DCV Name of the disease causing variant of interest i.e FAME1.chr8.119379052
#' @param geneticMap_DIR Directory path to genetic_map_HapMapII_GRCh37 location 
#' @param gen_allele_mismatch_rate Genotype and imputation error rate allowed, default is 0.1
#' @param MA_cutof Moving average threshold for allowing genotype and imputation errors (derived based on simulation studies), default is -0.4
#' @param meiosis Estimated number of meiosis between disease-test pair, default is 1
#' @return A vector with all the details of the IBD sharing to analyze later :
#' \enumerate{
#' \item FH_score (type \code{"numeric"})
#' \item left_LLR (type \code{"numeric"})
#' \item right_LLr (type \code{"numeric"})
#' \item total_cM_sharing (type \code{"numeric"})
#' \item total_left_cM_sharing (type \code{"numeric"})
#' \item total_right_cM_sharing (type \code{"numeric"})
#' \item number_of_allele_mismatches_in_the_markov_chain (type \code{"numeric"})
#' \item number_of_markers_in_the_markov_chain (type \code{"numeric"})
#' \item numer_of_haplotype_switches_in_the_markov_chain (type \code{"numeric"})
#' \item snp_density_in_data_file (type \code{"numeric"})
#' \item total_number_of_markers_in_data_file (type \code{"numeric"})
#' \item total_cM_span_of_data_file (type \code{"numeric"})
#' }
#' @export
#' @examples
#' library(vcfR)
#' library(data.table)
#' library(dplyr)
#' orig_DIR <- getwd()
#' temp_DIR <- tempdir()
#' setwd( temp_DIR )
#' write.vcf(FAME1_disease_cohort,paste0(temp_DIR,"/","FAME1_disease_cohort.vcf.gz"))
#' sample_info=data.frame(rbind(c("HG00362_1,HG00362_2","duo"),c("NA11920,Affected_parent_NA11920,Unaffected_parent_NA11920","trio"),c("HG00313_1,HG00313_2","duo")))
#' write.table(sample_info,paste0(temp_DIR,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' Phasing_by_pedigree(input_vcf = paste0(temp_DIR,"/FAME1_disease_cohort",".vcf.gz"),
#'                    dir_output = temp_DIR,
#'                    sample_info_file = paste0(temp_DIR,"/","sample_info.txt"))
#' disease_file <-fread(paste0(temp_DIR,"/",paste(c("NA11920","Affected_parent_NA11920","Unaffected_parent_NA11920"), collapse = ','),".vcf"), skip = "#CHROM")
#' write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' test_haplotype_info=Create_hap_VCF(haplotype_file=data.frame(cbind(FAME1_test_cohort@fix,FAME1_test_cohort@gt)))
#' test_haplotype_info_cM=Convert_bp_cM(haplotype_file=test_haplotype_info,DCV="FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
#' MAF=sapply(strsplit(disease_file$INFO,";",fixed=TRUE),"[[", 9)
#' MAF=sapply(strsplit(MAF,"=",fixed=TRUE),"[[", 2)
#' disease_file=as.data.frame(cbind(disease_file[,c("#CHROM","POS","REF","ALT")],MAF,disease_file[,"h1"]))
#' disease_file=disease_file[,c("#CHROM","POS","REF","ALT","MAF","h1")]
#' h1=substr(disease_file[,"h1"],0,1) # h1 is the disease haplotype. get the first allele of every marker in the VCF file. both alleles are the same in database_file
#' disease_file[,"h1"]=h1
#' colnames(disease_file)[ncol(disease_file)]=paste(c("NA11920","Affected_parent_NA11920","Unaffected_parent_NA11920"), collapse = ',') # add the name of the disease individual into the database_file
#' disease_file$MAF=as.numeric(disease_file$MAF)
#' disease_file=subset(disease_file,disease_file$MAF>0)
#' colnames(test_haplotype_info_cM)[1]="#CHROM"
#' common_markers=merge(disease_file, test_haplotype_info_cM, by=c("#CHROM","POS", "REF","ALT"))
#' fix=common_markers[,c("#CHROM","POS","REF","ALT")]
#' test_haplotype_info_cM=merge(test_haplotype_info_cM, fix, by=c("#CHROM","POS", "REF","ALT"))
#' disease_file=merge(disease_file, fix, by=c("#CHROM","POS", "REF","ALT"))
#' disease_file=Convert_bp_cM(haplotype_file=disease_file,DCV="FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
#' disease_file=disease_file[,which(!colnames(disease_file) %in% c("#CHROM","REF","ALT"))]
#' final_file=as.data.frame(cbind(disease_file,test_haplotype_info_cM[,c(11,11+1)]))
#' final_file=subset(final_file,final_file$MAF>0)
#' is.na(final_file) <- final_file=="."
#' final_file=na.omit(final_file)
#' final_file=distinct(final_file,position_cM,.keep_all= TRUE)
#' attach(final_file)
#' final_file <- final_file[order(POS),]
#' detach(final_file)
#' Final_IBD_score <-Calculate_IBD(final_file,"FAME1.chr8.119379052.",geneticMap_DIR=temp_DIR)
#' setwd( orig_DIR )

Calculate_IBD=function(data_file,DCV,geneticMap_DIR,gen_allele_mismatch_rate=0.01,MA_cutof-0.4,meiosis=1)
{
  # data_file is the main file with input genotype data
  
  chr=strsplit(as.character(DCV),".",fixed=TRUE)
  chr=sapply(chr, "[[", 2)
  
  recombination_map=read.delim(paste0(geneticMap_DIR,"/genetic_map_GRCh37_",chr,".txt"))
  colnames(recombination_map)=c("Chromosome","position_bp","Rate.cM.Mb.","position_cM")
  
  g1=gen_allele_mismatch_rate # replacing long variable name (g1=gen_allele_mismatch_rate)
  m=meiosis # replacing long variable name (m=meiosis)
  
  colnames(data_file)=c("position_bp","position_cM","Minor_allele_frequency",colnames(data_file)[4:6])
  
  data_file[,1]=as.numeric(data_file[,1])
  data_file[,2]=as.numeric(data_file[,2])
  data_file[,3]=as.numeric(data_file[,3])
  data_file[,4]=as.numeric(data_file[,4])
  data_file[,5]=as.numeric(data_file[,5])
  data_file[,6]=as.numeric(data_file[,6])
  
  data_file <- data_file[!is.na(data_file[,4]), ]
  data_file <- data_file[!is.na(data_file[,5]), ]
  data_file <- data_file[!is.na(data_file[,6]), ]
  
  data_file <- data_file[data_file[,4]<2, ] #further making sure only biallelic SNPs are in the data_file
  data_file <- data_file[data_file[,5]<2,]
  data_file <- data_file[data_file[,6]<2,]
  
  DCV_adjusted=strsplit(DCV, ".",fixed=TRUE)
  DCV_adjusted=as.vector(unlist(DCV_adjusted,recursive = FALSE))
  logical_vector=(data_file$position_bp>as.numeric(DCV_adjusted[3]))
  xR=table(logical_vector)["FALSE"]+1 # xR is the immediate right marker to the DCV locus
  xL=xR-2 # xL is the immediate left marker to the DCV locus
  
  #xL=xR-2 is taken to avoid the DCV locus in case if its in the data_file (when DCV is a point mutation like GEFS+)
  
  DCV_adjusted=as.numeric(DCV_adjusted[3]) # take the base pair position of the DCV locus
  
  #function to calculate cM distance using linear interpolation
  fun_to_cM <- approxfun(recombination_map$position_bp,recombination_map$position_cM,ties=mean)
  DCV_cM=fun_to_cM(DCV_adjusted)
  
  # IBD0 is the log of the joint probability of haplotypes not being IBD along the markers traversed
  # IBD1 is the log of the joint probability of haplotypes being IBD along the markers traversed
  
  # Initialize IBD values at DCV locus
  IBD0=0
  IBD1=0
  
  
  k=1 # k is used to traverse along markers; from 1 to end of sharing
  a=5 # a shows the haplotype of the test individual the markov chain is on, we start from the first haplotype which is on column 5.
  A=c(5,6) # A shows the set of both haplotypes of the test individual i.e column 5 and 6
  j=xL # j is used to traverse along markers; from xL to end of sharing
  # j and k store different values, j starts with 1 and k starts with the first marker to the left.
  i=4 # column index of the disease haplotype in data_file
  
  cumulative_IBD=vector(mode = "numeric",length=xL) # vector to save IBD values at each marker
  allele_mismatch=vector(mode = "numeric",length=xL) # vector to save allele mismatches at each marker; 0= Atleast one allele of the test sample matches with the disease haplotype, 1= both alleles do not match with the disease haplotype
  test_haplotype=vector(mode = "character",length=xL) # vector to save the haplotype the model is on
  IBD_difference=vector(mode = "numeric",length=xL) # vector to save the IBD difference between adjacent markers
  distance_from_the_DCV_cM=vector(mode = "numeric",length=xL) # vector to save the distance markov chain has traversed in cM from the DCV
  
  #To the left of the DCV locus ;Traverse from xL to end of sharing
  
  for(j in xL:1)
  {
    
    f1=data_file[j,3]         # f1=Minor/Alternate allele frequency
    if(is.na(f1))
    {
      next
      
    }
    
    Morgan_distance_between_adjacent_markers=(data_file[j+1,"position_cM"]-data_file[j,"position_cM"])/100
    d=Morgan_distance_between_adjacent_markers # replacing long variable name (d=Morgan_distance_between_adjacent_markers)
    
    r1=1-exp(-m*d) # r1 = Probability of recombination at the marker
    r0=exp(-m*d) # r0 = Probability of no recombination at the marker
    g0=1-g1 # g0 = Probability of no genotype/imputation error at the marker
    f0=1-f1 # f0 = Probability of major allele at the marker
    
    IBD0=IBD0+log(P_no_IBD)
    
    IBD1=IBD1+log(P_IBD)
    
    allele_mismatch[k]=0
    test_haplotype[k]=a
    distance_from_the_DCV_cM[k]=DCV_cM-data_file[j,"position_cM"]
    
    if(data_file[j,a]==data_file[j,i] && data_file[j,i]==0)
    {
      P_no_IBD= (g1*f1 + g0*f0)^2 # Probability of no IBD between disease haplotype and the test haplotype
      
      P_IBD=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1 # Probability of IBD between disease haplotype and the test haplotype
    }
    
    else if(data_file[j,a]==data_file[j,i] && data_file[j,i]==1)
    {
      P_no_IBD= (g0*f1 + g1*f0)^2
      
      P_IBD=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
    }
    
    else if(data_file[j,A[!A %in% a]]==data_file[j,i] && data_file[j,i]==0)
    {
      a=A[!A %in% a]

      P_no_IBD= (g1*f1 + g0*f0)^2
      
      P_IBD=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
    }
    
    else if(data_file[j,A[!A %in% a]]==data_file[j,i] && data_file[j,i]==1)
    {
      a=A[!A %in% a]
      
      P_no_IBD= (g0*f1 + g1*f0)^2
      
      P_IBD=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
    }
    
    else if(data_file[j,a]==1 && data_file[j,i]==0)
    {
      P_no_IBD= (g1*f1 + g0*f0) * (g0*f1 + g1*f0)
      
      P_IBD=g0*g1*f1*r0 + g0*g1*f0*r0 + g1*f1*g0*f1*r1 + g1*f1*g1*f0*r1 + g0*f0*g0*f1*r1 + g0*f0*g1*f0*r1
    }
    
    else if(data_file[j,a]==0 && data_file[j,i]==1)
    {
      P_no_IBD= (g0*f1 + g1*f0) * (g1*f1 + g0*f0)
      
      P_IBD=g0*r0*f1*g1 + g0*r0*f0*g1 + g0*g1*f1*f1*r1 + g0*f1*g0*f0*r1 + g1*f0*g1*f1*r1 + g1*f0*g0*f0*r1
    }
    
    cumulative_IBD[k]=-2*(IBD0 - IBD1)

    if(k>100)
    {
    IBD_difference[k-1]=cumulative_IBD[k]-cumulative_IBD[k-1]
    
      if(frollmean(IBD_difference[(k-100):(k-1)],100)<=MA_cutoff)    # 100 point moving average of IBD_difference
      {break}      #break the loop and stop if froll<=MA_cutoff
 
    }
    k=k+1
  }
  
  
  
  Max_left_IBD=max(cumulative_IBD) # the marker at which the markov chain gave the highest IBD to the left
  
  Max_left_cM=DCV_cM-data_file[xL-(which(cumulative_IBD==max(cumulative_IBD))[length(which(cumulative_IBD==max(cumulative_IBD)))]-1),"position_cM"] # the cM location at which the chain gave the highest IBD to the left
  
  number_of_allele_mismatches_in_the_left_markov_chain=sum(allele_mismatch[1:which(cumulative_IBD==max(cumulative_IBD))]) # total number of allele mismatches to the left
  number_of_markers_in_the_left_markov_chain=which(cumulative_IBD==max(cumulative_IBD)) # total number of markers traversed until end of left sharing
  
  numer_of_haplotype_switches_in_the_left_markov_chain=length(which(c(FALSE, tail(test_haplotype[1:which(cumulative_IBD==max(cumulative_IBD))[1]],-1) != head(test_haplotype[1:which(cumulative_IBD==max(cumulative_IBD))[1]],-1))))# number of times the chain switches between the two test haplotypes
  
  
  # in the same way we run the markov chain to the right of the DCV locus
  
  k=1 # k is used to traverse along markers; from 1 to end of sharing
  a=5 # a shows the haplotype of the test individual the markov chain is on, we start from the first haplotype which is on column 5.
  A=c(5,6) # A shows the set of both haplotypes of the test individual i.e column 5 and 6
  j=xR # j is used to traverse along markers; from xR to end of sharing
  # j and k store different values, j starts with 1 and k starts with the first marker to the left.
  i=4 # column index of the disease haplotype in data_file
  
  IBD0=0
  IBD1=0
  
  cumulative_IBD=vector(mode = "numeric",length=nrow(data_file)-xL) # vector to save IBD values at each marker
  allele_mismatch=vector(mode = "numeric",length=nrow(data_file)-xL) # vector to save allele mismatches at each marker; 0= Atleast one allele of the test sample matches with the disease haplotype, 1= both alleles do not match with the disease haplotype
  test_haplotype=vector(mode = "character",length=nrow(data_file)-xL) # vector to save the haplotype the model is on
  IBD_difference=vector(mode = "numeric",length=nrow(data_file)-xL) # vector to save the IBD difference between adjacent markers
  distance_from_the_DCV_cM=vector(mode = "numeric",length=nrow(data_file)-xL) # vector to save the distance markov chain has traversed in cM from the DCV
  
  #To the right of the DCV locus ;Traverse from xR to end of sharing

  for(j in xR:(nrow(data_file)-1))
  {
    
    f1=data_file[j,3]  # f1=Minor/Alternate allele frequency
    if(is.na(f1))
    {
      next
      
    }
    
    Morgan_distance_between_adjacent_markers=(data_file[j,"position_cM"]-data_file[j-1,"position_cM"])/100
    d=Morgan_distance_between_adjacent_markers # replacing long variable name (d=Morgan_distance_between_adjacent_markers)

    
    r1=1-exp(-m*d) # r1 = Probability of recombination at the marker
    r0=exp(-m*d) # r0 = Probability of no recombination at the marker
    g0=1-g1 # g0 = Probability of no genotype/imputation error at the marker
    f0=1-f1 # f0 = Probability of major allele at the marker
    
    IBD0=IBD0+log(P_no_IBD)
    
    IBD1=IBD1+log(P_IBD)
    
    allele_mismatch[k]=0
    test_haplotype[k]=a
    distance_from_the_DCV_cM[k]=data_file[j,"position_cM"]-DCV_cM
    
    if(data_file[j,a]==data_file[j,i] && data_file[j,i]==0){
      P_no_IBD= (g1*f1 + g0*f0)^2
      
      P_IBD=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1 #Probability of IBD
      
    }else if(data_file[j,a]==data_file[j,i] && data_file[j,i]==1){
      
      P_no_IBD= (g0*f1 + g1*f0)^2
      
      P_IBD=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
      
    }else if(data_file[j,A[!A %in% a]]==data_file[j,i] && data_file[j,i]==0){
      
      a=A[!A %in% a]
      
      P_no_IBD= (g1*f1 + g0*f0)^2
      
      P_IBD=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
      
    }else if(data_file[j,A[!A %in% a]]==data_file[j,i] && data_file[j,i]==1){
      
      a=A[!A %in% a]

      P_no_IBD= (g0*f1 + g1*f0)^2
      
      P_IBD=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1

    }else if(data_file[j,a]==1 && data_file[j,i]==0){

      P_no_IBD= (g1*f1 + g0*f0) * (g0*f1 + g1*f0)
      
      P_IBD=g0*g1*f1*r0 + g0*g1*f0*r0 + g1*f1*g0*f1*r1 + g1*f1*g1*f0*r1 + g0*f0*g0*f1*r1 + g0*f0*g1*f0*r1
      
    }else if(data_file[j,a]==0 && data_file[j,i]==1){

      P_no_IBD= (g0*f1 + g1*f0) * (g1*f1 + g0*f0)
      
      P_IBD=g0*r0*f1*g1 + g0*r0*f0*g1 + g0*g1*f1*f1*r1 + g0*f1*g0*f0*r1 + g1*f0*g1*f1*r1 + g1*f0*g0*f0*r1
      
    }
    
    cumulative_IBD[k]=-2*(IBD0 - IBD1)
    
    if(k>100)
    {
      IBD_difference[k-1]=cumulative_IBD[k]-cumulative_IBD[k-1]
      
      if(frollmean(IBD_difference[(k-100):(k-1)],100)<=MA_cutoff)    # 100 point moving average of IBD_difference
      {break}      #break the loop and stop if froll<=MA_cutoff
      
    }
    
    k=k+1
    
  }
  
  Max_right_IBD=max(cumulative_IBD) # the marker at which the markov chain gave the highest IBD to the right
  
  Max_right_cM=data_file[xR+ which(cumulative_IBD==max(cumulative_IBD))[length(which(cumulative_IBD==max(cumulative_IBD)))] -1,"position_cM"]-DCV_cM # the cM location at which the chain gave the highest IBD to the right
  
  
  number_of_allele_mismatches_in_the_right_markov_chain=sum(allele_mismatch[1:which(cumulative_IBD==max(cumulative_IBD))]) # total number of allele mismatches to the right
  number_of_markers_in_the_right_markov_chain=which(cumulative_IBD==max(cumulative_IBD)) # total number of markers traversed until end of right sharing
  
  numer_of_haplotype_switches_in_the_right_markov_chain=length(which(c(FALSE, tail(test_haplotype[1:which(cumulative_IBD==max(cumulative_IBD))[1]],-1) != head(test_haplotype[1:which(cumulative_IBD==max(cumulative_IBD))[1]],-1)))) # number of times the chain switches between the two test haplotypes
  
  
  Final_FH_score <- log(Max_left_IBD+Max_right_IBD) # FH_score is the log(Final_IBD_score)
  
  number_of_allele_mismatches_in_the_markov_chain=number_of_allele_mismatches_in_the_left_markov_chain+number_of_markers_in_the_right_markov_chain # Total number of allele mismatches
  number_of_markers_in_the_markov_chain=number_of_markers_in_the_left_markov_chain+number_of_markers_in_the_right_markov_chain #Total number of markers in the IBD sharing region
  
  numer_of_haplotype_switches_in_the_markov_chain <- numer_of_haplotype_switches_in_the_left_markov_chain+numer_of_haplotype_switches_in_the_right_markov_chain #Total number of haplotype switches

  snp_density_in_data_file=nrow(data_file)/(max(data_file$position_cM)-min(data_file$position_cM)) # number of markers in 1cM
  total_number_of_markers_in_data_file=nrow(data_file)
  total_cM_span_of_data_file=max(data_file$position_cM)-min(data_file$position_cM)

  IBD_report=paste(Final_FH_score,Max_left_IBD,Max_right_IBD,(Max_left_cM+Max_right_cM),Max_left_cM,Max_right_cM,number_of_allele_mismatches_in_the_markov_chain,number_of_markers_in_the_markov_chain,numer_of_haplotype_switches_in_the_markov_chain,snp_density_in_data_file,total_number_of_markers_in_data_file,total_cM_span_of_data_file,sep='\t')
  
  print("IBD results for the disease-test pair is")
  print(c(Final_FH_score,Max_left_IBD,Max_right_IBD,(Max_left_cM+Max_right_cM),Max_left_cM,Max_right_cM,number_of_allele_mismatches_in_the_markov_chain,number_of_markers_in_the_markov_chain,numer_of_haplotype_switches_in_the_markov_chain,snp_density_in_data_file,total_number_of_markers_in_data_file,total_cM_span_of_data_file))
  #return all the required details of the IBD sharing to analyze later
  
  return(IBD_report)
  
}













