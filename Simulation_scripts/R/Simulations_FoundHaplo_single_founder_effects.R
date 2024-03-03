#' This is a wrapper around all other functions to calculate FH IBD score in simulation to store them in a text file
#' Add the relevant paths inside the function to load other related R scrips.
#' bi-allelic MAF frequency taken from gnomAD data (EUR) from ANNOVAR directory
#' Using hg19 genome built
#' Simulating single founder effects and comparing the FH scores. FH score is the log(sharing_LRT) in the output, which is the column 20 of the output

#' seed_id= Integer to identify the founder scenario, values take 1-10
#' sim_ID= Integer to track each simulation (1-1320)
#' RE_loci= Name of the disease-causing variant of interest, i.e. FAME1.chr8.119379052. Use the OMIM abbreviation for the disease. Mostly repeat expansion disease are used in the simulation.
#' sharing_w= Expected sharing length in total around the disease variant (0.5,1,2 and 5 in cM)
#' gen_error= 1 if genotype errors should be simulated, 0 otherwise
#' S.Er= 1 if switch errors should be simulated, 0 otherwise
#' IBD.version= 1 for simulating and running FoundHaplo for cases and 1.1 for for simulating and running FoundHaplo for conttrols
#' gen_error.rate= Genotype and imputation error rate allowed; default is 0.01 (1%)
#' S.Er.rate= switch error rate allowed; default is 20.05 switches per mbp
#' path.to.save= Directory path to save the output of FoundHaplo IBD sharing for further analysis
#' meiosis= Estimated number of meiosis between disease-test pair; default is 1
#' root= The folder to 1000 Genomes VCF files in .RDS format. save all 1000 Genomes VCF files by disease_variant (example: FAME1.chr8.vcf.RDS) in RDS format in a single folder.
#' path_MIS_ref= The folder to MIS-1000 Genomes reference panel, which is used to filter in genetic markers based on 1000 Genomes reference panel in the Michigan Imputation Server. Make sure the reference files are in chr.vcf format for all the chromosomes
#' path_hapmap= Directory path to genetic_map_HapMapII_GRCh37 files, which are also FoundHaplo/input_files/public_data/genetic_map_HapMapII_GRCh37/
#' path_gnomad_frq= The path to .txt file that has gnomAD EUR frequency. path_gnomad_frq would be in "/path/hg19_EUR.sites.2015_08.txt" format (This MAF file was taken from ANNOVAR)


Simulations_FoundHaplo_single_founder_effects=function(seed_id,sim_ID,RE_loci,sharing_w,gen_error,S.Er,IBD.version,gen_error.rate,S.Er.rate,path.to.save,meiosis,root,path_MIS_ref,path_hapmap,path_gnomad_frq)
{
  ## Add the relevant paths
   source("/path/Simulations_loading_inputfiles.R")
   source("/path/Simulations_create_Hap_dataset.R")
   source("/path/Simulations_convert_bp_cM.R")
   source("/path/Simulations_add_MAF.R")
   source("/path/Simulations_geno_Error.R")
   source("/path/Simulations_sw_Er.R")
   source("/path/Simulations_calculate_IBD.R")
   source("/path/Simulations_moveme.R")
   source("/path/Simulations_force_sharing_disease_haplotypes.R")
   source("/path/Simulations_force_sharing_cases.R")
  ## Add the relevant paths
  

  
  
  library(tibble)
  library(pryr)
  library(data.table)
  library(taRifx)
  library(dplyr)
  library(stringi)
  library(stringr)
  library(R.utils)
  library(vcfR)
  library(foreach)
  library(doParallel)
  library(TeachingDemos)
  

  RE=RE_loci
  w=sharing_w
  ID=sim_ID
  

#  root="path/RDS_files"
  RE_adjusted=strsplit(RE, ".",fixed=TRUE)
  RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
  
  path1=paste0(root,"/",RE_adjusted[1],".",RE_adjusted[2],".vcf.RDS")
  RDS_file=Simulations_loading_inputfiles(path1) # load relevant .vcf.RDS
  RDS_Hap=Simulations_create_Hap_dataset(RDS_file) # Separate haplotypes into two columns
  
  ### replace the second haplotype of X chromosome on all the males
  RDS_Hap_IND=RDS_Hap[,10:ncol(RDS_Hap)]
  
  RDS_Hap_IND_males= RDS_Hap_IND[,colSums(RDS_Hap_IND=="",na.rm = T)/nrow(RDS_Hap_IND) > 0.5 | colSums(is.na(RDS_Hap_IND))/nrow(RDS_Hap_IND) > 0.5]
  col_males=which(colnames(RDS_Hap) %in% colnames(RDS_Hap_IND_males))
  RDS_Hap[,col_males]=RDS_Hap[,col_males-1]
  ### replace the second haplotype of X chromosome on all the males
  
  #### Reference sites from Michigan Imputation server
  # Filter in genetic markers based on 1000 Genomes reference panel in the Michigan Imputation Server
  # Make sure the reference files are in chr.vcf format for all the chromosomes
  
  ref1=fread(paste0(path_MIS_ref,"/",RE_adjusted[2],".vcf"),skip = "#CHROM")
  
  RDS_Hap=subset(RDS_Hap,RDS_Hap$POS %in% ref1$POS )
  #### Reference sites from Michigan Imputation server
  
  
  RDS_Hap=Simulations_convert_bp_cM(RDS_Hap,RE,built.pos="position_hg19",path_hapmap) # Add cM column based on HapMap genetic maps and linear interpolation in R

  RDS_Hap=Simulations_add_MAF(RDS_Hap,RE,built.pos="position_hg19",path_gnomad_frq)  #Annotate gnomAD EUR frequency
  RDS_Hap_original=RDS_Hap
  
  # simulate a founder, 10 disease haplotypes and 5 cases to share a region around the disease variant with each disease haplotype. 
  simulation_samples= colnames(RDS_file@gt)[2:length(colnames(RDS_file@gt))]
  set.seed(seed_id)
  founder_DH=sample(simulation_samples,61,replace = FALSE)
  hap=sample(c("a","b"),61,replace = TRUE)
  
  founder_DH=paste0(founder_DH,"_",hap)
  founder=founder_DH[1]
  DH=founder_DH[2:61]
  
  DH=matrix(DH, 10, 6, byrow=TRUE)
  DH_mutiple=DH
  
  DH_mutiple=cbind(rep(founder,nrow(DH_mutiple)),DH_mutiple)
 # simulate a founder, 10 disease haplotypes and 5 cases to share a region around the disease variant with each disease haplotype. 
  
  #simulated cases ready with  founder and disease haplotypes
  
  # create a list with both haplotypes for the selected samples
  
  simulated_set=strsplit(DH_mutiple,"_",fixed=TRUE)
  simulated_set=sapply(simulated_set, "[[", 1)
  simulated_set=as.vector(paste(simulated_set))
  
  simulated_set_a=paste0(simulated_set,"_","a")
  simulated_set_b=paste0(simulated_set,"_","b")
  
  z=data.frame(simulated_set_a,simulated_set_b)
  z
  z[] <- lapply(z, as.character)
  a1=vector("character",(2*nrow(z)))
  k=1
  for (i in 1:nrow(z))
  {
    for(j in 1:2)
    {
      a1[k]=z[i,j]
      k=k+1
    }
  }
  a1=unique(a1)  # a1 has the column names of founder, both haplotypes disease and simulated cases
  
  
  simulated_cases_set=strsplit(DH_mutiple[,3:7],"_",fixed=TRUE)
  simulated_cases_set=sapply(simulated_cases_set, "[[", 1)
  simulated_cases_set=as.vector(paste(simulated_cases_set))
  
  simulated_cases_set_a=paste0(simulated_cases_set,"_","a")
  simulated_cases_set_b=paste0(simulated_cases_set,"_","b")
  
  z=data.frame(simulated_cases_set_a,simulated_cases_set_b)
  z
  z[] <- lapply(z, as.character)
  a2=vector("character",(2*nrow(z)))
  k=1
  for (i in 1:nrow(z))
  {
    for(j in 1:2)
    {
      a2[k]=z[i,j]
      k=k+1
    }
  }
  
  a2=unique(a2) # a2 has the column names of both haplotypes of simulated cases
  
  
  simulated_cases_set=strsplit(DH_mutiple[,2],"_",fixed=TRUE)
  simulated_cases_set=sapply(simulated_cases_set, "[[", 1)
  simulated_cases_set=as.vector(paste(simulated_cases_set))
  
  simulated_cases_set_a=paste0(simulated_cases_set,"_","a")
  simulated_cases_set_b=paste0(simulated_cases_set,"_","b")
  
  z=data.frame(simulated_cases_set_a,simulated_cases_set_b)
  z
  z[] <- lapply(z, as.character)
  a3=vector("character",(2*nrow(z)))
  k=1
  for (i in 1:nrow(z))
  {
    for(j in 1:2)
    {
      a3[k]=z[i,j]
      k=k+1
    }
  }
  
  a3=unique(a3)
  # a3 has the column names of all disease haplotypes
  
  
  simulated_set=data.frame(RDS_Hap[,c("CHROM","POS","POS.cM","frq")],RDS_Hap[,a1]) # taking only the simulated founder, disease haplotypes and cases
  
  #########Force sharing the founder's genomic region around the selected disease variant into selected disease haplotypes
  
  DH_mutiple_sim_ID=cbind(seed_id,sim_ID,DH_mutiple)
  
  simulate_founder_DH <- function(y)
  {
    dis.h=strsplit(DH_mutiple_sim_ID[y,][3],"_",fixed=TRUE)
    dis.h=sapply(dis.h, "[[", 1)
    dis.h.a=paste0(dis.h,"_a")
    dis.h.b=paste0(dis.h,"_b")
    
    vec=strsplit(DH_mutiple_sim_ID[y,][4],"_",fixed=TRUE)
    vec=sapply(vec, "[[", 1)
    vec=as.vector(paste(vec))
    
    vec=paste0(vec,"_a")
    
    vec=which(colnames(simulated_set) %in% vec)
    foreach(x=vec) %do% {
      
      cases_sim_data<-data.frame(simulated_set[,c("CHROM","POS","POS.cM","frq",dis.h.a,dis.h.b,colnames(simulated_set)[x],colnames(simulated_set)[x+1])])
      
      DH_list=DH_mutiple_sim_ID[,4]
      rate_list=c(100/(10/2),100/(9/2),100/(8/2),100/(7/2),100/(6/2),100/(5/2),100/(4/2),100/(3/2),100/(2/2),100/(1/2)) 
      # The distances in Morgans from the DCV locus to the left or right breakpoints of the simulated ancestral segment in the descendant are distributed as independent exponential random variables with rate specified by "rate_list" variable here.
      # We simulate ten instances by taking ten different rates as given by "rate_list" to represent the expected length of the ancestral segments to be 1,2,3,4,5,6,7,8,9,10cM
     
      DH_rate=data.frame(cbind(DH_list,rate_list))
      DH_rate=subset(DH_rate,DH_rate$DH_list==DH_mutiple_sim_ID[y,4])
      tmp=char2seed(paste0(DH_mutiple_sim_ID[y,3],DH_mutiple_sim_ID[y,4],"f"),set=FALSE) 
      set.seed(tmp)       # create a temporary seed based on founder, disease haplotype combination
      rate=as.numeric(as.character(DH_rate$rate_list))
      x_L_founder=rexp(rate=rate,1)*100 # simulating left arm
      x_R_founder=rexp(rate=rate,1)*100 # simulating right arm
      
      chr=strsplit(as.character(RE),".",fixed=TRUE)
      chr=sapply(chr, "[[", 2)
      path=paste0(path_hapmap,"/genetic_map_GRCh37_",chr,".txt") # specify path to HapMap genetic maps
      recom.map=read.delim(path)

      
      RE_adjusted=strsplit(RE, ".",fixed=TRUE)
      RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
      RE_adjusted=as.numeric(RE_adjusted[3])
      
      fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean) # linear interpolation to predict cM values
      RE.cM=fun(RE_adjusted)
      
      if(RE.cM-min(cases_sim_data$POS.cM)<x_L_founder){x_L_founder=RE.cM-min(cases_sim_data$POS.cM)} # x_L_founder and x_R_founder can not exceed chromosomal end points
      if(max(cases_sim_data$POS.cM)<RE.cM+x_R_founder){x_R_founder=max(cases_sim_data$POS.cM)-RE.cM} # x_L_founder and x_R_founder can not exceed chromosomal end points

      
      cases_sim_data=Simulations_force_sharing_disease_haplotypes(x_L_founder,x_R_founder,sharing_w=(x_L_founder+x_R_founder),RE,cases_sim_data,DH_mutiple_sim_ID[y,][3:9],gen_error=1,S.Er=0,gen_error.rate=2.5 * 10^(-8),S.Er.rate,built.pos="position_hg19",path_hapmap)

      cases_sim_data[,5]=as.character(cases_sim_data[,5])
      cases_sim_data[,6]=as.character(cases_sim_data[,6])
      cases_sim_data[,7]=as.character(cases_sim_data[,7])
      cases_sim_data[,8]=as.character(cases_sim_data[,8])
      cases_sim_data
      
    }
  }
  
  output_simulate_founder_DH_sim_data<-sapply(1:10,simulate_founder_DH)
  pull_output_founder_DH_sim_data <- function(x) x[,7:8]
  x1=lapply(output_simulate_founder_DH_sim_data,pull_output_founder_DH_sim_data)
  
  output_founder <- as.data.frame(do.call(cbind,x1)) # simulated subset of the original data
  
  ##replace RDS_Hap with simulated data
  replace_RDS_Hap <- function(y)
  {
    RDS_Hap[,colnames(output_founder)[y]]<- output_founder[,colnames(output_founder)[y]]
    
  }
  output_RDS_Hap<-sapply(1:ncol(output_founder),replace_RDS_Hap)
  colnames(output_RDS_Hap)=colnames(output_founder)
  
  RDS_Hap_dummy=RDS_Hap[,!colnames(RDS_Hap) %in% colnames(output_founder)]
  RDS_Hap=data.frame(cbind(RDS_Hap_dummy,output_RDS_Hap))
  

  simulation_table=data.frame(matrix(ncol =33, nrow = 0))
  
  # Simulation for diseae haplotypes is completed....
  
  if(IBD.version=="1"){
    
    #########Force sharing the disease haplotype's genomic region around the selected disease variant into selected cases
    #Process is similar to the previous function
    simulated_set=data.frame(RDS_Hap[,c("CHROM","POS","POS.cM","frq")],RDS_Hap[,a1])
    
    simulate_cases <- function(y,sharing_w)
    {
      
      
      dis.h=strsplit(DH_mutiple_sim_ID[y,][4],"_",fixed=TRUE)
      dis.h=sapply(dis.h, "[[", 1)
      dis.h.a=paste0(dis.h,"_a")
      dis.h.b=paste0(dis.h,"_b")
      ########
      
      vec=strsplit(DH_mutiple_sim_ID[y,][5:9],"_",fixed=TRUE)
      vec=sapply(vec, "[[", 1)
      vec=as.vector(paste(vec))
      
      vec=paste0(vec,"_a")
      
      vec=which(colnames(simulated_set) %in% vec)
      foreach(x=vec) %do% {
        
        cases_sim_data<-data.frame(simulated_set[,c("CHROM","POS","POS.cM","frq",dis.h.a,dis.h.b,colnames(simulated_set)[x],colnames(simulated_set)[x+1])])
        cases_sim_data[,5]=as.character(cases_sim_data[,5])
        cases_sim_data[,6]=as.character(cases_sim_data[,6])
        cases_sim_data[,7]=as.character(cases_sim_data[,7])
        cases_sim_data[,8]=as.character(cases_sim_data[,8])
        
        
        tmp <- char2seed(paste0(sim_ID,sharing_w,colnames(cases_sim_data)[5],colnames(cases_sim_data)[7]),set=FALSE)
        set.seed(tmp)
        
        w=sharing_w
        x_L=rexp(rate=200/w,1)*100
        x_R=rexp(rate=200/w,1)*100
        
        #####
        chr=strsplit(as.character(RE),".",fixed=TRUE)
        chr=sapply(chr, "[[", 2)
        path=paste0(path_hapmap,"/genetic_map_GRCh37_",chr,".txt")
        recom.map=read.delim(path)
        
        RE_adjusted=strsplit(RE, ".",fixed=TRUE)
        RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
        
        RE_adjusted=as.numeric(RE_adjusted[3])
        
        fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
        RE.cM=fun(RE_adjusted)
        
        if(RE.cM-min(cases_sim_data$POS.cM)<x_L){x_L=RE.cM-min(cases_sim_data$POS.cM)}
        if(max(cases_sim_data$POS.cM)<RE.cM+x_R){x_R=max(cases_sim_data$POS.cM)-RE.cM}
        #####
        
        cases_sim_data=Simulations_force_sharing_cases(x_L,x_R,w,RE,cases_sim_data,DH_mutiple_sim_ID[y,4:9],gen_error=1,S.Er=1,gen_error.rate=0.01,S.Er.rate,built.pos="position_hg19",path_hapmap)
        cases_sim_data[,5]=as.character(cases_sim_data[,5])
        cases_sim_data[,6]=as.character(cases_sim_data[,6])
        cases_sim_data[,7]=as.character(cases_sim_data[,7])
        cases_sim_data[,8]=as.character(cases_sim_data[,8])
        cases_sim_data

      }
    }
    
    output_cases_sim_data<-sapply(1:10,simulate_cases,sharing_w)
    pull_output_cases_sim_data <- function(x) x[,7:8]
    x1=lapply(output_cases_sim_data,pull_output_cases_sim_data)
    output <- as.data.frame(do.call(cbind,x1))
    
    
    pull_output_cases_sim_data <- function(x) x[,5:6]
    x1=lapply(output_cases_sim_data,pull_output_cases_sim_data)
    
    output_DH <- as.data.frame(do.call(cbind,x1))
    
    output_DH=output_DH[,as.numeric(rbind(seq(1,100,10), t2=seq(2,100,10)))]
    
    a<-seq(1,ncol(output_DH),2)
    b<-ncol(output)/2 # Or some other number
    a<-sapply(a, function (x) rep(x,b))
    list1=as.vector(a)
    
    a=seq(1,ncol(output),2)
    list2=rep(a,ncol(output_DH)/2)
    
    ###
    list_total=paste0(list1,"_",list2)
    
    save_FH_score_test=function(z)  
    {
      print(z)
      y=as.numeric(sapply(strsplit(z,"_",fixed=TRUE),"[[", 1))
      x=as.numeric(sapply(strsplit(z,"_",fixed=TRUE),"[[", 2))
      
      sim_data<-data.frame(RDS_Hap[,c("CHROM","POS","POS.cM","frq")],output_DH[,c(y,y+1)],output[,c(x,x+1)])
      
      tmp <- char2seed(paste0(sim_ID,sharing_w,colnames(sim_data)[5],colnames(sim_data)[7]),set=FALSE) # identify the relevant seed to save
      
      set.seed(tmp)
      x_L=rexp(rate=200/sharing_w,1)*100
      x_R=rexp(rate=200/sharing_w,1)*100
      
      #####
      chr=strsplit(as.character(RE),".",fixed=TRUE)
      chr=sapply(chr, "[[", 2)
      path=paste0(path_hapmap,"/genetic_map_GRCh37_",chr,".txt")
      recom.map=read.delim(path)
      
      RE_adjusted=strsplit(RE, ".",fixed=TRUE)
      RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
      
      RE_adjusted=as.numeric(RE_adjusted[3])
      
      fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
      RE.cM=fun(RE_adjusted)
      

      if(RE.cM-min(sim_data$POS.cM)<x_L){x_L=RE.cM-min(sim_data$POS.cM)}
      if(max(sim_data$POS.cM)<RE.cM+x_R){x_R=max(sim_data$POS.cM)-RE.cM}

      
      DH_list=DH_mutiple_sim_ID[,4]
      rate_list=c(100/(10/2),100/(9/2),100/(8/2),100/(7/2),100/(6/2),100/(5/2),100/(4/2),100/(3/2),100/(2/2),100/(1/2))
      DH_rate=data.frame(cbind(DH_list,rate_list))
      
      DH_rate=DH_rate %>% filter_all(any_vars(. %in%  colnames(output_DH)[y:(y+1)]))
      
      DH_mutiple_sim_ID_dummy=as.data.frame(DH_mutiple_sim_ID) %>% filter_all(any_vars(. %in%  colnames(output_DH)[y:(y+1)]))
      
      tmp=char2seed(paste0(DH_mutiple_sim_ID_dummy$V3,DH_mutiple_sim_ID_dummy$V4,"f"),set=FALSE) # identify the relevant seed to save
      set.seed(tmp)
      rate=as.numeric(as.character(DH_rate$rate_list))
      x_L_founder=rexp(rate=rate,1)*100
      x_R_founder=rexp(rate=rate,1)*100
      
      #####
      chr=strsplit(as.character(RE),".",fixed=TRUE)
      chr=sapply(chr, "[[", 2)
      path=paste0(path_hapmap,"/genetic_map_GRCh37_",chr,".txt")
      recom.map=read.delim(path)
      
      RE_adjusted=strsplit(RE, ".",fixed=TRUE)
      RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
      
      RE_adjusted=as.numeric(RE_adjusted[3])
      
      fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
      RE.cM=fun(RE_adjusted)

      
      if(RE.cM-min(sim_data$POS.cM)<x_L_founder){x_L_founder=RE.cM-min(sim_data$POS.cM)}
      if(max(sim_data$POS.cM)<RE.cM+x_R_founder){x_R_founder=max(sim_data$POS.cM)-RE.cM}
      #####
      
      sim_data1=sim_data
      rm(sim_data)

      ###############extra filtering before running IBD

      sim_data1[,5]=as.character(sim_data1[,5])
      sim_data1[,6]=as.character(sim_data1[,6])
      sim_data1[,7]=as.character(sim_data1[,7])
      sim_data1[,8]=as.character(sim_data1[,8])
      
      #get markers with unique POS.cM 
      
      table1=sim_data1[,c("frq","POS.cM","POS")]
      
      t1=distinct(table1,POS.cM,.keep_all= TRUE)
      attach(t1)
      t1 <- t1[order(POS),]
      detach(t1)
      
      
      if(length(unique(t1$POS.cM))!=nrow(t1)){ stop(" distinct POS.cM error")}
      
      sim_data1=subset(sim_data1,sim_data1$POS %in% t1$POS)
      
      DH_mutiple_sim_ID_dummy=as.data.frame(DH_mutiple_sim_ID) %>% filter_all(any_vars(. %in%  colnames(output_DH)[y:(y+1)]))
      
      sim_data1<-sim_data1[,c("POS","POS.cM","frq",as.character(DH_mutiple_sim_ID_dummy$V4),as.character(colnames(output)[x]),colnames(output)[x+1])]
  
      ###########extra filtering before running IBD

      sharing_value<-Simulations_calculate_IBD(sim_data1,RE,gen_error.rate=gen_error.rate,S.Er.rate,MA.cutoff=-0.4,meiosis=1,maf=0,path_hapmap)

      
      if(sapply(strsplit(colnames(sim_data1)[6],"_",fixed=TRUE),"[[",1) %in% sapply(strsplit(as.matrix(DH_mutiple_sim_ID_dummy)[5:9],"_",fixed=TRUE),"[[",1))
      {
        sharing.value.v1<-cbind(seed_id,ID,RE,w,x_L,x_R,rate,x_L_founder,x_R_founder,as.matrix(DH_mutiple_sim_ID_dummy)[3],as.matrix(DH_mutiple_sim_ID_dummy)[4],as.matrix(DH_mutiple_sim_ID_dummy)[5],as.matrix(DH_mutiple_sim_ID_dummy)[6],as.matrix(DH_mutiple_sim_ID_dummy)[7],as.matrix(DH_mutiple_sim_ID_dummy)[8],as.matrix(DH_mutiple_sim_ID_dummy)[9],sharing_value,type="simulated_case")
        
      }else{sharing.value.v1<-cbind(seed_id,ID,RE,w,x_L,x_R,rate,x_L_founder,x_R_founder,as.matrix(DH_mutiple_sim_ID_dummy)[3],as.matrix(DH_mutiple_sim_ID_dummy)[4],as.matrix(DH_mutiple_sim_ID_dummy)[5],as.matrix(DH_mutiple_sim_ID_dummy)[6],as.matrix(DH_mutiple_sim_ID_dummy)[7],as.matrix(DH_mutiple_sim_ID_dummy)[8],as.matrix(DH_mutiple_sim_ID_dummy)[9],sharing_value,type="truth")
      
      }
      
      test.ind=strsplit(colnames(sim_data1)[5],"_",fixed=TRUE)
      test.ind=sapply(test.ind, "[[", 1)
      
      return(sharing.value.v1)

    }
    
    
    a=lapply(list_total,save_FH_score_test)
    
    final_file=data.frame(do.call("rbind",a))
    dim(final_file)
    
    path1<-paste0(path.to.save,"/",ID,".txt")
    
    final_file=dplyr::distinct(final_file,disease.haplotype,Test.individual, .keep_all= TRUE)
    
    write.table(final_file,path1,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
  }
  
  
  #### running the model for controls
  
  else if(IBD.version=="1.1"){
    
    
    RDS_Hap$POS=as.numeric(paste(RDS_Hap$POS))
    RDS_Hap$POS.cM=as.numeric(paste(RDS_Hap$POS.cM))
    
    RDS_Hap=data.frame(RDS_Hap[,c("CHROM","POS","POS.cM","frq")],RDS_Hap[,11:ncol(RDS_Hap)])
     
    # Simulate 1% of error rate for controls 
    if(gen_error==1){
      
      sim_data=Simulations_geno_Error(RDS_Hap,gen_error.rate=0.01)
      
    }

    rm(RDS_Hap)
  
    
    save_FH_score_test=function(z)
    {
      print(z)
      y=as.numeric(sapply(strsplit(z,"_",fixed=TRUE),"[[", 1))
      x=as.numeric(sapply(strsplit(z,"_",fixed=TRUE),"[[", 2))
      
      sim_data1=sim_data
      DH_list=DH_mutiple_sim_ID[,4]
      rate_list=c(100/(10/2),100/(9/2),100/(8/2),100/(7/2),100/(6/2),100/(5/2),100/(4/2),100/(3/2),100/(2/2),100/(1/2))
      DH_rate=data.frame(cbind(DH_list,rate_list))
      
      DH_rate=DH_rate %>% filter_all(any_vars(. %in%  colnames(sim_data1)[y:(y+1)]))
      
      DH_mutiple_sim_ID_dummy=as.data.frame(DH_mutiple_sim_ID) %>% filter_all(any_vars(. %in%  colnames(sim_data1)[y:(y+1)]))
      
      tmp=char2seed(paste0(DH_mutiple_sim_ID_dummy$V3,DH_mutiple_sim_ID_dummy$V4,"f"),set=FALSE)
      set.seed(tmp)
      rate=as.numeric(as.character(DH_rate$rate_list))
      x_L_founder=rexp(rate=rate,1)*100
      x_R_founder=rexp(rate=rate,1)*100
      
      #####
      chr=strsplit(as.character(RE),".",fixed=TRUE)
      chr=sapply(chr, "[[", 2)
      path=paste0(path_hapmap,"/genetic_map_GRCh37_",chr,".txt")
      recom.map=read.delim(path)
      
      RE_adjusted=strsplit(RE, ".",fixed=TRUE)
      RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
      
      RE_adjusted=as.numeric(RE_adjusted[3])
      
      fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
      RE.cM=fun(RE_adjusted)
      
      
      if(RE.cM-min(sim_data1$POS.cM)<x_L_founder){x_L_founder=RE.cM-min(sim_data1$POS.cM)}
      if(max(sim_data1$POS.cM)<RE.cM+x_R_founder){x_R_founder=max(sim_data1$POS.cM)-RE.cM}
      #####
  
      if(x==y){
        # save a default value when comparing disease haplotype with it self default is 111111 for IBD value
        
        DH_mutiple_sim_ID_dummy=as.data.frame(DH_mutiple_sim_ID) %>% filter_all(any_vars(. %in%  colnames(sim_data1)[y:(y+1)]))
        disease.haplotype=as.character(DH_mutiple_sim_ID_dummy$V4)
        Test.individual=colnames(sim_data1)[x]
        Test.individual=sapply(strsplit(Test.individual,"_",fixed=TRUE),"[[",1)
        
        sharing_value<-matrix(c(RE,disease.haplotype,Test.individual,rep(111111,11)),nrow=1,ncol=length(c(RE,disease.haplotype,Test.individual,rep(111111,11))))
        colnames(sharing_value)=c("RE","disease.haplotype","Test.individual","sharing_LRT","stop.left.IBD","stop.right.IBD","sharing_left","sharing_right","num.gen.errors","num.markers.sharing","maf","snp.density","num.snps","cM_length")
        run_time=0
      }else{
        sim_data1<-data.frame(sim_data1[,c("CHROM","POS","POS.cM","frq")],sim_data1[,c(y,y+1)],sim_data1[,c(x,x+1)])
        if(S.Er==1){
          sim_data1=Simulations_sw_Er(sim_data1,colnames(sim_data1)[5],S.Er.rate)
          
        }

        DH_mutiple_sim_ID_dummy=as.data.frame(DH_mutiple_sim_ID) %>% filter_all(any_vars(. %in%  colnames(sim_data)[y:(y+1)]))
        disease.haplotype=as.character(DH_mutiple_sim_ID_dummy$V4)
        sim_data1<-sim_data1[,c("POS","POS.cM","frq", disease.haplotype,colnames(sim_data1)[7],colnames(sim_data1)[8])]
        
        sharing_value=Simulations_calculate_IBD(sim_data1,RE,gen_error.rate=gen_error.rate,S.Er.rate,MA.cutoff=-0.4,meiosis=1,maf=0,path_hapmap)

    colnames(sharing_value)=c("RE","disease.haplotype","Test.individual","sharing_LRT","stop.left.IBD","stop.right.IBD","sharing_left","sharing_right","num.gen.errors","num.markers.sharing","maf","snp.density","num.snps","cM_length")
        
        
      }
      
      sharing.value.v1<-cbind(seed_id,ID,RE,w,0,0,rate,x_L_founder,x_R_founder,as.matrix(DH_mutiple_sim_ID_dummy)[3],as.matrix(DH_mutiple_sim_ID_dummy)[4],as.matrix(DH_mutiple_sim_ID_dummy)[5],as.matrix(DH_mutiple_sim_ID_dummy)[6],as.matrix(DH_mutiple_sim_ID_dummy)[7],as.matrix(DH_mutiple_sim_ID_dummy)[8],as.matrix(DH_mutiple_sim_ID_dummy)[9],sharing_value,"control")
      
      
      
      test.ind=strsplit(colnames(sim_data1)[5],"_",fixed=TRUE)
      test.ind=sapply(test.ind, "[[", 1)
      
      sharing.value.v1

    }
    
    DH_list=DH_mutiple_sim_ID[,4]
    a<-which(colnames(sim_data) %in% paste0(sapply(strsplit(DH_list,"_",fixed=TRUE),"[[", 1),"_","a"))
    b<-length(5:ncol(sim_data))/2 # Or some other number
    a<-sapply(a, function (x) rep(x,b))
    list1=as.vector(a)
    
    a=seq(5,ncol(sim_data),2)
    list2=rep(a,length(which(colnames(sim_data) %in% paste0(sapply(strsplit(DH_list,"_",fixed=TRUE),"[[", 1),"_","a"))))
    
    
    list_total=paste0(list1,"_",list2)
    
    a=lapply(list_total,save_FH_score_test)
    
    final_file=data.frame(do.call("rbind",a))
    dim(final_file)
    
    path1<-paste0(path.to.save,"/",ID,".txt")
    
    final_file=dplyr::distinct(final_file,disease.haplotype,Test.individual, .keep_all= TRUE)
    
    write.table(final_file,path1,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
    
    
    
  }
  
}
