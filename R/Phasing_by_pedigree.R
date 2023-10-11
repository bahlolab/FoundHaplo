#' Phase samples with a disease variant by pedigrees
#'
#' @description
#' Takes a VCF file with known disease individuals and create disease haplotypes based on trio/duo, related samples or recessive variants .
#' @details
#' Make sure the input_vcf chromosome column doesn't have the "chr" prefix.
#' @param input_vcf File path to a VCF file  (type \code{"character"})
#' @param output_DIR Directory to save the output VCF files. Must have a dedicated directory for the output (type \code{"character"})
#' @param sample_info_file File path to a tab delimited .txt file with sample names and type of phasing to be used. Each phasing is included in a new line, include sample names as in the VCF file in below mentioned order.
#' For the type "trio": sample_info_file should have, affected-offspring,affected-parent,unaffected-parent trio.
#' For the type "duo" or "related": sample_info_file should have, affected-sample_1,affected-sample_2 duo.
#' For the type "recessive", a single affected sample is phased around a recessive variant using haplotype homozygosity : sample_info_file should have, affected-sample, recessive.
#' The function works similarly for "duo" and "related".
#' @param n.cores Number of cores to parallelize (type \code{"numeric"}).
#' @return Write a valid new VCF file in output_DIR "h1" column contains the disease haplorype saved as homozygous genotypes.
#' @import vcfR
#' @import foreach
#' @import doParallel
#' @import parallel
#' @importFrom plyr join
#' @importFrom dplyr mutate_all
#' @importFrom stringr str_count
#' @export
#' @examples
#' orig_DIR <- getwd()
#' temp_DIR <- tempdir()
#' setwd(temp_DIR)
#' library(vcfR)
#' write.vcf(FAME1_disease_cohort,paste0(temp_DIR,"/","FAME1_disease_cohort.vcf.gz"))
#' sample_info=data.frame(rbind(c("HG00362_1_HG00362_1,HG00362_2_HG00362_2","duo"),c("NA11920_NA11920,Affected_parent_NA11920_Affected_parent_NA11920,Unaffected_parent_NA11920_Unaffected_parent_NA11920","trio"),c("HG00313_1_HG00313_1,HG00313_2_HG00313_2","duo")))
#' write.table(sample_info,paste0(temp_DIR,"/","sample_info.txt"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
#' Phasing_by_pedigree(input_vcf = paste0(temp_DIR,"/FAME1_disease_cohort",".vcf.gz"),
#'                          output_DIR = temp_DIR,
#'                          sample_info_file = paste0(temp_DIR,"/","sample_info.txt"),n.cores=1)
#' setwd( orig_DIR )

Phasing_by_pedigree=function(input_vcf,output_DIR,sample_info_file,n.cores)
{
  if(!file.exists(input_vcf)){stop("input_vcf does not exist")}
  if(!dir.exists(output_DIR)){stop("output_DIR does not exist")}
  if(!file.exists(sample_info_file)){stop("sample_info_file does not exist")}

  sample_info=read.delim(sample_info_file,header=FALSE)

  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "PSOCK"
  )

  print(my.cluster)

  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)

  #check if it is registered (optional)
  foreach::getDoParRegistered()
  ## [1] TRUE
  #how many workers are available? (optional)
  foreach::getDoParWorkers()

  system.time(foreach(line = 1:nrow(sample_info),.packages=c('foreach','dplyr','vcfR','FoundHaplo','stringi','stringr')) %dopar% {
    if(sample_info[line,2] == "trio")
    {
      affected.o=sapply(strsplit(sample_info[line,1],",",fixed="TRUE"),"[[",1)
      affected.p=sapply(strsplit(sample_info[line,1],",",fixed="TRUE"),"[[",2)
      unaffected.p=sapply(strsplit(sample_info[line,1],",",fixed="TRUE"),"[[",3)

      input_file_vcf=read.vcfR(file=input_vcf,verbose = FALSE)

      input_file=as.data.frame(cbind(input_file_vcf@fix,input_file_vcf@gt))

      input_file=input_file[, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",affected.o,affected.p,unaffected.p)]

      haplotype_file=Create_hap_VCF(input_file) # seperates haplotypes to two columns as :a and :b

      haplotype_file$POS=as.numeric(haplotype_file$POS)

      haplotype_file_1=data.frame(apply(haplotype_file, 2, function(x) as.numeric(as.character(x))))

      colnames(haplotype_file_1)=colnames(haplotype_file)

      haplotype_file_1$REF=haplotype_file$REF[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ALT=haplotype_file$ALT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$INFO=haplotype_file$INFO[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$CHROM=haplotype_file$CHROM[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FORMAT=haplotype_file$FORMAT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$QUAL=haplotype_file$QUAL[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ID=haplotype_file$ID[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FILTER=haplotype_file$FILTER[haplotype_file$POS %in% haplotype_file_1$POS]


      dummy_sum=haplotype_file_1[,seq(10,ncol(haplotype_file_1),2)]+haplotype_file_1[,seq(11,ncol(haplotype_file_1),2)]
      colnames(dummy_sum)=c("affected.o","affected.p","unaffected.p")
      haplotype_file_1=data.frame(cbind(haplotype_file_1,dummy_sum))

      # remove markers that can not be distinguished
      haplotype_file_1=haplotype_file_1[which(haplotype_file_1$affected.o!=1 | haplotype_file_1$affected.p!=1|haplotype_file_1$unaffected.p!=1),]

      haplotype_file_2 = haplotype_file_1[,!(names(haplotype_file_1) %in% c("affected.o","affected.p","unaffected.p"))]

      completeFun <- function(data, desiredCols) {
        completeVec <- complete.cases(data[, desiredCols])
        return(data[completeVec, ])
      }

      data_file=completeFun(haplotype_file_2, c(10:ncol(haplotype_file_2))) # remove rows with NA in haplotype columns

      h1=vector(mode = "numeric",length=2)
      h2=vector(mode = "numeric",length=2)

      for(j in 1:nrow(data_file))
      {
        if((data_file[j,10]==data_file[j,12] | data_file[j,10]==data_file[j,13]) & (data_file[j,11]==data_file[j,14] | data_file[j,11]==data_file[j,15])) {
          h1[j]=data_file[j,10] # disease haplotype of the offspring
          h2[j]=data_file[j,11] # other haplotype of the offspring
        } else if((data_file[j,11]==data_file[j,12] | data_file[j,11]==data_file[j,13]) & (data_file[j,10]==data_file[j,14] | data_file[j,10]==data_file[j,15])) {
          h1[j]=data_file[j,11]
          h2[j]=data_file[j,10]

        } else {  # when the disease haplotype alleles can not be distinguished
          h1[j]=NA
          h2[j]=NA
        }
      }

    }

    if(sample_info[line,2] %in% c("duo","related"))
    {
      sample_1=sapply(strsplit(sample_info[line,1],",",fixed="TRUE"),"[[",1)
      sample_2=sapply(strsplit(sample_info[line,1],",",fixed="TRUE"),"[[",2)

      input_file_vcf=read.vcfR(file=input_vcf,verbose = FALSE)

      input_file=as.data.frame(cbind(input_file_vcf@fix,input_file_vcf@gt))

      input_file=input_file[, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample_1,sample_2)]

      haplotype_file=Create_hap_VCF(input_file)

      haplotype_file$POS=as.numeric(haplotype_file$POS)

      haplotype_file_1=data.frame(apply(haplotype_file, 2, function(x) as.numeric(as.character(x))))

      colnames(haplotype_file_1)=colnames(haplotype_file)

      haplotype_file_1$REF=haplotype_file$REF[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ALT=haplotype_file$ALT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$INFO=haplotype_file$INFO[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$CHROM=haplotype_file$CHROM[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FORMAT=haplotype_file$FORMAT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$QUAL=haplotype_file$QUAL[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ID=haplotype_file$ID[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FILTER=haplotype_file$FILTER[haplotype_file$POS %in% haplotype_file_1$POS]

      dummy_sum=haplotype_file_1[,seq(10,ncol(haplotype_file_1),2)]+haplotype_file_1[,seq(11,ncol(haplotype_file_1),2)]
      colnames(dummy_sum)=c("sample_1","sample_2")
      haplotype_file_1=data.frame(cbind(haplotype_file_1,dummy_sum))
      # remove markers that can not be distinguished
      haplotype_file_1=haplotype_file_1[which(haplotype_file_1$sample_1!=1 | haplotype_file_1$sample_2!=1),]

      haplotype_file_2 = haplotype_file_1[,!(names(haplotype_file_1) %in% c("sample_1","sample_2"))]

      completeFun <- function(data, desiredCols) {
        completeVec <- complete.cases(data[, desiredCols])
        return(data[completeVec, ])
      }

      data_file=completeFun(haplotype_file_2, c(10:ncol(haplotype_file_2)))

      h1=vector(mode = "numeric",length=2)
      h2=vector(mode = "numeric",length=2)

      for(j in 1:nrow(data_file))
      {

        if((data_file[j,10]+data_file[j,11]+data_file[j,12]+data_file[j,13])==0) {
          h1[j]=0
          h2[j]=data_file[j,10]+data_file[j,11]-h1[j]
        } else if((data_file[j,10]+data_file[j,11]+data_file[j,12]+data_file[j,13])==1) {
          h1[j]=0
          h2[j]=data_file[j,10]+data_file[j,11]-h1[j]

        } else if((data_file[j,10]+data_file[j,11]+data_file[j,12]+data_file[j,13])==3) {
          h1[j]=1
          h2[j]=data_file[j,10]+data_file[j,11]-h1[j]

        } else if((data_file[j,10]+data_file[j,11]+data_file[j,12]+data_file[j,13])==4) {
          h1[j]=1
          h2[j]=data_file[j,10]+data_file[j,13]-h1[j]

        } else { print(j)
          h1[j]=NA
          h2[j]=NA
        }
      }
    }

    if(sample_info[line,2] == "recessive")
    {
      sample1=sample_info[line,1]

      input_file_vcf=read.vcfR(file=input_vcf,verbose = FALSE)

      input_file=as.data.frame(cbind(input_file_vcf@fix,input_file_vcf@gt))

      input_file=input_file[, c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT",sample1)]

      haplotype_file=Create_hap_VCF(input_file) # seperates haplotypes to two columns as :a and :b

      haplotype_file$POS=as.numeric(haplotype_file$POS)

      haplotype_file_1=data.frame(apply(haplotype_file, 2, function(x) as.numeric(as.character(x))))

      colnames(haplotype_file_1)=colnames(haplotype_file)

      haplotype_file_1$REF=haplotype_file$REF[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ALT=haplotype_file$ALT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$INFO=haplotype_file$INFO[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$CHROM=haplotype_file$CHROM[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FORMAT=haplotype_file$FORMAT[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$QUAL=haplotype_file$QUAL[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$ID=haplotype_file$ID[haplotype_file$POS %in% haplotype_file_1$POS]
      haplotype_file_1$FILTER=haplotype_file$FILTER[haplotype_file$POS %in% haplotype_file_1$POS]

      completeFun <- function(data, desiredCols) {
        completeVec <- complete.cases(data[, desiredCols])
        return(data[completeVec, ])
      }

      data_file=completeFun(haplotype_file_1, c(10:ncol(haplotype_file_1))) # remove rows with NA in haplotype columns

      h1=vector(mode = "numeric",length=2)
      h2=vector(mode = "numeric",length=2)

      for(j in 1:nrow(data_file))
      {
        if((data_file[j,10]==data_file[j,11])) {
          h1[j]=data_file[j,10] # disease haplotype
          h2[j]=data_file[j,10] # other haplotype
        }  else {  # when the disease haplotype alleles can not be distinguished
          h1[j]=NA
          h2[j]=NA
        }
      }

    }


    data_file$h1=h1

    data_file=data_file[!is.na(data_file$h1),]

    data_file$h1=paste0(data_file$h1, "/", data_file$h1)
    colnames(data_file)[1]="#CHROM"

    data_file[["ID"]][is.na(data_file[["ID"]])] <- "."
    data_file[["QUAL"]][is.na(data_file[["QUAL"]])] <- "."
    data_file[["FILTER"]][is.na(data_file[["FILTER"]])] <- "."

    #add dummy FORMAT
    number_of_FORMAT_fields=str_count(unique(data_file$FORMAT)[1], ":")
    dummy1=rep(paste(replicate(number_of_FORMAT_fields, ":100"), collapse = ""),nrow(data_file))
    dummy2=data_file$h1
    data_file$h1=paste0(dummy2,dummy1)


    write.table(input_file_vcf@meta,paste0(output_DIR,"/",sample_info[line,1],".vcf"),sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
    write.table(data_file,paste0(output_DIR,"/",sample_info[line,1],".vcf"),sep ="\t",quote=FALSE, row.names=FALSE,col.names = TRUE,append=TRUE)

    print(paste0("line is ",line))
  })
  parallel::stopCluster(my.cluster)
}

