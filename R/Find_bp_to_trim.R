#' Outputs the base pair positions corresponding to +-20cM from a given locus
#'
#' @description
#' Takes a text file with disease causing variants written in "FAME1.chr8.119379052." format, variants added in new line or a vector of variants in the same format and calculate base pair values corresponding to +-20cM from the disease variant locus
#' Genome built should be GRCh37 and genetic map in "dir_geneticMap" must have chromosomes named with "chr" prefix, add the prefix if needed
#' @param input_vector A vector with the variants (type \code{"character"})
#' @param dir_geneticMap Directory to genetic_map_HapMapII_GRCh37 location (type \code{"character"})
#' @param output_file Path to save the resutls (type \code{"character"})
#' @param size_to_trim_cM Size of the region allowed around the disease variant, default is 10cM either side (type \code{"numeric"})
#' @return Return a data frame with three columns 1.Original variant 2. base pair value of 10cM to the left of the variant 3. base pair value of 10cM to the right of the variant
#' @export
#' @examples
#' orig_dir <- getwd()
#' temp_dir <- tempdir()
#' setwd( temp_dir )
#' write.table(genetic_map_GRCh37_chr8,"genetic_map_GRCh37_chr8.txt",sep = "\t",quote=FALSE, row.names=FALSE,col.names = TRUE)
#' Find_bp_to_trim(input_vector=c("FAME1.chr8.119379052."),dir_geneticMap=temp_dir,output_file=paste0(temp_dir,"/DCV_bp.txt"))
#' setwd( orig_dir )



Find_bp_to_trim=function(input_vector,dir_geneticMap,output_file,size_to_trim_cM=10)
{

    positions=list()
    for(j in 1:length(input_vector))
    {

        DCV_adjusted=strsplit(input_vector[j], ".",fixed=TRUE)
        DCV_adjusted=as.vector(unlist(DCV_adjusted,recursive = FALSE))
        chr=DCV_adjusted[2]

        path=paste0(dir_geneticMap,"/genetic_map_GRCh37_",chr,".txt") # geneticMap files are in input_files/public_data/genetic_map_HapMapII_GRCh37
        recombination_map=read.delim(path)
        colnames(recombination_map)=c("Chromosome","position_bp","Rate.cM.Mb.","position_cM")

        fun_to_cM <- approxfun(recombination_map$position_bp,recombination_map$position_cM,ties=mean)
        DCV_cM=fun_to_cM(DCV_adjusted[3])


        DCV_cM_left=DCV_cM-size_to_trim_cM
        DCV_cM_right=DCV_cM+size_to_trim_cM


        temp=(recombination_map$position_cM>as.numeric(DCV_cM_left))
        POS_DCV_cM_left=table(temp)["FALSE"]+1
        POS_DCV_cM_left=POS_DCV_cM_left-1

        temp=(recombination_map$position_cM>as.numeric(DCV_cM_right))
        POS_RE_cM_right=table(temp)["FALSE"]+1


        a1=recombination_map$position_bp[POS_DCV_cM_left]
        a2=recombination_map$position_bp[POS_RE_cM_right]
        if(is.na(a1))
        {
            a1=min(recombination_map$position_bp)
            left_length=DCV_cM-fun_to_cM(a1)
            print(paste0("warning : Chromosome ends before +",size_to_trim_cM,"cM to the left"))
            print(paste0("allowing only ",left_length,"cM to the left"))
        }
        if(is.na(a2))
        {
            a2=max(recombination_map$position_bp)
           right_length=fun_to_cM(a2)-DCV_cM
            print(paste0("warning : Chromosome ends before +",size_to_trim_cM,"cM to the right"))
            print(paste0("allowing only ",right_length,"cM to the right"))
        }

        positions[[j]]=c(a1,a2)
    }

    positions_start=sapply(positions, "[[", 1)
    positions_end=sapply(positions, "[[",2)

    final_file=as.data.frame(cbind(input_vector,positions_start,positions_end))

    write.table(final_file,output_file,sep = "\t",quote=FALSE, row.names=FALSE,col.names = FALSE)
    return(final_file)
}
