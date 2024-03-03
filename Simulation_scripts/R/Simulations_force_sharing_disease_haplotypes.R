Simulations_force_sharing_disease_haplotypes=function(x_L,x_R,sharing_w,RE_loci,RDS_Hap,simulate_I,gen_error,S.Er,gen_error.rate,S.Er.rate,built.pos="position_hg19",path_hapmap)
{
    
    
    print(dim(RDS_Hap))
    
    w=sharing_w
    RE=RE_loci
    
    x=strsplit(as.character(RE),".",fixed=TRUE)
    x=sapply(x, "[[", 2)
    
    
    if(built.pos=="position_hg19")
    {
        path=paste0(path_hapmap,"/genetic_map_GRCh37_",x,".txt")
        recom.map=read.delim(path)
    }
    else { stop("built is not hg19")
    }
    
    x=RDS_Hap
    
    print(dim(x))
    DBT=simulate_I
    
    
    
    x$POS=as.numeric(paste(x$POS))
    x$POS.cM=as.numeric(paste(x$POS.cM))
    print("start force sharing")
    print(dim(x))
    
    DB=DBT[1]
    
    T1=DBT[2][DBT[2] %in% colnames(x)]
    
    
    
    RE_adjusted=strsplit(RE, ".",fixed=TRUE)
    RE_adjusted=as.vector(unlist(RE_adjusted))
    
    RE_adjusted=as.numeric(RE_adjusted[3])
    
    
    
    fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
    RE.cM=fun(RE_adjusted)
    
    R1=RE.cM-x_L
    R2=RE.cM+x_R
    
    
    xL=which.max(x$POS.cM>=R1)
    xR=findInterval(R2,x$POS.cM)
    
    print(xL)
    print(xR)
    
    
    
    to=x

    to[xL:xR,T1]=to[xL:xR,DB]
    
    
    
    f=to
    
    
    
    
    
    if(gen_error==1){
        print("start error loop")
        
        f=Simulations_geno_Error(to,gen_error.rate)
        
        
    }
    
    if(S.Er==1){
        
        print("start switch")
        f1=Simulations_sw_Er(f,DBT,S.Er.rate)
        f=f1
        
    }
    
    
    return(f)
}

