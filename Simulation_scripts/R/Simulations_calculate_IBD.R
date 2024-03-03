
# Make sure to use recombination map to approximate the cM
Simulations_calculate_IBD=function(data.file,RE_loci,gen_error.rate=0.01,S.Er.rate=20.05,MA.cutoff=-0.4,meiosis=1,maf,path_hapmap)
{

    RE=RE_loci
    g1=gen_error.rate
    m=meiosis
    
    
    x=strsplit(as.character(RE),".",fixed=TRUE)
    x=sapply(x, "[[", 2)
    path=paste0(path_hapmap,"/genetic_map_GRCh37_",x,".txt")
    recom.map=read.delim(path)
    
  
    colnames(data.file)=c("POS","POS.cM","frq",colnames(data.file)[4:6])
    
    
    
    data.file[,1]=as.numeric(data.file[,1])
    data.file[,2]=as.numeric(data.file[,2])
    data.file[,3]=as.numeric(data.file[,3])
    data.file[,4]=as.numeric(as.character(data.file[,4]))
    data.file[,5]=as.numeric(as.character(data.file[,5]))
    data.file[,6]=as.numeric(as.character(data.file[,6]))
    
    data.file <- data.file[!is.na(data.file[,4]), ]
    data.file <- data.file[!is.na(data.file[,5]), ]
    data.file <- data.file[!is.na(data.file[,6]), ]
    
    data.file <- data.file[data.file[,4]<2, ]
    data.file <- data.file[data.file[,5]<2,]
    data.file <- data.file[data.file[,6]<2,]
    
    RE_adjusted=strsplit(RE, ".",fixed=TRUE)
    RE_adjusted=as.vector(unlist(RE_adjusted,recursive = FALSE))
    r=(data.file$POS>as.numeric(RE_adjusted[3]))
    xR=table(r)["FALSE"]+1
    xL=xR-1
    
    #
    RE_adjusted=as.numeric(RE_adjusted[3])
    
    fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
    RE.cM=fun(RE_adjusted)
    
    IBD0=0
    IBD1=0
    

    
    k=1
    a=5
    A=c(5,6)
    j=xL
    i=4
    
    
    IBD_T=vector(mode = "numeric",length=2)
    error=vector(mode = "numeric",length=2)
    order=vector(mode = "character",length=2)
    diff=vector(mode = "numeric",length=2)
    dist=vector(mode = "numeric",length=2)
    
    
    for(j in xL:1)
    {
        
        f1=data.file[j,3]
        if(is.na(f1))
        {
            next
            
        }
        
        cum_distance=data.file[j+1,"POS.cM"]-data.file[j,"POS.cM"]
        cum_distance=cum_distance/100
        
        #d=cumulative distance in Morgans
        d=cum_distance
        
        
        
        if(data.file[j,a]==data.file[j,i] && data.file[j,i]==0)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g1*f1 + g0*f0)^2
            
            
            
            tp1=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
            
        }
        
        else if(data.file[j,a]==data.file[j,i] && data.file[j,i]==1)
        {
            
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g0*f1 + g1*f0)^2
            
            tp1=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
            
            
        }
        
        else if(data.file[j,A[!A %in% a]]==data.file[j,i] && data.file[j,i]==0)
        {
            
            
            a=A[!A %in% a]
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g1*f1 + g0*f0)^2
            
            
            
            tp1=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
            
        }
        
        else if(data.file[j,A[!A %in% a]]==data.file[j,i] && data.file[j,i]==1)
        {
            
            
            a=A[!A %in% a]
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g0*f1 + g1*f0)^2
            
            tp1=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
        }
        
        else if(data.file[j,a]==1 && data.file[j,i]==0)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            
            tp0= (g1*f1 + g0*f0) * (g0*f1 + g1*f0)
            
            tp1=g0*g1*f1*r0 + g0*g1*f0*r0 + g1*f1*g0*f1*r1 + g1*f1*g1*f0*r1 + g0*f0*g0*f1*r1 + g0*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=1
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
            
        }
        
        else if(data.file[j,a]==0 && data.file[j,i]==1)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            
            tp0= (g0*f1 + g1*f0) * (g1*f1 + g0*f0)
            
            tp1=g0*r0*f1*g1 + g0*r0*f0*g1 + g0*g1*f1*f1*r1 + g0*f1*g0*f0*r1 + g1*f0*g1*f1*r1 + g1*f0*g0*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=1
            order[k]=a
            dist[k]=RE.cM-data.file[j,"POS.cM"]
            
        }
        
        IBD_T[k]=-2*(IBD0 - IBD1)
        
        
        
        
        if(k>2)
        {

            diff[k-1]=IBD_T[k]-IBD_T[k-1]
            
            # rollmean of single vector and single window
            froll=frollmean(diff, 100)
            FirstNegative<-which(froll<=MA.cutoff)[1]
            #make threshold a parameter
            
            if(!is.na(FirstNegative))
            {break}
        }
        
        
        k=k+1
        
        
    }
    
    
    
    stop.left.IBD=max(IBD_T)
    
    left.end=RE.cM-data.file[xL-(which(IBD_T==max(IBD_T))[length(which(IBD_T==max(IBD_T)))]-1)+1,"POS.cM"]
    
    num.gen.errors.left=sum(error[1:which(IBD_T==max(IBD_T))])
    left.end.markers=which(IBD_T==max(IBD_T))
    
    
    
    k=1
    j=xR
    IBD0=0
    IBD1=0
    
    
    IBD_T=vector(mode = "numeric",length=2)
    error=vector(mode = "numeric",length=2)
    order=vector(mode = "character",length=2)
    diff=vector(mode = "numeric",length=2)
    dist=vector(mode = "numeric",length=2)
    
    for(j in xR:(nrow(data.file)-1))
    {
        
        f1=data.file[j,3]
        if(is.na(f1))
        {
            next
            
        }
        
        cum_distance=data.file[j,"POS.cM"]-data.file[j-1,"POS.cM"]
        cum_distance=cum_distance/100
        
        #d=cumulative distance in Morgans
        d=cum_distance
        
        
        if(data.file[j,a]==data.file[j,i] && data.file[j,i]==0)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g1*f1 + g0*f0)^2
            
            
            
            tp1=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
            
        }
        
        else if(data.file[j,a]==data.file[j,i] && data.file[j,i]==1)
        {
            
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g0*f1 + g1*f0)^2
            
            tp1=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
            
            
        }
        
        else if(data.file[j,A[!A %in% a]]==data.file[j,i] && data.file[j,i]==0)
        {
            
            
            a=A[!A %in% a]
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g1*f1 + g0*f0)^2
            
            
            
            tp1=g1*r0*f1*g1 + g0*r0*f0*g0 + f1*g1*r1*f1*g1 + f1*g1*r1*f0*g0 + f0*g0*g1*f1*r1 + f0*g0*f0*g0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
            
        }
        
        else if(data.file[j,A[!A %in% a]]==data.file[j,i] && data.file[j,i]==1)
        {
            
            
            a=A[!A %in% a]
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            tp0= (g0*f1 + g1*f0)^2
            
            tp1=g0*r0*f1*g0 + g1*r0*f0*g1 + g0*g0*f1*f1*r1 + f1*g0*r1*f0*g1 + f0*g1*g0*f1*r1 + g1*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=0
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
        }
        
        else if(data.file[j,a]==1 && data.file[j,i]==0)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            
            tp0= (g1*f1 + g0*f0) * (g0*f1 + g1*f0)
            
            tp1=g0*g1*f1*r0 + g0*g1*f0*r0 + g1*f1*g0*f1*r1 + g1*f1*g1*f0*r1 + g0*f0*g0*f1*r1 + g0*f0*g1*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=1
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
            
        }
        
        else if(data.file[j,a]==0 && data.file[j,i]==1)
        {
            
            r1=1-exp(-m*d)
            r0=exp(-m*d)
            g0=1-g1
            f0=1-f1
            
            
            tp0= (g0*f1 + g1*f0) * (g1*f1 + g0*f0)
            
            tp1=g0*r0*f1*g1 + g0*r0*f0*g1 + g0*g1*f1*f1*r1 + g0*f1*g0*f0*r1 + g1*f0*g1*f1*r1 + g1*f0*g0*f0*r1
            
            
            IBD0=IBD0+log(tp0)
            
            IBD1=IBD1+log(tp1)
            
            error[k]=1
            order[k]=a
            dist[k]=data.file[j,"POS.cM"]-RE.cM
            
        }
        
        IBD_T[k]=-2*(IBD0 - IBD1)
        
        
        if(k>2)
        {

            diff[k-1]=IBD_T[k]-IBD_T[k-1]
            
            # rollmean of single vector and single window
            froll=frollmean(diff, 100)
            FirstNegative<-which(froll<=MA.cutoff)[1]
            #make threshold a parameter
            
            if(!is.na(FirstNegative))
            {break}
        }
        
        k=k+1
        
    }
    
    stop.right.IBD=max(IBD_T)
    
    right.end=data.file[xR+ which(IBD_T==max(IBD_T))[length(which(IBD_T==max(IBD_T)))] -2,"POS.cM"]-RE.cM
    
    
    num.gen.errors.right=sum(error[1:which(IBD_T==max(IBD_T))])
    right.end.markers=which(IBD_T==max(IBD_T))
    
   
    sharing_LRT <- stop.left.IBD+stop.right.IBD
    
    num.gen.errors=num.gen.errors.left+num.gen.errors.right
    num.markers.sharing=left.end.markers+right.end.markers
    
    
    sharing_left=left.end
    sharing_right=right.end
    

    
    disease.haplotype=colnames(data.file)[i]

    x1=unlist(strsplit(colnames(data.file)[a], "_",fixed=TRUE))
    Test.individual=x1[1]
    
    snp.density=nrow(data.file)/(max(data.file$POS.cM)-min(data.file$POS.cM))
    num.snps=nrow(data.file)
    total_cM_length=max(data.file$POS.cM)-min(data.file$POS.cM)


 sharing.value=cbind(RE,disease.haplotype,Test.individual,sharing_LRT,stop.left.IBD,stop.right.IBD,sharing_left,sharing_right,num.gen.errors,num.markers.sharing,maf,snp.density,num.snps,total_cM_length)
    

    return(sharing.value)
    
    
    
    
    
    
    
}












