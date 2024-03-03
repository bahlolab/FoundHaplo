

#make sure CHROM is chr...
#path_gnomad_frq would be in "/path/humandb/hg19_EUR.sites.2015_08.txt" format


Simulations_add_MAF=function(RDS_Hap,RE,built.pos="position_hg19",path_gnomad_frq)
{
   
        if(built.pos=="position_hg19")
        {
            ALL=read.delim(path_gnomad_frq, stringsAsFactors=FALSE,header=FALSE) # give the path to EUR frequency file.txt
        }else{
            stop("not hg19")
        }
        
        
        chr=RDS_Hap$CHROM[1]
        ALL.chr=subset(ALL,ALL$V1 %in% chr)
        #filter biallelic SNPs
        ALL.chr=ALL.chr[!(duplicated(ALL.chr$V2) | duplicated(ALL.chr$V2, fromLast = TRUE)), ]
        colnames(ALL.chr)[2:4]=c("POS","REF","ALT")
        ALL.chr=subset(ALL.chr,nchar(ALL.chr$REF)<2 & nchar(ALL.chr$ALT)<2)
        ALL.chr=ALL.chr[,c(2,5)]
        colnames(ALL.chr)=c("POS","ALL")
        
        
        ALL.f=subset(ALL,ALL$V1 %in% chr)
        
        
        ALL.f=ALL.f[!(duplicated(ALL.f$V2) | duplicated(ALL.f$V2, fromLast = TRUE)), ]
        
        ALL.f=subset(ALL.f,nchar(ALL.f$V3)<2 & nchar(ALL.f$V4)<2)
        
        
        
        colnames(ALL.f)[2:4]=c("POS","REF","ALT")
        ALL.alleles=ALL.f[,c(2,3,4)]
        
        ALL.alleles.REF=ALL.alleles[,c(1,2)]
        colnames(ALL.alleles.REF)[2]="REF.1"
        
        ALL.alleles.ALT=ALL.alleles[,c(1,3)]
        colnames(ALL.alleles.ALT)[2]="ALT.1"
        
        
        RDS_Hap=RDS_Hap[!duplicated(RDS_Hap[c(2)]),]
        
        RDS_Hap$POS=as.numeric(RDS_Hap$POS)
        
        RDS_Hap=merge(RDS_Hap, ALL.chr, by="POS")
        RDS_Hap = merge(RDS_Hap, ALL.alleles.REF, by = "POS")
        RDS_Hap = merge(RDS_Hap, ALL.alleles.ALT, by = "POS")
        
        
        
        RDS_Hap=subset(RDS_Hap,!is.na(RDS_Hap$ALT.1))
        RDS_Hap=subset(RDS_Hap,!is.na(RDS_Hap$REF.1))
        #double check if the original REF/ALT is missing
        
        RDS_Hap=subset(RDS_Hap,RDS_Hap$REF==RDS_Hap$REF.1)
        RDS_Hap=subset(RDS_Hap,RDS_Hap$ALT==RDS_Hap$ALT.1)
        RDS_Hap=subset(RDS_Hap,!RDS_Hap$ALT.1 %in% 1)
        RDS_Hap$frq=RDS_Hap$ALL
        
        RDS_Hap=RDS_Hap[ , -which(names(RDS_Hap) %in% c("ALL","REF.1","ALT.1"))]
     
        attach(RDS_Hap)
        RDS_Hap <- RDS_Hap[order(POS),]
        detach(RDS_Hap)
        
       return(RDS_Hap)
    }
    

    
    
    
    
