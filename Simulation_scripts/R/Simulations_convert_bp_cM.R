Simulations_convert_bp_cM=function(RDS_Hap,RE,built.pos="position_hg19",path_hapmap)
{

x=strsplit(as.character(RE),".",fixed=TRUE)
x=sapply(x, "[[", 2)


 if(built.pos=="position_hg19")
        {
             path=paste0(path_hapmap,"/genetic_map_GRCh37_",x,".txt")
    		 recom.map=read.delim(path)
        }else{
            stop("not hg19")
        }
        


Position.bp.=RDS_Hap[,"POS"]
Map.cM.=as.data.frame(matrix(0,ncol=1,nrow=nrow(RDS_Hap)))
file1=cbind(Position.bp.,Map.cM.)
colnames(file1)=c("Position.bp.","Map.cM.")



file1$Map.cM. <- recom.map$Map.cM.[match(file1$Position.bp., recom.map$Position.bp.)]
file1$Position.bp.=as.numeric(paste(file1$Position.bp.))

file2=subset(file1,file1$Map.cM.>0)

fun <- approxfun(recom.map$Position.bp.,recom.map$Map.cM.,ties=mean)
POS.cM=fun(file1$Position.bp.)

RDS_Hap=add_column(RDS_Hap,POS.cM, .after = 2)
RDS_Hap=subset(RDS_Hap,!is.na(RDS_Hap$POS.cM))

return(RDS_Hap)
                            

}


