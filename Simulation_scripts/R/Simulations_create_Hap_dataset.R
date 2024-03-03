#' output is haplotype seperated ready to use file

Simulations_create_Hap_dataset=function(RDS_file)
{
  x1=RDS_file
  x2=x1@gt
  fix=x1@fix

  to=colnames(x2)[-1]
  x2=as.data.frame(x2[,-1])
  colnames(x2)=to

t1=data.frame(matrix(NA,nrow=nrow(x2),ncol=(ncol(x2)*2)))

n1=to
n2=to
n1=paste(n1,sep="_","a")
n2=paste(n2,sep="_","b")
z=data.frame(n1,n2)
z
z[] <- lapply(z, as.character)
a=vector("character",(2*nrow(z)))
k=1
for (i in 1:nrow(z))
{
  for(j in 1:2)
  {
    a[k]=z[i,j]
    k=k+1
  }
}

colnames(t1)=a



for(i in 1:ncol(x2)){
  t1[,(2*i-1)]=substr(x2[,i],0,1)
  t1[,(2*i)]=substr(x2[,i],3,3)
}

t2=as.data.frame(cbind(fix,t1))


#t2=remove.factors(t2)
t2=t2 %>% mutate_all(as.character)

t3=subset(t2,nchar(t2$REF)<2)
t4=subset(t3,nchar(t3$ALT)<2)
print(dim(t4))


frq=rep(0,nrow(t4))

t4=add_column(t4,frq, .after = 2)
return(t4)

}


#RDS_Hap=create_Hap_dataset(RDS_file)
