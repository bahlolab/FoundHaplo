
Simulations_geno_Error=function(sim_data,gen_error.rate=0.01)
{


  Add.geno <-lapply(seq(5,ncol(sim_data),by=2),
                    FUN = function(x){
                      
                      #

                      r1=rbinom(nrow(sim_data),1,gen_error.rate/2)
                      table1=as.data.frame(cbind(seq(1:nrow(sim_data)),r1))
                      colnames(table1)=c("posi","rbinorm")
                      x2=table1$posi[table1$rbinorm==1]

                      t1=sim_data[,x][x2]
                      zero <- which(t1==0)
                      one <- which(t1==1)

                      t1[zero]=1
                      t1[one]=0

                      sim_data[,x][x2]=t1
                      ###################

                      r2=setdiff(seq(1:nrow(sim_data)),x2)

                      r1=rbinom(length(r2),1,gen_error.rate/2)
                      table1=as.data.frame(cbind(r2,r1))
                      colnames(table1)=c("posi","rbinorm")
                      x2=table1$posi[table1$rbinorm==1]

                      t1=sim_data[,(x+1)][x2]
    

                      zero <- which(t1==0)
                      one <- which(t1==1)

                      t1[zero]=1
                      t1[one]=0

                      sim_data[,(x+1)][x2]=t1


                      cbind.data.frame(sim_data[,x],sim_data[,(x+1)])
                    }

  )


  output <- matrix(unlist(Add.geno), nrow = nrow(sim_data), ncol = (ncol(sim_data)-4), byrow = FALSE)

  f1=cbind.data.frame(sim_data[,1:4],output)
  colnames(f1)=colnames(sim_data)


  return(f1)
}







