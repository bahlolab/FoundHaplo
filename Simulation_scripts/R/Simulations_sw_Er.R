Simulations_sw_Er=function(sim_data,simulate_I,S.Er.rate=20.05)
{


  DB=simulate_I[1]
  x1=strsplit(DB,"_",fixed=TRUE)
  x1=sapply(x1, "[[", 1)

  n1=paste(x1,sep="_","a")


  DB.i=which(colnames(sim_data)==n1)
  
  if(identical(DB.i, integer(0))){DB.i=0
  }

library(LaplacesDemon)



C_l=max(sim_data$POS)
x=sim_data$POS[1]
switch=vector()
y=vector()
i=1
while(x<C_l)
{
  switch[i]=rexp(1,(1/(S.Er.rate/2)))

  x=x+switch[i]

  y[i]=x
  i=i+1
}


pos <-lapply(seq(1,(length(y)-1),by=2),
                   FUN = function(x){
                     c(y[x]:y[x+1])
                   }



)


positionsv1=unlist(pos)

positionsv1=round(positionsv1)
positions=which(sim_data$POS %in% positionsv1)


#swaping for all the individuals
swap_ind <-lapply(seq(5,ncol(sim_data),by=2),
                  FUN = function(x){
               
                    #
                    if(x==DB.i)
                    {
                      cbind.data.frame(sim_data[,x],sim_data[,(x+1)])

                    }

                    else{
                    C_l=max(sim_data$POS)
                    z=sim_data$POS[1]
                    switch=vector()
                    y=vector()
                    i=1
                    while(z<C_l)
                    {
                     switch[i]=rexp(1,(1/(S.Er.rate/2)))

                      z=z+switch[i]

                      y[i]=z
                      i=i+1
                    }


                    pos <-lapply(seq(1,(length(y)-1),by=2),
                                 FUN = function(z){
                                   c(y[z]:y[z+1])
                                 }



                    )


                    positionsv1=unlist(pos)

                    positionsv1=round(positionsv1)
                    positions=which(sim_data$POS %in% positionsv1)
                    #
                    sim_data[,x]=as.vector(sim_data[,x])
                    sim_data[,(x+1)]=as.vector(sim_data[,(x+1)])
                    p1 <- sim_data[,x][positions]
                    sim_data[,x][positions] <- sim_data[,(x+1)][positions]
                    sim_data[,(x+1)][positions] <- p1
                    cbind.data.frame(sim_data[,x],sim_data[,(x+1)])
                    }
                  }
)

#swapped data
output <- matrix(unlist(swap_ind), nrow = nrow(sim_data), ncol = (ncol(sim_data)-4), byrow = FALSE)

f1=cbind.data.frame(sim_data[,1:4],output)
colnames(f1)=colnames(sim_data)


return(f1)
}



