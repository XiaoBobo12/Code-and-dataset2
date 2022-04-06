rm(list=ls())
result<-data.frame(matrix(nrow=9,ncol=8))
colnames(result)<-c("Allc.alpha","Interim.Time","LRT.nK.EP","LRT.ED","LRT.PES","MERT.nK.EP","MERT.ED","MERT.PES")
Allc.seq=c(2,1,3)
Inter.seq=c(2,1,3)
for(i in 4:6)
{
  for(j in 1:3)
  {
      print(paste("i=",i,"j=",j))
      if(j==2)
      {
        order.num=Inter.seq[i-3]
      }else
      {
        order.num=Inter.seq[i-3]*10+Allc.seq[j]
      }
      load(paste("C:/Users/pc/Desktop/Code and data/Table 3/Simulation Data/Example",order.num,".RData",sep=""))
      result$Allc.alpha[3*(i-4)+j]=alpha.allc[1]
      result$Interim.Time[3*(i-4)+j]=tau.vector[1]
      result$LRT.nK.EP[3*(i-4)+j]=paste(nc+nt,"(",simul.esti$Power,")",sep="")
      result$LRT.ED[3*(i-4)+j]=simul.esti$Events.average.simulation
      result$LRT.PES[3*(i-4)+j]=simul.esti$Power.cumulative[1]
  }
}
for(i in 1:3)
{
  for(j in 1:3)
  {
    print(paste("i=",i,"j=",j))
    if(j==2)
    {
      order.num=i
    }else
    {
      order.num=Inter.seq[i]*10+Allc.seq[j]
    }
    load(paste("C:/Users/pc/Desktop/Code and data/Table 3/Simulation Data/Example",order.num,".RData",sep=""))
    result$MERT.nK.EP[3*(i-1)+j]=paste(nc+nt,"(",simul.esti$Power,")",sep="")
    result$MERT.ED[3*(i-1)+j]=simul.esti$Events.average.simulation
    result$MERT.PES[3*(i-1)+j]=simul.esti$Power.cumulative[1]
  }
}
write.csv(x=result,file="D:/Table2.csv")