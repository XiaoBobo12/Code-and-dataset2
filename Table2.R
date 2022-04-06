rm(list=ls())
data.results<-data.frame(matrix(nrow=40,ncol=9))
colnames(data.results)<-c("Allc","Harzard.Ratio","MERT.t1.t2","TRUE.t1.t2","l.t1","l.t2","t1.t2.1.1","t1.t2.1.2","t1.t2.2.1")
for(i in 2:5)
{
  for(j in 1:5)
  {
    for(c in 1:10)
    {
      print(paste("i=",i,"j=",j,"c=",c))
      if(j==1)
      {
        if(c<10)
        {
          order.num=i*10+c
        }
        else
        {
          order.num=i*100+c
        }
      }
      else
      {
        order.num=i*100+(j-1)*10+c
      }
      load(paste("C:/Users/pc/Desktop/Code and data/Table 2/Simulation Data/Scenario",order.num,".RData",sep=""))
      if(j==3)
      {
        data.results[(i-2)*10+c,1]=Allc
        data.results[(i-2)*10+c,2]=theta
        data.results[(i-2)*10+c,3]=paste("(",simul.esti$Parameters$t1,",",simul.esti$Parameters$t2,")",sep="")
        if(t1.true==t2.true)
        {
          data.results[(i-2)*10+c,4]=t1.true
        }
        else
        {
          data.results[(i-2)*10+c,4]=paste("(",simul.esti$Parameters$t1.true,",",simul.esti$Parameters$t2.true,")",sep="")
        }
      }
      data.results[(i-2)*10+c,j+4]=paste(sample.est$maximum.sample.size,"(",simul.esti$Power,")",sep="")
    }
  }
}
write.csv(x=data.results,file="D:/Table2.csv")