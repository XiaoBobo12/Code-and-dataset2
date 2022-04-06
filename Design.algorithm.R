library(markovchain)
##################################
## survival functions and hazards functions
##################################

hazard.function.control<-function(t,k,h0)
  ## the hazard function of the control group
  ## t:the survival time
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
{
  return(k*h0*(h0*t)^(k-1))
}

hazard.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the hazard function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t:the survival time
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(t1==t2)
  {
    h0t=k*h0*(h0*t)^(k-1)
    h1t=theta*h0t
    return((t<t1)*h0t+(t>=t1)*h1t)
  }
  else
  {
    h0t=k*h0*(h0*t)^(k-1)
    t.std=(t-t1)/(t2-t1)
    lt=0*(t<=t1)+pbeta(t.std,a,b)*(t1<t && t<=t2)+(t>t2)
    h1t=(1-lt+theta*lt)*h0t
  }
  return(h1t)
}

cum.hazard.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the cumulative hazard function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  cum.hazard.treatment.unit<-integrate(f=hazard.function.treatment,lower=0,upper=t,k,h0,t1,t2,theta,a,b)$value
  return(cum.hazard.treatment.unit)
}

cum.hazard.function.treatment.ref<-function(t,k,h0,t1,t2,theta,a,b,Cum.Hazd)
  ## the standardialized cumulative hazard function of the treatment group, which is used to produce the survival time
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  return(cum.hazard.function.treatment(t,k,h0,t1,t2,theta,a,b)-Cum.Hazd)
}

survival.function.control<-function(t,k,h0)
  ## the survival function of the control group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
{
  return(exp(-(h0*t)^k))
}

survival.function.treatment<-function(t,k,h0,t1,t2,theta,a,b)
  ## the survival function of the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(length(t)==1)
  {
    survival.proportion=exp(-(cum.hazard.function.treatment(t,k,h0,t1,t2,theta,a,b)))
    return(survival.proportion)
  }
  else
  {
    n=length(t)
    survival.proportion<-rep(NA,n)
    for(i in 1:n)
    {
      survival.proportion[i]=exp(-(cum.hazard.function.treatment(t[i],k,h0,t1,t2,theta,a,b)))
    }
    return(survival.proportion)
  }
}

survival.time.control<-function(n,k,h0)
  ## -n:number of observations. If length(n) > 1, the length is taken to be the number required
  ## -rate0:rates before the time of change
  ## -k:the weibull shape parameter
{
  U=runif(n)
  radnm=(-log(U))^(1/k)/h0
  return(radnm)
}

survival.time.treatment<-function(n,k,h0,t1,t2,theta,a,b)
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  if(theta==1)
  {
    U=runif(n)
    survival.time=(-log(U))^(1/k)/h0
  }
  else
  {
    U=runif(n)
    Cum.Hazd=-log(U)
    survival.time<-rep(NA,n)
    for(i in 1:n)
    {
      survival.time[i]=uniroot(f=cum.hazard.function.treatment.ref,interval=c(0,10000),k,h0,t1,t2,theta,a,b,Cum.Hazd[i])$root
    }
  }
  return(list(survival.time=survival.time,survival.proportion=U))
}


rpoipross<-function(n,rate)
  ## produce the entry time points accroding to the Poisson process
  ## -n:the number of time points
  ## -rate:intensity
{
  if(length(n)!=1)
    stop("n is not one parameter.")
  if(length(rate)!=1)
    stop("rate is not one parameter.")
  if(is.integer(n))
  {
    stop("n is not integer!")
  }
  if(n<=0)
  {
    stop("n is not bigger than 0")
  }
  if(!is.double(rate))
  {
    stop("rate is not a decimal or integer")
  }
  if(rate<=0)
  {
    stop("rate is not bigger than 0")
  }
  t=0
  S=rep(0,n)
  for(i in c(1:n))
  {
    U<-runif(1)
    t=t-log(U)/rate
    S[i]=t
  }
  return (S)
}
########################################################
#########Simulation algorithms##########################
########################################################
dataset.produce<-function(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1,t2,theta,a,b)
  ## produce the dataset used for the hypothesis in the Monte-Carlo simulation-procedure
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.num: the splitted number of an unit period for an Markov chain
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  datasetc<-matrix(nrow=nc,ncol=7)
  datasett<-matrix(nrow=nt,ncol=7)
  datasetc[,1]=rep(1,nc)## "1" represents the control group 
  datasett[,1]=rep(0,nt)## "0" represents the treatment group 
  h.loss.control=-log(1-loss.control)## the hazard constant of loss to follow-up in the control group
  h.loss.treatment=-log(1-loss.treatment)## the hazard constant of loss to follow-up in the treatment group
  ##produce the initial survival time of two groups
  datasetc[,2]=survival.time.control(nc,k,h0)
  datasett[,2]=survival.time.treatment(nt,k,h0,t1,t2,theta,a,b)$survival.time
  ##produce the recruitment time of two groups
  if(accrual.type[1]==0)
  {
    ##produce the entry time of control group and treatment group by poisson process
    entry=rpoipross(nc+nt,(nc+nt)/A)
    entry.total.order=c(1:(nc+nt))## original order number of two groups
    entry.control.order=sample(entry.total.order,nc)## select the order number corresponding to the entry time of the control group randomly
    entry.treatment.order=entry.total.order[-match(entry.control.order,entry.total.order)]## the rest is the order number of the treatment group
    entry.treatment.order=entry.treatment.order[sample(1:nt)]## permutate the order number of the treatment group
    datasetc[,3]=entry[entry.control.order]## get the entry number of the control group
    datasett[,3]=entry[entry.treatment.order]## get the entry number of the treatment group
  }
  else
  {
    a=accrual.type[2]
    datasetc[,3]=A*(runif(nc,0,1))^(1/a)
    datasett[,3]=A*(runif(nt,0,1))^(1/a)
  }
  ## produce the random censoring time
  datasetc[,4]=rexp(nc,h.loss.control)
  datasett[,4]=rexp(nt,h.loss.treatment)
  ## if the survial time is longer than the  censor time,the subject is censored;uncersored vice versa
  datasetc[,5]= datasetc[,2]<datasetc[,4]
  datasett[,5]= datasett[,2]<datasett[,4]
  ## obtain the observed time since the recruitment
  datasetc[,6]=apply(cbind(datasetc[,2],datasetc[,4]),1,min)
  datasett[,6]=apply(cbind(datasett[,2],datasett[,4]),1,min)
  ## obtain the observed time since the beginning of the clinical trial
  datasetc[,7]=datasetc[,6]+datasetc[,3]
  datasett[,7]=datasett[,6]+datasett[,3]
  dataset.whole=matrix(nrow=nc+nt,ncol=8)## the whole dataset of the control arm and the treatment arm 
  dataset.cb0=rbind(datasett,datasetc)## combine two datasets
  dataset.cb=dataset.cb0[order(dataset.cb0[,7]),]## sort the matrix by the column,which represents "the end time of actual visit"
  dataset.whole[,1]=seq(1,nc+nt)## the order number of two groups
  dataset.whole[,2:8]=dataset.cb##the data need to store
  colnames(dataset.whole)=c("order.number","group","survival.time","entry.time","censorring.time","censor.or.death","actual.observed.time.since.recruitment","actual.observed.time.since.begining")
  return(dataset.whole)
}

log.rank.test<-function(dataset.input)
  #This function is used to calculate the value of the log-rank test statistic
{
  dataset.analysis<-dataset.input[,c(2,6,7)]## group, censor.or.death, actual.observed.time.since.recruitment
  dataset.analysis<-dataset.analysis[order(dataset.analysis[,3]),]
  dataset.uncensored<-dataset.analysis[dataset.analysis[,2]==1,]## remove the censored patients
  number.uncensored<-length(dataset.uncensored[,1])## get the number of uncensored subject
  nc<-sum(dataset.analysis[,1]==1)
  ne<-sum(dataset.analysis[,1]==0)
  ncvector=rep(0,nc+ne+1)
  pcvector=rep(0,nc+ne+1)
  nevector=rep(0,nc+ne+1)
  pevector=rep(0,nc+ne+1)
  ncvector[1]=nc
  nevector[1]=ne
  pcvector[1]=nc/(nc+ne)
  pevector[1]=ne/(nc+ne)
  Sw=matrix(0,nrow=number.uncensored,ncol=3)
  ## the data used to calculate log-rank statistics 
  ## the first row represents the order number of the death after the onset of effect
  ## the second row represents the  numerator of the log-rank statistics
  ## the third row represents the square of the denumerator of the log-rank statistics
  iterator=0## an iterator used to calculate the log-rank statistics,which represents the number of deaths used in the calculation of log-rank test statistics
  for(i in 1:(nc+ne+1))
  {
    if(dataset.analysis[i,1]==1)## the number of the patiens in control group minus 1,while the other remains unchanged
    {
      ncvector[i+1]=ncvector[i]-1
      nevector[i+1]=nevector[i]
    }
    if(dataset.analysis[i,1]==0)## the number of the patiens in treatment group minus 1,while the other remains unchanged
    {
      ncvector[i+1]=ncvector[i]
      nevector[i+1]=nevector[i]-1
    }
    if(ncvector[i+1]==0 & nevector[i+1]==0)## once the number of either group is 0,stop the calculation of log-rank statistics
      break
    pcvector[i+1]=ncvector[i+1]/(ncvector[i+1]+nevector[i+1])
    pevector[i+1]=nevector[i+1]/(ncvector[i+1]+nevector[i+1])
    if(dataset.analysis[i,2]==1)
    {
      ## uncesored
      iterator=iterator+1##iterator add 1
      Sw[iterator,1]=iterator
      if(iterator!=1)
      {
        Sw[iterator,2]=Sw[iterator-1,2]+(dataset.analysis[i,1]- pcvector[i])
        Sw[iterator,3]=Sw[iterator-1,3]+pcvector[i]*pevector[i]
      }
      else
      {
        Sw[iterator,2]=(dataset.analysis[i,1]- pevector[i])
        Sw[iterator,3]=pcvector[i]*pevector[i]
      }
    }
  }
  ##get the value of the log-rank statistics corresponding to interim analysis and final analysis
  Sw_used=Sw[iterator,]
  Sw_value<-Sw_used[2]/sqrt(Sw_used[3])## calculate the value of the piecewise weitghted log-rank statistics
  return(list(Sw_value=Sw_value,death.all=number.uncensored,num.recru=(nc+ne),Sw.ustd=Sw_used[2],SW.Var.H0=Sw_used[3]))
}

Maximin.efficiency.robust.test<-function(dataset.input,t.low,t.upper)
  ## This function is used to calculate the value of the Maximin efficiency robust test
  ## dataset.input: the dataset used to survival analysis
  ## t.low: the minimum delay time
  ## t.upper: the maximum delay time
{
  dataset.analysis<-dataset.input[,c(2,6,7)]## group, censor.or.death, actual.observed.time.since.recruitment
  dataset.analysis<-dataset.analysis[order(dataset.analysis[,3]),]
  number.uncensored<-sum(dataset.analysis[,2]==1)## get the number of uncensored subject
  nc<-sum(dataset.analysis[,1]==1)
  ne<-sum(dataset.analysis[,1]==0)
  nc.tir<-nc
  ne.tir<-ne
  Calc.matrix<-matrix(nrow=nc+ne,ncol=11)
  colnames(Calc.matrix)<-c("pc","pe","pc.pe","psi","censor.or.not","t1_t2","greater.t2","Group","weight","Ustd.unit","Variance.unit")
  Calc.matrix[1,1]=nc.tir/(nc.tir+ne.tir)
  Calc.matrix[1,2]=ne.tir/(nc.tir+ne.tir)
  Calc.matrix[1,3]=nc.tir*ne.tir/(nc.tir+ne.tir)^2
  Calc.matrix[1,4]=Calc.matrix[1,3]*dataset.analysis[1,2]
  Calc.matrix[1,5]=dataset.analysis[1,2]
  Time.obs=dataset.analysis[1,3]
  Calc.matrix[1,6]=Time.obs>=t.low & Time.obs<=t.upper
  Calc.matrix[1,7]=Time.obs>t.upper
  Calc.matrix[1,8]=dataset.analysis[1,1]
  nc.tir=nc.tir-dataset.analysis[1,1]
  ne.tir=ne.tir-(1-dataset.analysis[1,1])
  for(i in 2:(nc+ne))
  {
    Calc.matrix[i,1]=nc.tir/(nc.tir+ne.tir)
    Calc.matrix[i,2]=ne.tir/(nc.tir+ne.tir)
    Calc.matrix[i,3]=nc.tir*ne.tir/(nc.tir+ne.tir)^2
    Calc.matrix[i,4]=Calc.matrix[i,3]*dataset.analysis[i,2]+Calc.matrix[i-1,4]
    Calc.matrix[i,5]=dataset.analysis[i,2]
    Time.obs=dataset.analysis[i,3]
    Calc.matrix[i,6]=Time.obs>=t.low & Time.obs<=t.upper
    Calc.matrix[i,7]=Time.obs>t.upper
    Calc.matrix[i,8]=dataset.analysis[i,1]
    nc.tir=nc.tir-dataset.analysis[i,1]
    ne.tir=ne.tir-(1-dataset.analysis[i,1])
  }
  Calc.matrix[,4]=Calc.matrix[,4]/(nc+ne)
  psi.tau<-Calc.matrix[nc+ne,4]
  if(sum(Calc.matrix[,6])==0)
  {
    psi.t1=0
    psi.t2=0
  }
  else
  {
    psi.t1<-ifelse(min(which(Calc.matrix[,6]==1))==1,0,Calc.matrix[min(which(Calc.matrix[,6]==1))-1,4])
    psi.t2<-Calc.matrix[max(which(Calc.matrix[,6]==1)),4]
  }
  Calc.matrix[,9]=sqrt((psi.tau-psi.t1)/(psi.tau-Calc.matrix[,4]*Calc.matrix[,6]))*Calc.matrix[,6]+2*sqrt((psi.tau-psi.t1)/(psi.tau-psi.t2))*Calc.matrix[,7]
  Calc.matrix[,10]=Calc.matrix[,9]*(Calc.matrix[,8]-Calc.matrix[,1])*Calc.matrix[,5]
  Calc.matrix[,11]=Calc.matrix[,9]^2*Calc.matrix[,3]*Calc.matrix[,5]
  Value.unstd<-sum(Calc.matrix[,10])
  var.H0<-sum(Calc.matrix[,11])
  Sw_value<-Value.unstd/sqrt(var.H0)
  return(list(Sw_value=Sw_value,death.all=number.uncensored,num.recru=(nc+ne),Sw.ustd=Value.unstd,SW.Var.H0=var.H0,Calc.matrix=Calc.matrix,psi.t2=psi.t2,psi.t1=psi.t1,psi.tau=psi.tau))
}

maximum.duration.trial.Monte.Carlo.simulation<-function(nc,nt,accrual.type,A,tau.vector,bound.list,loss.control,loss.treatment,k,h0,t1,t2,t1.true,t2.true,theta,a,b,sides,alpha,test.type,simnum)
  ## this function is used to simulate the maximum duration trial
  ## nc,nt:the sample sizes of the control and treatment groups
  ## accrual.type:select the way of how patients enter ths study
  ## c(0):patients enter the study according to the poisson process
  ## c(1,a):patients enter the study according to the distribution of F(x)=(x/A)^a
  ## A: the recruitment time(restricted to an integer)
  ## tau.vector: the calendar time of every analysis 
  ## bound.list: (lowerbound,upperbound) which is determined by the multivariate distribution
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used in the MRET test)
  ## t1.true,t2.true:the minimum and maximum delay time (used in the generation of the dataset)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time  
  ## sides: 1:right-sided test;2:two-sided test
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## simnum: the number of repetitions in the Monte-Carlo simulation procedure
{
  numberVisits=length(tau.vector)
  ## [lowerbound,upperbound] corresponding to the acceptance region
  lowerbound<-bound.list$lowerbound
  upperbound<-bound.list$upperbound
  test.statistics.simulation=matrix(0,nrow=simnum,ncol=numberVisits)
  Sw.ustd.simulation=matrix(0,nrow=simnum,ncol=numberVisits)
  SW.Var.H0.simulation=matrix(0,nrow=simnum,ncol=numberVisits)
  Decision.simulation=matrix(0,nrow=simnum,ncol=numberVisits)
  Events.number.total.simulation=matrix(0,nrow=simnum,ncol=numberVisits)
  Total.number.simulation=matrix(0,nrow=simnum,ncol=numberVisits)## store the cumulative number of enrolled patients corresponding to every analysis
  if(test.type==1)
  {
    ## log-rank test
    for(l in 1:simnum)
    {
      print(paste(l,"/",simnum))
      dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
      Sw_value=rep(0,numberVisits)
      for(i in 1:numberVisits)
      {
        analysis.time=tau.vector[i]
        dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
        dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
        dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
        dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
        test.statistic=log.rank.test(dataset.analysis)
        Sw_value[i]=test.statistic$Sw_value
        Sw.ustd.simulation[l,i]=test.statistic$Sw.ustd
        SW.Var.H0.simulation[l,i]=test.statistic$SW.Var.H0
        Events.number.total.simulation[l,i]=test.statistic$death.all
        Total.number.simulation[l,i]=test.statistic$num.recru
      }
      test.statistics.simulation[l,]=Sw_value
      reject.decision=!(Sw_value<=upperbound&lowerbound<=Sw_value)## the rejection situation for respective situation,1 represents rejection,0 represents acceptance
      for(i in 1:(numberVisits-1))## once the hypothesis is rejected in one interim ananlysis,the following interim analyses and final analysis won't be conducted,so let them be 0
      {
        if(reject.decision[i]==1)
        {
          for(j in c((i+1):numberVisits))
          {
            reject.decision[j]=0
            Events.number.total.simulation[l,j]=0
            Total.number.simulation[l,j]=0
          }
          break
        } 
      }
      Decision.simulation[l,]=reject.decision
    }
  }
  else
  {
    ## MRET test
    for(l in 1:simnum)
    {
      print(paste(l,"/",simnum))
      dataset.whole<-dataset.produce(nc,nt,accrual.type,A,tau,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
      Sw_value=rep(0,numberVisits)
      for(i in 1:numberVisits)
      {
        analysis.time=tau.vector[i]
        dataset.analysis<-dataset.whole[dataset.whole[,4]<=analysis.time,]
        dataset.analysis[dataset.analysis[,8]>analysis.time,6]=0
        dataset.analysis[dataset.analysis[,8]>analysis.time,7]=analysis.time-dataset.analysis[dataset.analysis[,8]>analysis.time,4]
        dataset.analysis[dataset.analysis[,8]>analysis.time,8]=analysis.time
        test.statistic=Maximin.efficiency.robust.test(dataset.input=dataset.analysis,t.low=t1,t.upper=t2)
        Sw_value[i]=test.statistic$Sw_value
        Sw.ustd.simulation[l,i]=test.statistic$Sw.ustd
        SW.Var.H0.simulation[l,i]=test.statistic$SW.Var.H0
        Events.number.total.simulation[l,i]=test.statistic$death.all
        Total.number.simulation[l,i]=test.statistic$num.recru
      }
      test.statistics.simulation[l,]=Sw_value
      reject.decision=!(Sw_value<=upperbound&lowerbound<=Sw_value)## the rejection situation for respective situation,1 represents rejection,0 represents acceptance
      for(i in 1:(numberVisits-1))## once the hypothesis is rejected in one interim ananlysis,the following interim analyses and final analysis won't be conducted,so let them be 0
      {
        if(reject.decision[i]==1)
        {
          for(j in c((i+1):numberVisits))
          {
            reject.decision[j]=0
            Events.number.total.simulation[l,j]=0
            Total.number.simulation[l,j]=0
          }
          break
        } 
      }
      Decision.simulation[l,]=reject.decision
    }
  }
  times.reject.respective<-apply(Decision.simulation,2,sum)
  times.reject<-sum(times.reject.respective)
  Power=times.reject/simnum## total exprical power
  Power.respective=times.reject.respective/simnum## emprical power correspongding to respective stage
  Power.cumulative=rep(0,k)## cumulative emprical power correspongding to respective stage
  Power.cumulative[1]=Power.respective[1]
  for(i in 2:numberVisits)
  {
    Power.cumulative[i]=Power.cumulative[i-1]+Power.respective[i]
  }
  ## calculate the empirical means of the cumulative numbers of the events at every analysis
  ## calculate the empirical means of the cumulative numbers of the recruited patients at every analysis
  Events.mean.cumulative.total.simulation=rep(0,numberVisits)
  Total.mean.cumulative.simulation=rep(0,numberVisits)
  for(i in 1:numberVisits)
  {
    Vector.Events.number.total.simulation=Events.number.total.simulation[,i]
    Vector.Total.number.simulation=Total.number.simulation[,i]
    Events.mean.cumulative.total.simulation[i]=mean(Vector.Events.number.total.simulation[Vector.Events.number.total.simulation!=0])
    Total.mean.cumulative.simulation[i]=mean(Vector.Total.number.simulation[Vector.Total.number.simulation!=0])
  }
  ## calculate the average number of actual deaths
  ## calculate the average number of actual enrolled patients
  Events.total.simulation=rep(0,simnum)
  Total.patients.number.simulation=rep(0,simnum)
  for(l in 1:simnum)
  {
    Events.total.simulation[l]=max(Events.number.total.simulation[l,])
    Total.patients.number.simulation[l]=max(Total.number.simulation[l,])
  }
  Events.average.simulation=mean(Events.total.simulation)
  Average.patients.number.simulation=mean(Total.patients.number.simulation)
  ## calculate the empirical distribution parameters of the test statistics
  mean.test.statistics.std=apply(test.statistics.simulation,2,mean)
  mean.test.statistics.ustd=apply(Sw.ustd.simulation,2,mean)
  Var.test.statistics.matrix.std=matrix(nrow=numberVisits,ncol=numberVisits)
  Var.test.statistics.matrix.H1.ustd=matrix(nrow=numberVisits,ncol=numberVisits)
  Var.test.statistics.matrix.H0.ustd=rep(0,numberVisits)
  for(i in 1:numberVisits)
  {
    Var.test.statistics.matrix.H0.ustd[i]=mean(SW.Var.H0.simulation[,i])
    for(j in 1:numberVisits)
    {
      Var.test.statistics.matrix.std[i,j]=cov(test.statistics.simulation[,i],test.statistics.simulation[,j])
      Var.test.statistics.matrix.H1.ustd[i,j]=cov(Sw.ustd.simulation[,i],Sw.ustd.simulation[,j])
    }
  }
  Var.Var.H0.ustd=var(SW.Var.H0.simulation)
  Parameters<-list(nc=nc,nt=nt,accrual.type=accrual.type,A=A,tau.vector=tau.vector,bound.list=bound.list,loss.control=loss.control,loss.treatment=loss.treatment,k=k,h0=h0,t1=t1,t2=t2,t1.true=t1.true,t2.true=t2.true,theta=theta,a=a,b=b,sides=sides,alpha=alpha,test.type=test.type,simnum=simnum)
  Distribution<-list(mean.test.statistics.std=mean.test.statistics.std,mean.test.statistics.ustd=mean.test.statistics.ustd,Var.test.statistics.matrix.std=Var.test.statistics.matrix.std,Var.test.statistics.matrix.H1.ustd=Var.test.statistics.matrix.H1.ustd,Var.test.statistics.matrix.H0.ustd=Var.test.statistics.matrix.H0.ustd)
  Original.Data=list(SW.Var.H0.simulation=SW.Var.H0.simulation,Sw.ustd.simulation=Sw.ustd.simulation,test.statistics.simulation=test.statistics.simulation)
  return(list(Power=Power,Power.respective=Power.respective,Power.cumulative=Power.cumulative,Events.mean.cumulative.total.simulation=Events.mean.cumulative.total.simulation,Total.mean.cumulative.simulation=Total.mean.cumulative.simulation,Events.average.simulation=Events.average.simulation,Average.patients.number.simulation=Average.patients.number.simulation,Parameters=Parameters,Distribution=Distribution,Original.Data=Original.Data))
}


####################################################
###the functions which are based on the Markov chain
####################################################
#####the following function is used to generate the Markov chain
marcov.survival.states.chain.simple<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t2,theta,a,b)
  ## the function used to simulate the state proportions at different time points
  ## Note: this marcov model is not applicable to the switch treatments under the delayed treatment effect pattern
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the number of the further splitted small intervals for the interval which needs further
  ## interval.period: the duration of the interval which will be further splitted
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## a,b:the shape parameter of the distribution of the delay time
{
  split.num<-1/interval.period*split.interval.num## the splitted number of the unit time interval for an Markov chain
  ##this alogorithm restricts that split.interval.num is exactly divided by interval.period
  interval.num=tau/interval.period## the total number of the intervals which are initially splitted
  unit.interval.period=interval.period/split.interval.num## the period of the smallest splitted interval
  num<-split.num*tau## the total number of the splitted intervals
  loss.control.split<-1-(1-loss.control)^(1/split.num)
  loss.treatment.split<-1-(1-loss.treatment)^(1/split.num)
  Control.matrix<-matrix(0,nrow=(num+1),ncol=4)
  colnames(Control.matrix)<-c("Time","Loss","Event","Surv.cont")
  Control.matrix[1,4]=1
  Treatment.matrix<-matrix(0,nrow=(num+1),ncol=4)
  colnames(Treatment.matrix)<-c("Time","Loss","Event","Surv.Trt")
  Treatment.matrix[1,4]=1
  Trans.matrix.list.control<-list()
  Trans.matrix.list.treatment<-list()
  for(l in 1:interval.num)
  {
    t.low=(l-1)*interval.period
    t.high=l*interval.period
    cum.hazard.control<-integrate(f=hazard.function.control,lower=t.low,upper=t.high,k,h0)$value
    Prop.death.control.split<-1-exp(-cum.hazard.control/split.interval.num)
    cum.hazard.treatment.split<-integrate(f=hazard.function.treatment,lower=t.low,upper=t.high,k,h0,t1,t2,theta,a,b)$value
    Prop.death.treatment.split<-1-exp(-cum.hazard.treatment.split/split.interval.num)
    for(j in 1:split.interval.num)
    {
      if(t.high>tau-A)
      {
        Prob.censor=1/round(num+1-t.low*split.num-j)
      }
      else
      {
        Prob.censor=0
      }
      treatment.resd<-1-loss.treatment.split-Prop.death.treatment.split
      control.resd<-1-loss.control.split-Prop.death.control.split
      trans.matrix.uncensored.treatment<-rbind(c(1,0,0),c(0,1,0),c(loss.treatment.split,Prop.death.treatment.split,treatment.resd))
      trans.matrix.uncensored.control<-rbind(c(1,0,0),c(0,1,0),c(loss.control.split,Prop.death.control.split,control.resd))
      trans.matrix.admin.censor<-rbind(c(1,0,0),c(0,1,0),c(Prob.censor,0,1-Prob.censor))
      trans.matrix.treatment=trans.matrix.uncensored.treatment%*%trans.matrix.admin.censor
      trans.matrix.control=trans.matrix.uncensored.control%*%trans.matrix.admin.censor
      Trans.matrix.treatment=new("markovchain",states=c("Loss","Event","Surv.Trt"),transitionMatrix=trans.matrix.treatment,name="Survival")
      Trans.matrix.control=new("markovchain",states=c("Loss","Event","Surv.cont"),transitionMatrix=trans.matrix.control,name="Survival")
      i=(l-1)*split.interval.num+j
      Trans.matrix.list.control[[i]]=Trans.matrix.control
      Trans.matrix.list.treatment[[i]]=Trans.matrix.treatment
      Time=t.low+j*unit.interval.period
      Control.vector<-Control.matrix[i,2:4]
      Treatment.vector<-Treatment.matrix[i,2:4]
      Control.matrix[i+1,]<-c(Time,Control.vector*Trans.matrix.control)
      Treatment.matrix[i+1,]<-c(Time,Treatment.vector*Trans.matrix.treatment)
    }
  }
  Survival.markovchainList.control<- new("markovchainList", markovchains =Trans.matrix.list.control,name = "Survival Transition list")
  Survival.markovchainList.treatment<- new("markovchainList", markovchains =Trans.matrix.list.treatment,name = "Survival Transition list")
  Survival.markovchainList<-list(Survival.markovchainList.control=Survival.markovchainList.control,Survival.markovchainList.treatment=Survival.markovchainList.treatment)
  Control.matrix=data.frame(Control.matrix)
  Treatment.matrix=data.frame(Treatment.matrix)
  return(list(Control.matrix=Control.matrix,Treatment.matrix=Treatment.matrix))
}

###the following function is used to determine the distribution parameters of a weighted log-rank test statistic
marcov.survival.stat.distribution.chain.simple<-function(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  ## the function used to calculate the distribution parameters for the test statistics
  ## A: the recruitment time(restricted to an integer)
  ## tau:the duration of the whole clinical trial(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time(used for the calculation of Phi and the calculation of MERT)
  ## t1.true,t2.true:the minimum and maximum delay time(used in l(t) for the generation of the Markov chain)
  ## theta: the hazard ratio after complete onset of the treatment effect
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## a,b:the shape parameter of the distribution of the delay time
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  list.survival.states.chain.simple<-marcov.survival.states.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1.true,t2.true,theta,a,b)
  ##obtain the proportions of the states at different time points by the Markov chain
  Pc<-1/(1+Allc)
  Pt<-Allc/(1+Allc)
  Control.matrix<-list.survival.states.chain.simple$Control.matrix
  Control.matrix[,2:4]<-Pc*Control.matrix[,2:4]
  Treatment.matrix<-list.survival.states.chain.simple$Treatment.matrix
  Treatment.matrix[,2:4]<-Pt*Treatment.matrix[,2:4]
  num=tau*split.interval.num/interval.period## the total number of the splited intervals
  ##Calc.matrix: the matrix is used to calculate the distribution parameters of the test statistics; it is made up of the following 16 elements
  Calc.matrix<-data.frame(matrix(0,nrow=num,ncol=21))
  colnames(Calc.matrix)<-c("d.cont","d.trt","d.all","Surv.cont.init","Surv.trt.init","Hazard.cont","Hazard.trt","Hazard.ratio","Surv.cont.end","Surv.trt.end","Cont.tot.exit","Trt.tot.exit","Ratio.Surv","Phi","Survival.time.init","Survival.time.end","Weights","d.porp","gamma","eta","omega")
  #d.cont: the death rate of the control group in a tiny interval
  #d.trt: the death rate of the treatment group in a tiny interval
  #d.all: the death rate in a tiny interval
  #Surv.cont: the proportion of the alive subjects in the control group
  #Surv.trt: the proportion of the alive subjects in the treatment group
  #Hazard.cont: the hazard rate of the control group
  #Hazard.trt: the hazard rate of the treatment group
  #Hazard.ratio: the ratio of the hazard rate of  the control group to that of the treatment group, also denoted by theta
  #Ratio.Surv: the ratio of the proportion of the alive subjects in the control group to that in the treatment group, also denoted by phi
  #Phi:a parameter which is used to calculate the weights of the Maximin Efficiency Robust Test
  #Survival.time: the survival time corresponding to the observation in the matrix
  #Weights:weights function
  #d.porp: the death proportion in the interval
  #gamma:phi*theta/(1+phi*theta)-phi/(1+phi)
  #eta:phi/(1+phi)^2
  #omega:phi*theta/(1+phi*theta)^2
  unit.interval<-interval.period/split.interval.num
  Death.cum=0
  if(test.type==1)
  {
    t1=0
    t2=0
  }
  if(t1==0)
  {
    Phi.t1=0
  }
  if(t2==0)
  {
    Phi.t2=0
  }
  if(Phi.type==1)
  {
    for(i in 1:num)
    {
      #Calc.matrix[i,1]=Control.matrix[i+1,3]-Control.matrix[i,3]
      Calc.matrix$d.cont[i]=Control.matrix$Event[i+1]-Control.matrix$Event[i]
      #Calc.matrix[i,2]=Treatment.matrix[i+1,3]-Treatment.matrix[i,3]
      Calc.matrix$d.trt[i]=Treatment.matrix$Event[i+1]-Treatment.matrix$Event[i]
      #Calc.matrix[i,3]=Calc.matrix[i,1]+Calc.matrix[i,2]
      Calc.matrix$d.all[i]=Calc.matrix$d.cont[i]+Calc.matrix$d.trt[i]
      Death.cum=Death.cum+Calc.matrix$d.all[i]
      #Calc.matrix[i,4]=Control.matrix[i,5]
      Calc.matrix$Surv.cont.init[i]=Control.matrix$Surv.cont[i]
      Calc.matrix$Surv.cont.end[i]=Control.matrix$Surv.cont[i+1]
      #Calc.matrix[i,5]=Treatment.matrix[i,5]
      Calc.matrix$Surv.trt.init[i]=Treatment.matrix$Surv.Trt[i]
      Calc.matrix$Surv.trt.end[i]=Treatment.matrix$Surv.Trt[i+1]
      #Calc.matrix[i,6]=-log(1-Calc.matrix[i,1]/Calc.matrix[i,4])/unit.interval
      Calc.matrix$Hazard.cont[i]=-log(1-Calc.matrix$d.cont[i]/Calc.matrix$Surv.cont.init[i])/unit.interval
      #Calc.matrix[i,7]=-log(1-Calc.matrix[i,2]/Calc.matrix[i,5])/unit.interval
      Calc.matrix$Hazard.trt[i]=-log(1-Calc.matrix$d.trt[i]/Calc.matrix$Surv.trt.init[i])/unit.interval
      #Calc.matrix[i,8]=Calc.matrix[i,6]/Calc.matrix[i,7]
      Calc.matrix$Hazard.ratio[i]=Calc.matrix$Hazard.cont[i]/Calc.matrix$Hazard.trt[i]
      ##Calc.matrix[i,9]=Calc.matrix[i,4]/Calc.matrix[i,5]
      Calc.matrix$Ratio.Surv[i]=Calc.matrix$Surv.cont.init[i]/Calc.matrix$Surv.trt.init[i]
      Calc.matrix$Trt.tot.exit[i]=Treatment.matrix$Loss[i+1]
      Calc.matrix$Cont.tot.exit[i]=Control.matrix$Loss[i+1]
      if(i==1)
      {
        integrand.low<-hazard.function.control(0,k,h0)*Pc*Pt
        integrand.high.numerator<-Calc.matrix$Hazard.cont[i]*Calc.matrix$Surv.cont.end[i]*(Pt-Calc.matrix$Trt.tot.exit[i])
        integrand.high.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i])+(Pc-Calc.matrix$Cont.tot.exit[i])
        integrand.high<-integrand.high.numerator/integrand.high.denominator
        Calc.matrix$Phi[i]=(integrand.low+integrand.high)/2*unit.interval
      }
      else
      {
        integrand.low.numerator<-Calc.matrix$Hazard.cont[i-1]*Calc.matrix$Surv.cont.end[i-1]*(Pt-Calc.matrix$Trt.tot.exit[i-1])
        integrand.low.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i-1])+(Pc-Calc.matrix$Cont.tot.exit[i-1])
        integrand.low<-integrand.low.numerator/integrand.low.denominator
        integrand.high.numerator<-Calc.matrix$Hazard.cont[i]*Calc.matrix$Surv.cont.end[i]*(Pt-Calc.matrix$Trt.tot.exit[i])
        integrand.high.denominator<-(Pt-Calc.matrix$Trt.tot.exit[i])+(Pc-Calc.matrix$Cont.tot.exit[i])
        integrand.high<-integrand.high.numerator/integrand.high.denominator
        Calc.matrix$Phi[i]=Calc.matrix$Phi[i-1]+(integrand.low+integrand.high)/2*unit.interval
      }
      Survival.time.init=Control.matrix[i,1]
      Survival.time.end=Control.matrix[i+1,1]
      Calc.matrix$Survival.time.init[i]=Survival.time.init
      Calc.matrix$Survival.time.end[i]=Survival.time.end
      ##remove the influences of machine errors, which results from the summation of the survival time
      if(round(Survival.time.end,8)==t1)
      {
        Phi.t1=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==t2)
      {
        Phi.t2=Calc.matrix$Phi[i]
      }
      #if(round(Survival.time.end,8)==tau)
      #{
      #  Phi.tau=Calc.matrix$Phi[i]
      #}
      if(i==num-1)
      {
        Phi.tau=Calc.matrix$Phi[i]
      }
    }
  }
  else
  {
    for(i in 1:num)
    {
      #Calc.matrix[i,1]=Control.matrix[i+1,3]-Control.matrix[i,3]
      Calc.matrix$d.cont[i]=Control.matrix$Event[i+1]-Control.matrix$Event[i]
      #Calc.matrix[i,2]=Treatment.matrix[i+1,3]-Treatment.matrix[i,3]
      Calc.matrix$d.trt[i]=Treatment.matrix$Event[i+1]-Treatment.matrix$Event[i]
      #Calc.matrix[i,3]=Calc.matrix[i,1]+Calc.matrix[i,2]
      Calc.matrix$d.all[i]=Calc.matrix$d.cont[i]+Calc.matrix$d.trt[i]
      Death.cum=Death.cum+Calc.matrix$d.all[i]
      #Calc.matrix[i,4]=Control.matrix[i,5]+Control.matrix[i,6]
      Calc.matrix$Surv.cont.init[i]=Control.matrix$Surv.cont[i]
      Calc.matrix$Surv.cont.end[i]=Control.matrix$Surv.cont[i+1]
      #Calc.matrix[i,5]=Treatment.matrix[i,5]+Treatment.matrix[i,6]
      Calc.matrix$Surv.trt.init[i]=Treatment.matrix$Surv.Trt[i]
      Calc.matrix$Surv.trt.end[i]=Treatment.matrix$Surv.Trt[i+1]
      #Calc.matrix[i,6]=-log(1-Calc.matrix[i,1]/Calc.matrix[i,4])/unit.interval
      Calc.matrix$Hazard.cont[i]=-log(1-Calc.matrix$d.cont[i]/Calc.matrix$Surv.cont.init[i])/unit.interval
      #Calc.matrix[i,7]=-log(1-Calc.matrix[i,2]/Calc.matrix[i,5])/unit.interval
      Calc.matrix$Hazard.trt[i]=-log(1-Calc.matrix$d.trt[i]/Calc.matrix$Surv.trt.init[i])/unit.interval
      #Calc.matrix[i,8]=Calc.matrix[i,6]/Calc.matrix[i,7]
      Calc.matrix$Hazard.ratio[i]=Calc.matrix$Hazard.cont[i]/Calc.matrix$Hazard.trt[i]
      ##Calc.matrix[i,9]=Calc.matrix[i,4]/Calc.matrix[i,5]
      Calc.matrix$Ratio.Surv[i]=Calc.matrix$Surv.cont.init[i]/Calc.matrix$Surv.trt.init[i]
      Calc.matrix$Trt.tot.exit[i]=Treatment.matrix$Loss[i+1]
      Calc.matrix$Cont.tot.exit[i]=Control.matrix$Loss[i+1]
      if(i==1)
      {
        #Calc.matrix[i,10]=Calc.matrix[i,3]*Calc.matrix[i,4]*Calc.matrix[i,5]/(Calc.matrix[i,4]+Calc.matrix[i,5])^2
        Calc.matrix$Phi[i]=Calc.matrix$d.all[i]*Calc.matrix$Surv.cont.init[i]*Calc.matrix$Surv.trt.init[i]/(Calc.matrix$Surv.cont.init[i]+Calc.matrix$Surv.trt.init[i])^2
      }
      else
      {
        #Calc.matrix[i,10]=Calc.matrix[i-1,10]+Calc.matrix[i,3]*Calc.matrix[i,4]*Calc.matrix[i,5]/(Calc.matrix[i,4]+Calc.matrix[i,5])^2
        Calc.matrix$Phi[i]=Calc.matrix$Phi[i-1]+Calc.matrix$d.all[i]*Calc.matrix$Surv.cont.init[i]*Calc.matrix$Surv.trt.init[i]/(Calc.matrix$Surv.cont.init[i]+Calc.matrix$Surv.trt.init[i])^2
      }
      Survival.time.init=Control.matrix[i,1]
      Survival.time.end=Control.matrix[i+1,1]
      Calc.matrix$Survival.time.init[i]=Survival.time.init
      Calc.matrix$Survival.time.end[i]=Survival.time.end
      ##remove the influences of machine errors, which results from the summation of the survival time
      if(round(Survival.time.end,8)==t1)
      {
        Phi.t1=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==t2)
      {
        Phi.t2=Calc.matrix$Phi[i]
      }
      if(round(Survival.time.end,8)==tau)
      {
        Phi.tau=Calc.matrix$Phi[i]
      }
    }
  }
  for(i in 1:num)
  {
    Survival.time.end=Calc.matrix$Survival.time.end[i]
    ##Caution:<=t1 or <t1,please determine by debug
    if(round(Survival.time.end,8)<t1)
    {
      Calc.matrix$Weights[i]=0
    }
    else if(round(Survival.time.end,8)>=t1 & round(Survival.time.end,8)<=t2)
    {
      ##Caution:>=t1 or >t1,please determine by debug
      Calc.matrix$Weights[i]=((Phi.tau-Calc.matrix$Phi[i])/(Phi.tau-Phi.t1))^(-0.5)
      #Note:tau is commonly assumed to be greated than t2; therefore, don't worry that the above value is Inf
    }
    else
    {
      Calc.matrix$Weights[i]=2*((Phi.tau-Phi.t2)/(Phi.tau-Phi.t1))^(-0.5)
    }
    Calc.matrix$d.porp[i]=Calc.matrix$d.all[i]/Death.cum
  }
  ##The parameters for the number of deaths
  Mean.unstd.d.unit=0
  Var.H0.unstd.d.unit=0
  Var.H1.unstd.d.unit=0
  if(test.type==1)
  {
    Calc.matrix$Weights=Calc.matrix$Weights/2
    ##log-rank test
    for(i in 1:num)
    {
      phi=Calc.matrix$Ratio.Surv[i]
      theta=Calc.matrix$Hazard.ratio[i]
      gamma=phi*theta/(1+phi*theta)-phi/(1+phi)
      #Calc.matrix[i,14]=gamma
      Calc.matrix$gamma[i]=gamma
      eta=phi/(1+phi)^2
      #Calc.matrix[i,15]=eta
      Calc.matrix$eta[i]=eta
      omega=phi*theta/(1+phi*theta)^2
      #Calc.matrix[i,16]=omega
      Calc.matrix$omega[i]=omega
      #rho=Calc.matrix[i,13]
      rho=Calc.matrix$d.porp[i]
      Mean.unstd.d.unit=Mean.unstd.d.unit+gamma*rho
      Var.H0.unstd.d.unit=Var.H0.unstd.d.unit+eta*rho
      Var.H1.unstd.d.unit=Var.H1.unstd.d.unit+omega*rho
    }
  }
  else
  {
    ## Maximin efficiency robust test
    for(i in 1:num)
    {
      weight=Calc.matrix$Weights[i]
      phi=Calc.matrix$Ratio.Surv[i]
      theta=Calc.matrix$Hazard.ratio[i]
      gamma=phi*theta/(1+phi*theta)-phi/(1+phi)
      #Calc.matrix[i,14]=gamma
      Calc.matrix$gamma[i]=gamma
      eta=phi/(1+phi)^2
      #Calc.matrix[i,15]=eta
      Calc.matrix$eta[i]=eta
      omega=phi*theta/(1+phi*theta)^2
      #Calc.matrix[i,16]=omega
      Calc.matrix$omega[i]=omega
      #rho=Calc.matrix[i,13]
      rho=Calc.matrix$d.porp[i]
      Mean.unstd.d.unit=Mean.unstd.d.unit+gamma*rho*weight
      Var.H0.unstd.d.unit=Var.H0.unstd.d.unit+rho*eta*weight^2
      Var.H1.unstd.d.unit=Var.H1.unstd.d.unit+rho*omega*weight^2
    }
  }
  Mean.unstd.d.unit=as.numeric(Mean.unstd.d.unit)
  Var.H0.unstd.d.unit=as.numeric(Var.H0.unstd.d.unit)
  Var.H1.unstd.d.unit=as.numeric(Var.H1.unstd.d.unit)
  return(list(Mean.unstd.d.unit=Mean.unstd.d.unit,Var.H0.unstd.d.unit=Var.H0.unstd.d.unit,Var.H1.unstd.d.unit=Var.H1.unstd.d.unit,Calc.matrix=Calc.matrix,Death.cum=as.numeric(Death.cum),list.survival.states.chain.simple=list.survival.states.chain.simple))
}

group.sequential.boundary.by.markov.chain<-function(tau.vector,A,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t2,Allc,test.type,Phi.type,alpha,alpha.allc,sides)
  ## This function is used to predict the information time for every analysis
  ## A: the recruitment time(restricted to an integer)
  ## tau.vector:the calendar time of every analysis(restricted to an integer)
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time (used for the calculation of Phi)
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
  ## alpha: the significance level
  ## alpha.allc: the allocation fraction of the significance level
  ## sides:1:the right-sided test; 2:the left-sided test
{
  num.analysis=length(tau.vector)
  ## obtain the markov chain for every analysis
  list.marcov.survival.stat.distribution.chain.simple<-list()
  for(i in 1:num.analysis)
  {
    tau=tau.vector[i]
    list.marcov.survival.stat.distribution.chain.simple[[i]]<-marcov.survival.stat.distribution.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true=0,t2,t2.true=0,theta=1,a=0,b=0,Allc,test.type,Phi.type)
  }
  ## obtain the multivariate distribution of the MERT statistics at all the analyses under the null hypothesis
  mean.vector=rep(0,num.analysis)
  covariance.matrix<-matrix(0,nrow=num.analysis,ncol=num.analysis)
  for(i in 1:num.analysis)
  {
    covariance.matrix[i,i]=1
  }
  for(i in 1:(num.analysis-1))
  {
    for(j in (i+1):num.analysis)
    {
      list.pred<-list.marcov.survival.stat.distribution.chain.simple[[i]]
      list.latt<-list.marcov.survival.stat.distribution.chain.simple[[j]]
      Calc.matrix.pred<-list.pred$Calc.matrix
      Calc.matrix.latt<-list.latt$Calc.matrix
      Death.cum.pred<-list.pred$Death.cum
      Var.H0.unstd.d.unit.pred<-list.pred$Var.H0.unstd.d.unit
      Death.cum.latt<-list.latt$Death.cum
      Var.H0.unstd.d.unit.latt<-list.latt$Var.H0.unstd.d.unit
      time.pred<-Calc.matrix.pred$Survival.time.init
      num.pred<-length(time.pred)
      #weights.latter<-Calc.matrix.latt$Weights[Calc.matrix.latt$Survival.time.init==time.pred]
      weights.latter<-Calc.matrix.latt$Weights[1:num.pred]
      Calc.matrix.pred<-cbind(Calc.matrix.pred,weights.latter)
      Cov.H0.d.unit=0
      for(l in 1:num.pred)
      {
        weight.pred=Calc.matrix.pred$Weights[l]
        weight.latter=Calc.matrix.pred$weights.latter[l]
        phi=Calc.matrix.pred$Ratio.Surv[l]
        eta=phi/(1+phi)^2
        rho=Calc.matrix.pred$d.porp[l]
        Cov.H0.d.unit=Cov.H0.d.unit+rho*eta*weight.pred*weight.latter
      }
      Cor.H0=(Cov.H0.d.unit*Death.cum.pred)/sqrt(Death.cum.pred*Var.H0.unstd.d.unit.pred)/sqrt(Death.cum.latt*Var.H0.unstd.d.unit.latt)
      covariance.matrix[j,i]=Cor.H0
      covariance.matrix[i,j]=Cor.H0
    }
  }
  alpha.vector<-alpha*alpha.allc
  alpha.cum<-cumsum(alpha.vector)
  upperbound<-rep(NA,num.analysis)
  lowerbound<-rep(NA,num.analysis)
  library(mvtnorm)
  if(sides==1)
  {
    lowerbound<-rep(-Inf,num.analysis)
    upperbound[1]=qnorm(1-alpha.vector[1])
    for(i in 2:num.analysis)
    {
      mean.used=mean.vector[1:i]
      covariance.matrix.used=covariance.matrix[1:i,1:i]
      upper.search=qnorm(1-alpha)
      init=0
      repeat
      {
        upper=c(upperbound[1:(i-1)],upper.search)
        alpha.cum.itr=1-pmvnorm(upper=upper,mean=mean.used,sigma=covariance.matrix.used)
        if(init==1)
        {
          if(alpha.cum.itr.pre<alpha.cum[i]&alpha.cum.itr>=alpha.cum[i])
          {
            break
          }
          if(alpha.cum.itr.pre>=alpha.cum[i]&alpha.cum.itr<alpha.cum[i])
          {
            break
          }
        }
        if(alpha.cum.itr>alpha.cum[i])
        {
          alpha.cum.itr.pre=alpha.cum.itr
          upper.search.pre=upper.search
          upper.search=upper.search*1.1
          init=1
        }
        else
        {
          alpha.cum.itr.pre=alpha.cum.itr
          upper.search.pre=upper.search
          upper.search=upper.search*0.9
          init=1
        }
      }
      repeat
      {
        upper.search.mid=(upper.search+upper.search.pre)/2
        upper=c(upperbound[1:(i-1)],upper.search.mid)
        alpha.cum.itr=1-pmvnorm(upper=upper,mean=mean.used,sigma=covariance.matrix.used)
        if(abs(alpha.cum.itr-alpha.cum[i])<=0.000000000001)
        {
          upperbound[i]=upper.search.mid
          break
        }
        if(alpha.cum.itr>alpha.cum[i])
        {
          upper.rep=max(upper.search,upper.search.pre)
          upper.search.pre=upper.rep
          upper.search=upper.search.mid
        }
        else
        {
          upper.rep=min(upper.search,upper.search.pre)
          upper.search.pre=upper.rep
          upper.search=upper.search.mid
        }
      }
    }
  }
  else
  {
    lowerbound[1]=-qnorm(1-alpha.vector[1]/2)
    upperbound[1]=qnorm(1-alpha.vector[1]/2)
    for(i in 2:num.analysis)
    {
      mean.used=mean.vector[1:i]
      covariance.matrix.used=covariance.matrix[1:i,1:i]
      upper.search=qnorm(1-alpha/2)
      lower.search=-qnorm(1-alpha/2)
      init=0
      repeat
      {
        upper=c(upperbound[1:(i-1)],upper.search)
        lower=c(lowerbound[1:(i-1)],lower.search)
        alpha.cum.itr=1-pmvnorm(lower=lower,upper=upper,mean=mean.used,sigma=covariance.matrix.used)
        if(init==1)
        {
          if(alpha.cum.itr.pre<alpha.cum[i]&alpha.cum.itr>=alpha.cum[i])
          {
            break
          }
          if(alpha.cum.itr.pre>=alpha.cum[i]&alpha.cum.itr<alpha.cum[i])
          {
            break
          }
        }
        if(alpha.cum.itr>alpha.cum[i])
        {
          alpha.cum.itr.pre=alpha.cum.itr
          upper.search.pre=upper.search
          lower.search.pre=lower.search
          upper.search=upper.search*1.1
          lower.search=lower.search*1.1
          init=1
        }
        else
        {
          alpha.cum.itr.pre=alpha.cum.itr
          upper.search.pre=upper.search
          lower.search.pre=lower.search
          upper.search=upper.search*0.9
          lower.search=lower.search*0.9
          init=1
        }
      }
      repeat
      {
        upper.search.mid=(upper.search+upper.search.pre)/2
        lower.search.mid=(lower.search+lower.search.pre)/2
        upper=c(upperbound[1:(i-1)],upper.search.mid)
        lower=c(lowerbound[1:(i-1)],lower.search)
        alpha.cum.itr=1-pmvnorm(lower=lower,upper=upper,mean=mean.used,sigma=covariance.matrix.used)
        if(abs(alpha.cum.itr-alpha.cum[i])<=10^(-20))
        {
          lowerbound[i]=lower.search.mid
          upperbound[i]=upper.search.mid
          break
        }
        if(alpha.cum.itr>alpha.cum[i])
        {
          upper.rep=max(upper.search,upper.search.pre)
          lower.rep=min(lower.search,lower.search.pre)
          upper.search.pre=upper.rep
          lower.search.pre=lower.rep
          upper.search=upper.search.mid
          lower.search=lower.search.mid
        }
        else
        {
          upper.rep=min(upper.search,upper.search.pre)
          lower.rep=max(lower.search,lower.search.pre)
          upper.search.pre=upper.rep
          lower.search.pre=lower.rep
          upper.search=upper.search.mid
          lower.search=lower.search.mid
        }
      }
    }
  }
  bound.list<-list(lowerbound=lowerbound,upperbound=upperbound)
  return(list(bound.list=bound.list,covariance.matrix=covariance.matrix,list.marcov.survival.stat.distribution.chain.simple=list.marcov.survival.stat.distribution.chain.simple))
}

power.calculation.group.sequential.design.by.marcov.chain<-function(n,A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  ## This function is used to predict the information time for every analysis
  ## This function assumes that interim analysis is projected after the termination of the recruitment
  ## n: the maximum sample sizes
  ## A: the recruitment time(restricted to an integer)
  ## tau.vector:the calendar time of every analysis(restricted to an integer)
  ## boundary.list:made up of the lowerbound and upperbound
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2: the minimum and maximum delay time(which is used in the MERT)
  ## t1.true,t2.true:the minimum and maximum delay time(which is used in the delay function l(t) and the generation of the Markov chain)  ## theta: the hazard ratio after complete onset of the treatment effect
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  num.analysis=length(tau.vector)
  ##obtain the markov chain for every analysis
  list.marcov.survival.stat.distribution.chain.simple<-list()
  for(i in 1:num.analysis)
  {
    tau=tau.vector[i]
    list.marcov.survival.stat.distribution.chain.simple[[i]]<-marcov.survival.stat.distribution.chain.simple(A,tau,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  }
  num.power<-length(n)
  if(num.power==1)
  {
    ## only one maximum sample size is obtained
    ## obtain the multivariate distribution of the MERT statistics at all the analyses under the alternative hypothesis
    mean.vector=rep(NA,num.analysis)
    covariance.matrix<-matrix(NA,nrow=num.analysis,ncol=num.analysis)
    Death.num.cum=rep(NA,num.analysis)
    Death.prop.cum=rep(NA,num.analysis)
    for(i in 1:num.analysis)
    {
      Mean.unstd.d.unit<-list.marcov.survival.stat.distribution.chain.simple[[i]]$Mean.unstd.d.unit
      Var.H0.unstd.d.unit<-list.marcov.survival.stat.distribution.chain.simple[[i]]$Var.H0.unstd.d.unit
      Var.H1.unstd.d.unit<-list.marcov.survival.stat.distribution.chain.simple[[i]]$Var.H1.unstd.d.unit
      Death.prop.cum[i]<-list.marcov.survival.stat.distribution.chain.simple[[i]]$Death.cum
      Death.num.cum[i]=n*min(1,tau.vector[i]/A)*Death.prop.cum[i]
      ## Note that the subjects may not be wholly recruited at the interim analysis
      covariance.matrix[i,i]=Var.H1.unstd.d.unit/Var.H0.unstd.d.unit
      mean.vector[i]=sqrt(Death.num.cum[i])*Mean.unstd.d.unit/sqrt(Var.H0.unstd.d.unit)
    }
    for(i in 1:(num.analysis-1))
    {
      for(j in (i+1):num.analysis)
      {
        list.pred<-list.marcov.survival.stat.distribution.chain.simple[[i]]
        list.latt<-list.marcov.survival.stat.distribution.chain.simple[[j]]
        Calc.matrix.pred<-list.pred$Calc.matrix
        Calc.matrix.latt<-list.latt$Calc.matrix
        Var.H0.unstd.d.unit.pred<-list.pred$Var.H0.unstd.d.unit
        Var.H0.unstd.d.unit.latt<-list.latt$Var.H0.unstd.d.unit
        time.pred<-Calc.matrix.pred$Survival.time.init
        num.pred<-length(time.pred)
        #weights.latter<-Calc.matrix.latt$Weights[Calc.matrix.latt$Survival.time.init==time.pred]
        weights.latter<-Calc.matrix.latt$Weights[1:num.pred]
        Calc.matrix.pred<-cbind(Calc.matrix.pred,weights.latter)
        Cov.H0.d.unit=0
        for(l in 1:num.pred)
        {
          weight.pred=Calc.matrix.pred$Weights[l]
          weight.latter=Calc.matrix.pred$weights.latter[l]
          phi=Calc.matrix.pred$Ratio.Surv[l]
          theta=Calc.matrix.pred$Hazard.ratio[l]
          omega=phi*theta/(1+phi*theta)^2
          rho=Calc.matrix.pred$d.porp[l]
          Cov.H0.d.unit=Cov.H0.d.unit+rho*omega*weight.pred*weight.latter
        }
        Cor.H0=(Cov.H0.d.unit*Death.num.cum[i])/sqrt(Death.num.cum[i]*Var.H0.unstd.d.unit.pred)/sqrt(Death.num.cum[j]*Var.H0.unstd.d.unit.latt)
        covariance.matrix[j,i]=Cor.H0
        covariance.matrix[i,j]=Cor.H0
      }
    }
    lowerbound=boundary.list$lowerbound
    upperbound=boundary.list$upperbound
    Power=1-pmvnorm(lower=lowerbound,upper=upperbound,mean=mean.vector,sigma=covariance.matrix)
    Distribution.list=list(mean.vector=mean.vector,covariance.matrix=covariance.matrix)
  }
  else
  {
    ## a lot of maximum sample sizes are obtained
    mean.matrix=matrix(nrow=num.power,ncol=num.analysis)
    list.covariance.matrix<-list()
    for(l in 1:num.power)
    {
      list.covariance.matrix[[l]]=matrix(nrow=num.analysis,ncol=num.analysis)
    }
    Death.num.cum.matrix=matrix(nrow=num.power,ncol=num.analysis)
    Death.prop.cum.matrix=matrix(nrow=num.power,ncol=num.analysis)
    lowerbound=boundary.list$lowerbound
    upperbound=boundary.list$upperbound
    Power=rep(NA,num.power)
    for(l in 1:num.power)
    {
      for(i in 1:num.analysis)
      {
        Mean.unstd.d.unit=list.marcov.survival.stat.distribution.chain.simple[[i]]$Mean.unstd.d.unit
        Var.H0.unstd.d.unit=list.marcov.survival.stat.distribution.chain.simple[[i]]$Var.H0.unstd.d.unit
        Var.H1.unstd.d.unit=list.marcov.survival.stat.distribution.chain.simple[[i]]$Var.H1.unstd.d.unit
        Death.prop.cum.matrix[l,i]=list.marcov.survival.stat.distribution.chain.simple[[i]]$Death.cum
        Death.num.cum.matrix[l,i]=n[l]*min(1,tau.vector[i]/A)*Death.prop.cum.matrix[l,i]
        list.covariance.matrix[[l]][i,i]=Var.H1.unstd.d.unit/Var.H0.unstd.d.unit
        mean.matrix[l,i]=sqrt(Death.num.cum.matrix[l,i])*Mean.unstd.d.unit/sqrt(Var.H0.unstd.d.unit)
      }
      for(i in 1:(num.analysis-1))
      {
        for(j in (i+1):num.analysis)
        {
          list.pred<-list.marcov.survival.stat.distribution.chain.simple[[i]]
          list.latt<-list.marcov.survival.stat.distribution.chain.simple[[j]]
          Calc.matrix.pred<-list.pred$Calc.matrix
          Calc.matrix.latt<-list.latt$Calc.matrix
          Var.H0.unstd.d.unit.pred<-list.pred$Var.H0.unstd.d.unit
          Var.H0.unstd.d.unit.latt<-list.latt$Var.H0.unstd.d.unit
          time.pred<-Calc.matrix.pred$Survival.time.init
          num.pred<-length(time.pred)
          #weights.latter<-Calc.matrix.latt$Weights[Calc.matrix.latt$Survival.time.init==time.pred]
          weights.latter<-Calc.matrix.latt$Weights[1:num.pred]
          Calc.matrix.pred<-cbind(Calc.matrix.pred,weights.latter)
          Cov.H0.d.unit=0
          for(h in 1:num.pred)
          {
            weight.pred=Calc.matrix.pred$Weights[h]
            weight.latter=Calc.matrix.pred$weights.latter[h]
            phi=Calc.matrix.pred$Ratio.Surv[h]
            theta=Calc.matrix.pred$Hazard.ratio[h]
            omega=phi*theta/(1+phi*theta)^2
            rho=Calc.matrix.pred$d.porp[h]
            Cov.H0.d.unit=Cov.H0.d.unit+rho*omega*weight.pred*weight.latter
          }
          Cor.H0=(Cov.H0.d.unit*Death.num.cum.matrix[l,i])/sqrt(Death.num.cum.matrix[l,i]*Var.H0.unstd.d.unit.pred)/sqrt(Death.num.cum.matrix[l,j]*Var.H0.unstd.d.unit.latt)
          list.covariance.matrix[[l]][j,i]=Cor.H0
          list.covariance.matrix[[l]][i,j]=Cor.H0
        }
      }
      Power[l]=1-pmvnorm(lower=lowerbound,upper=upperbound,mean=mean.matrix[l,],sigma=list.covariance.matrix[[l]]) 
      Distribution.list=list(mean.matrix=mean.matrix,list.covariance.matrix=list.covariance.matrix)
    }
  }
  Parameters.list<-list(n=n,A=A,tau.vector=tau.vector,boundary.list=boundary.list,split.interval.num=split.interval.num,interval.period=interval.period,loss.control=loss.control,loss.treatment=loss.treatment,k=k,h0=h0,t1=t1,t1.true=t1.true,t2=t2,t2.true=t2.true,theta=theta,a=a,b=b,Allc=Allc,test.type=test.type,Phi.type=Phi.type)
  return(list(Power=Power,Distribution.list=Distribution.list,Parameters.list=Parameters.list))
}

sample.size.reg<-function(Power.sample,Power)
{
  Power.sample<-data.frame(Power.sample)
  names(Power.sample)<-c("sample.size","power","qnorm.power","sqrtn")
  fit<-lm(formula=sqrtn~qnorm.power,data=Power.sample)
  print(summary(fit))
  Power.point=data.frame(qnorm.power=qnorm(Power))
  Predict.sample<-predict(fit,Power.point,interval="prediction",level=0.95)
  sample.size<-round(Predict.sample[1]^2)
  return(sample.size)
}

maximum.sample.size.search.by.marcov.chain<-function(power,A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  ## This function is used to determine the required maximum sample size by the iterative regressions between the square root of the maximum sample sizes and the normal quantile of the maximum sample size
  ## This function assumes that interim analysis is projected after the termination of the recruitment
  ## power: the normal power
  ## A: the recruitment time(restricted to an integer)
  ## tau.vector:the calendar time of every analysis(restricted to an integer)
  ## boundary.list:made up of the lowerbound and upperbound
  ## split.interval.num: the splitted number of the smallest interval for an Markov chain
  ## interval.period: the duration of the smallest interval
  ## loss.control: the rate of loss to follow-up in an unit period for the control group
  ## loss.treatment: the rate of loss to follow-up in an unit period for the treatment group
  ## k,h0: the survival function of the control group is assumed to be exp[-(h0*t)^k]
  ## t1,t2:the minimum and maximum delay time (used for the calculation of Phi)
  ## Allc: the allocation ratio of the sample sizes of the treatment and control groups
  ## test.type:1:log-rank test;2:Maximin efficiency robustness test
  ## Phi.type:1:Phi obtained by numeric integration;2:estimated Phi
{
  Power.sample<-matrix(nrow=4,ncol=4)
  colnames(Power.sample)<-c("sample.size","power","qnorm.power","sqrtn")
  n.init<-c(50,200,350,500)
  for(i in 1:4)
  {
    Power.sample[i,1]=n.init[i]
    Power.sample[i,2]=power.calculation.group.sequential.design.by.marcov.chain(n.init[i],A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)$Power
    Power.sample[i,3]=qnorm(Power.sample[i,2])
    Power.sample[i,4]=sqrt(n.init[i])
  }
  itr.num=4
  repeat
  {
    n.tir<-sample.size.reg(Power.sample,power)
    Power.itr<-power.calculation.group.sequential.design.by.marcov.chain(n.tir,A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)$Power
    itr.add=c(n.tir,Power.itr,qnorm(Power.itr),sqrt(n.tir))
    Power.sample<-rbind(Power.sample,itr.add)
    itr.num=itr.num+1
    if(abs(Power.sample[itr.num,1]-Power.sample[itr.num-1,1])<=1)
    {
      break
    }
  }
  if(Power.itr<power)
  {
    repeat
    {
      n.tir=n.tir+1
      Power.itr<-power.calculation.group.sequential.design.by.marcov.chain(n.tir,A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
      if(Power.itr$Power>=power)
      {
        break
      }
    }
    maximum.sample.size=n.tir
    analytic.power=Power.itr
  }
  else
  {
    maximum.sample.size=Power.sample[itr.num,1]
    analytic.power=power.calculation.group.sequential.design.by.marcov.chain(maximum.sample.size,A,tau.vector,boundary.list,split.interval.num,interval.period,loss.control,loss.treatment,k,h0,t1,t1.true,t2,t2.true,theta,a,b,Allc,test.type,Phi.type)
  }
  nc=round(maximum.sample.size/(Allc+1))
  nt=maximum.sample.size-nc
  return(list(maximum.sample.size=maximum.sample.size,analytic.power=analytic.power,nc=nc,nt=nt))
}
