# Estimate F from tag-recapture data

#NOTES: this script implements McAuley et al's 2007 method

library(fishmethods)
rm(list=ls(all=TRUE))
setwd("C:/Matias/Analyses/Conventional tagging/F_estimation")



#---DATA SECTION-------

#1. tagging data
source("C:/Matias/Analyses/SOURCE_SCRIPTS/Source_conventional_data.R")





#---PARAMETERS SECTION-------

#Whiskery shark
  #Max Age
A.max=20

  #parturition
Part="1 September"

  #Growth parameters
    #female (Simfendorfer et al 2000)
Linf.w=120.7   #ageing
K.w=0.369
to.w=-0.544
sd.w=10  #assumed SD of vonB function (no SD provided by Simpfendorfer et al 2000)







#---PROCEDURE SECTION-------

#1. Create useful vars
SPECIES=c("TK","BW","WH","GM")
Species.names=c("sandbar shark","dusky shark","whiskery shark","gummy shark")
names(Species.names)=SPECIES

#Time at liberty
Tagging$DaysAtLarge=as.numeric(round((Tagging$DATE_CAPTR-Tagging$"RELEASE DATE")/(60*60*24),0))



#Predict age from observed length
fn.age.len=function(A.max,Linf,K,to,sd,obs.length)
{
  age=0:A.max
  
  #von B function
  mean.length=Linf*(1-exp(-K*(age-to)))
  
  sd=rep(sd,length(age))
  
  #prob of age given length
  prob=function(data)
  {
    dnorm(data,mean.length,sd)
  }
  pred.len=sapply(obs.length,prob)
  
  id=rep(NA,ncol(pred.len))
  for (i in 1:ncol(pred.len))  id[i]=which(pred.len[,i]==max(pred.len[,i]))
  
  pred.age=age[id]
  
  plot(age,mean.length,type='l',lwd=2)
  points(pred.age,obs.length,pch=19,col=2)
  
  return(pred.age=pred.age)
}

Tagging.WH=subset(Tagging,Species=="WH" & !(is.na(Rel_FL)))

obs.len=Tagging.WH$Rel_FL

Tagging.WH$Pred.Age=fn.age.len(A.max,Linf.w,K.w,to.w,sd.w,obs.len)

#drop those in condition 3
table(Tagging.WH$CONDITION)
Tagging.WH=subset(Tagging.WH,!CONDITION==3)

#drop those at liberty <90 days
#Tagging.WH=subset(Tagging.WH,!CONDITION==3)


#Reporting/unreporting rates                
Yrs=unique(Tagging.WH$Year)
Yrs=Yrs[!(is.na(Yrs))]
AREA=unique(Tagging.WH$Areas)

nonRep.Rate=expand.grid(Areas=AREA,Year=Yrs)  #DUMMY, GET FROM RORY!!!
nonRep.Rate$nr.rate=0.2


#Estimate Z

#1. "fishmethods" package
#Hoenig
# irm_cr
# irm_h
# tag_model_avg





#2. McAuley et al 2007 implementation
Estim.F.at.age=function(DATA,NR.Rate) #INCOMPLETE!!!
{
  #1. Correct reported tagged catches
  DATA$Reported=1
  Tab=aggregate(Reported~Pred.Age+Areas+Month+Year,DATA,sum)
  Tab=merge(Tab,NR.Rate,by=c("Areas","Year"),all.x=T)
    
  Tab$Corrected=with(Tab,Reported/(1-nr.rate))
  
  #2. Numbers tagged at start of month
  
  
}
#Estim.F.at.age(Tagging.WH,nonRep.Rate)



#--------------------Second alternative
#   Script for calculating exploitation rates from tag-recaptures as done by Simpfendorfer et al (1999)
#   and McAuley et al 2005 and 2007



# Simpfendorfer et al (1999)
#DATA
#number of regions
n.reg=3
#number of months
n.mon=12
#number of tagged cohorts
n.coho=2

#number of tags reported from cohort a in region i during month t
#dummy
n.rep=round(runif(n.reg*n.mon*n.coho,0,20))
R=array(n.rep,c(n.mon,n.reg,n.coho))
dimnames(R)[[2]]=c('Region1','Region2','Region3')
dimnames(R)[[3]]=c('cohort1','cohort2')


#Actual data            #Simpfendorfer et al (1999)
R[,,1][,1]=c(3,11,13,1,12,8,1,6,1,5,5,3)
R[,,1][,2]=rep(0,12)
R[,,1][,3]=c(0,0,0,4,4,8,4,2,3,0,2,0)
R[,,2][,1]=c(5,3,3,4,2,0,0,0,1,2,2,0)
R[,,2][,2]=rep(0,12)
R[,,2][,3]=c(2,0,2,2,0,0,0,0,0,1,1,0)


#number of tagged sharks at start of year
#Dummy
N.start=c(300,400)

#non-reporting rate per region and month
D=round(matrix(runif(n.reg*n.mon,0,0.5),ncol=n.reg,nrow=n.mon),2)

#annual tag shedding
S=0.0358  #Simpfendorfer et al (1999)


#PROCEDURE

#1.Apply non-reporting rate to each region
R.c=R
R.c[,,1]=R.c[,,1]/(1-D)
R.c[,,2]=R.c[,,2]/(1-D)


#2.Calculate number of tags of cohort a caught in month t corrected by non-reporting and tag shedding
R.a.t=apply(R.c,1,colSums)/exp(-S/12)



#3. Calculate exploitation rate
P=rowSums(R.a.t)



#4. Add exploitation rate to survival probability
l.x=l.x(1-P)exp-M





#McAuley et al (2005, 2007)           #missing:   USE REAL DATA AND ADAPT TO THESE
#DATA
#number of regions
n.reg=3
#number of months
n.mon=12
#number of tagged cohorts
n.coho=2

#number of tags reported from cohort a in region i during month t
#dummy
n.rep=round(runif(n.reg*n.mon*n.coho,0,20))
R=array(n.rep,c(n.mon,n.reg,n.coho))
dimnames(R)[[2]]=c('Region1','Region2','Region3')
dimnames(R)[[3]]=c('cohort1','cohort2')

#number of new tags of cohort a at start of month in each region
#Dummy
N.start=array(round(runif(n.reg*n.mon*n.coho,50,100)),c(n.mon,n.reg,n.coho))
dimnames(N.start)=dimnames(R)

#totla number of new tags of cohort a at month t 
N.start=apply(N.start,1,colSums)


#non-reporting rate per region and month
#Dummy
D=round(matrix(runif(n.reg*n.mon,0,0.5),ncol=n.reg,nrow=n.mon),2)

#annual tag shedding
S=0.0358  #Simpfendorfer et al (1999)

#natural mortality
M=0.1

#PROCEDURE

#1.Apply non-reporting rate to each region
R.c=R
R.c[,,1]=R.c[,,1]/(1-D)
R.c[,,2]=R.c[,,2]/(1-D)

R.a.t=apply(R.c,1,colSums)

#2.Calculate number of tagged sharks of cohort a present in population at start of month, corrected for shedding
N.t=R.a.t
N.t[,]=NA

N.t[1:2,1]=N.start[1:2,1]

for(j in 1:n.coho)
{
  for ( i in 2:n.mon)
  {
    N.t[j,i]=(N.t[j,i-1]-R.a.t[j,i-1])*exp(-(M+S)/12)+N.start[j,i]
  }
}

#3. Calculate instantaneous F rate by cohort a in year T by solving Baranov catch equation       #INCOMPLETE!!!!!!
#note: catch equation has no analytical solution, hence use numerical methods (e.g. Martell 2010)
C.pred=(N*F*(1-exp(-(M+F))))/(M+F)

epsilon=(R.a.t-C.pred)^2


#4. Use F in proportion surviving in demographic analysis
l.x=lx-1(1-Fx-1)e(-Mx-1)