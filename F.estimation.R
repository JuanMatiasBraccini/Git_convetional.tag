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