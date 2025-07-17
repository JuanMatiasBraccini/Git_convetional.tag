# Estimate F from tag-recapture data

#NOTES: this script implements a range of methods for estimating F from mark-recapture
#note: only estimate F for each release year to minimize effect of movement out of release area
#      only use animals released within fishing ground
#      only use animals released in conditions 1 and 2 to minimize mortality due to tagging

rm(list=ls(all=TRUE))
library(fishmethods)
library(tidyverse)
library(matrixStats)
library(triangle,quietly=T)
library(reshape)

if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')

#---DATA SECTION-------

#1. tagging data
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_conventional_data.R"))

LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))
Selectivity.pars=data.frame(SPECIES=c("TK","BW","WH","GM"),
                            theta1=c(135.5,130.1,173.7,184.3),
                            theta2=c(117695,29237,26415,29739))

Catch.proportion.by.zone=data.frame(Species=c("TK","BW","WH","GM"),  #from daily logbook data, all years combined
                                    South=c(0.0597,0.248,0.648,0.931),
                                    South.west=c(0.197,0.626,0.303,0.0643),
                                    West= c(0.743,0.126,0.0492,0.00435))

Reported.F.McAuley.2007=list(TK=data.frame(Finyear=rep(2001:2003,times=13),
                                           Age=rep(c(seq(0,15,3),18,25:30),each=3))%>%
                               arrange(Finyear)%>%
                               mutate(f=c(.045,.07,.19,.07,.04,.025,rep(0,7),
                                          .03,.1,.09,.01,.015,.014,.013,.013,.011,.01,.005,.004,0,
                                          0,.055,.28,.075,.04,.03,.039,.03,.028,.02,.01,.005,0),
                                      Age.mid.point=rep(c(1.5,4.5,7.5,10.5,13.5,16.5,21,25:30),times=3)),
                             BW=data.frame(Finyear=rep(1994:1995,times=11),
                                           Age=rep(0:10,each=2))%>%
                               arrange(Finyear)%>%
                               mutate(f=c(.21,.145,.05,.025,.025,.022,.02,.015,.03,.01,0,
                                          .17,.09,.04,.025,.02,.012,.01,.012,.015,.008,0),
                                      Age.mid.point=Age))
  
  
#---PARAMETERS SECTION-------
min.days.liberty=30
sd.vonB=10  #assumed SD of vonB function (no SD provided by Simpfendorfer et al 2000)

released.in.fishing.grounds=recaptured.in.fishing.grounds=TRUE
tag.shedding=0.0358 #McAuley et al 2007 (yrs-1)
tag.non.reporting=data.frame(Finyear=c(1994:2003,2001:2003),
                             Species=c(rep('BW',10),rep('TK',3)),
                             South=c(.16,.14,.29,.34,.24,.18,.25,.17,.22,.4,rep(NA,3)),
                             South.west=c(.39,.34,.54,.61,.44,.54,.49,.58,.76,.82,.23,.13,.09),
                             West=c(.15,.1,.25,.43,.46,.6,.25,.38,.32,.52,.5,.47,.67),
                             North=c(rep(NA,10),.39,.55,.84))  #McAuley et al 2007 
# tag.non.reporting=tag.non.reporting%>%
#                 mutate(Keep=ifelse((Finyear%in%1994:2003 & Species%in%c('BW','GM','WH'))|
#                                      (Finyear%in%2001:2003 & Species%in%c('TK')),'yes','no'))%>%
#                 filter(Keep=='yes')

#assume whiskery and gummy have same non-reporting rate of duskies (conservative)
assumed.tag.non.reporting=tag.non.reporting%>%filter(Species=='BW')
nn=nrow(assumed.tag.non.reporting)
assumed.tag.non.reporting=rbind(assumed.tag.non.reporting,
                                assumed.tag.non.reporting)%>%
                        mutate(Species=c(rep('WH',nn),rep('GM',nn)))
tag.non.reporting=rbind(tag.non.reporting,assumed.tag.non.reporting)  

#export reporting and shedding
Dis.Sp=c("BW","TK","WH","GM")
names(Dis.Sp)=c('Dusky shark','Sandbar shark','Whiskery shark','Gummy shark')
hndl.dat.out="C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/Data_outs"
for(i in 1:length(Dis.Sp))
{
  a=tag.non.reporting%>%filter(Species==Dis.Sp[i])
  NmS=names(Dis.Sp)[i]
  write.csv(a,paste0(hndl.dat.out,'/',NmS,'/',paste0(NmS,"_","Con_tag_non_reporting_from_F.estimation.R_.csv")),row.names=F)
  write.csv(tag.shedding,paste0(hndl.dat.out,'/',NmS,'/',paste0(NmS,"_","Con_tag_shedding_from_F.estimation.R_.csv")),row.names = F)
}



#released years used in analysis
used.released.years=tag.non.reporting%>%
                    mutate(Keep=ifelse((Finyear%in%1994:1995 & Species%in%c('BW','GM','WH'))|
                                       (Finyear%in%2001:2003 & Species%in%c('TK')),'yes','no'))%>%
                  filter(Keep=='yes')
these.acous.yrs=2012:2014

SPECIES=c("TK","BW","WH","GM")
names(SPECIES)=c(18007,18003,17003,17001)
Species.names=c("sandbar shark","dusky shark","whiskery shark","gummy shark")
names(Species.names)=SPECIES

M.overall_minM=c(0.05,0.04,0.23,0.14) #from Assessment_all.species.R
M.at.sel.peak_minM=c(0.06,0.05,0.22,0.16) 
M.overall_meanM=c(0.15,0.11,0.34,0.28) 
M.at.sel.peak_meanM=c(0.19,0.2,0.31,0.28) 
names(M.overall_minM)=names(M.at.sel.peak_minM)=names(M.overall_meanM)=names(M.at.sel.peak_meanM)=SPECIES

MESH=16.5  #mesh used for selectivity schedule, in cm (6.5 inch)

F.type='Overall'

#---PROCEDURE SECTION-------
setwd(handl_OneDrive("Analyses/Conventional tagging/F_estimation"))

#1. Create useful objects
LH.data=LH.data%>%filter(SPECIES %in%c(18007,18003,17003,17001))

#keep only individuals tagged with conventional tags and released in conditions 1 or 2
Tagging.acous=Tagging%>%
  mutate(DATE_REL=as.Date( paste( Yr.rel,Mn.rel , Day.rel , sep = "-" ), format = "%Y-%m-%d" ),
         DATE_CAPTR=as.Date( paste( Yr.rec,Mn.rec , Day.rec , sep = "-" ), format = "%Y-%m-%d" ),
         DaysAtLarge=as.numeric(round((DATE_CAPTR-DATE_REL),0)),
         DaysAtLarge=ifelse(DaysAtLarge<0,NA,DaysAtLarge))%>%
  mutate(obs.len=Rel_FL)%>%
  filter(Yr.rel>=1993 & CONDITION%in%c(1,2))%>%
  filter(Tag.type=='acoustic' & Species%in%c('GM','WH'))%>%
  mutate(Fin.yr.rel=case_when(Mn.rel>=7~Yr.rel,
                              Mn.rel<7~Yr.rel-1),
         Fin.yr.rec=case_when(Mn.rec>=7~Yr.rec,
                              Mn.rec<7~Yr.rec-1))
#with(Tagging.acous,table(Fin.yr.rel,Fin.yr.rec,Species,useNA = 'ifany'))

Tagging=Tagging%>%
  mutate(DATE_REL=as.Date( paste( Yr.rel,Mn.rel , Day.rel , sep = "-" ), format = "%Y-%m-%d" ),
          DATE_CAPTR=as.Date( paste( Yr.rec,Mn.rec , Day.rec , sep = "-" ), format = "%Y-%m-%d" ),
          DaysAtLarge=as.numeric(round((DATE_CAPTR-DATE_REL),0)),
          DaysAtLarge=ifelse(DaysAtLarge<0,NA,DaysAtLarge))%>%
  mutate(obs.len=Rel_FL)%>%
  filter(Yr.rel>=1993 & CONDITION%in%c(1,2))%>%
  filter(Tag.type=='conventional' & Species%in%SPECIES)%>%
  mutate(Fin.yr.rel=case_when(Mn.rel>=7~Yr.rel,
                              Mn.rel<7~Yr.rel-1),
         Fin.yr.rec=case_when(Mn.rec>=7~Yr.rec,
                              Mn.rec<7~Yr.rec-1))
  

#Put data in list and predict age from observed length
List.tagging=vector('list',length(Species.names))
names(List.tagging)=Species.names
fn.age.len=function(A.max,Linf,K,to,LF_o,sd,obs.length,MAIN)
{
  age=0:A.max
  
  #von B function
  mean.length.3par=Linf*(1-exp(-K*(age-to)))
  mean.length=LF_o+(Linf-LF_o)*(1-exp(-K*age))
  
  sd=rep(sd,length(age))
  
  #prob of age given length
  prob=function(data)  dnorm(data,mean.length,sd)
  pred.len=sapply(obs.length,prob)
  
  id=rep(NA,ncol(pred.len))
  for (i in 1:ncol(pred.len))  id[i]=which.min(abs(pred.len[,i] - max(pred.len[,i])))  
  
  pred.age=age[id]
  
  plot(age,mean.length,type='l',lwd=2,main=MAIN)
  points(pred.age,obs.length,pch=19,col=2)
  
  return(pred.age=pred.age)
}
for(s in 1:length(List.tagging))
{
  d=Tagging%>%
    filter(Species==names(Species.names)[s])%>%
    filter(!is.na(Rel_FL))%>%
    filter((Rec.method!='Research longline')%>% replace_na(TRUE))   #remove records recaptured by scientific survey
  id=as.numeric(names(SPECIES)[match(names(Species.names)[s],SPECIES)])
  lh=LH.data%>%filter(SPECIES==id)
  lh$to=ifelse(is.na(lh$to),-1,lh$to)
  d$pred.age=fn.age.len(A.max=lh$Max_Age,
                        Linf=lh$FL_inf,
                        K=lh$K,
                        to=lh$to,
                        LF_o=lh$LF_o,
                        sd=sd.vonB,
                        obs.length=d$obs.len,
                        MAIN=Species.names[s])
  List.tagging[[s]]=d
}

#NSF (sandbar only)
List.tagging.NSF=list()
List.tagging.NSF[[1]]=List.tagging$`sandbar shark`%>%
                    filter(Lat.rels>(-26) & Long.rels>114)  #remove Ningaloo closure
                    
names(List.tagging.NSF)=names(List.tagging)[1]

#Acoustic tagging (gummy and whiskery only)
SPECIES.acous=subset(SPECIES,SPECIES%in%c("WH", "GM"))
List.tagging.acoustic=vector('list',length(SPECIES.acous))
names(List.tagging.acoustic)=Species.names[match(SPECIES.acous,names(Species.names))]

for(s in 1:length(List.tagging.acoustic))
{
  d=Tagging.acous%>%
    filter(Species==SPECIES.acous[s])%>%
    filter(!is.na(Rel_FL))%>%
    filter((Rec.method!='Research longline')%>% replace_na(TRUE))   #remove records recaptured by scientific survey
  id=as.numeric(names(SPECIES)[match(SPECIES.acous[s],SPECIES)])
  lh=LH.data%>%filter(SPECIES==id)
  lh$to=ifelse(is.na(lh$to),-1,lh$to)
  d$pred.age=fn.age.len(A.max=lh$Max_Age,
                        Linf=lh$FL_inf,
                        K=lh$K,
                        to=lh$to,
                        LF_o=lh$LF_o,
                        sd=sd.vonB,
                        obs.length=d$obs.len,
                        MAIN=Species.names[s])
  List.tagging.acoustic[[s]]=d
}



#TDGDLF fishing ground
if(released.in.fishing.grounds)
{
  for(s in 1:length(List.tagging))
  {
    List.tagging[[s]]=List.tagging[[s]]%>%
                        filter(Lat.rels<=(-26))
  }
}

#Recaptured in fishing ground
if(recaptured.in.fishing.grounds)
{
  for(s in 1:length(List.tagging))
  {
    List.tagging[[s]]=List.tagging[[s]]%>%
      filter((Lat.rec<=(-26))%>% replace_na(TRUE))
  }
}
#Visualize spatial distribution
fn.plot.spatial=function(d,Tit,yrs.selected)
{
  p=d%>%
    mutate(Labs=ifelse(Yr.rel%in%yrs.selected,paste(Yr.rel,'(used for F estim.)'),Yr.rel))%>%
    ggplot(aes(Long.rels,Lat.rels))+
    geom_point(color=1)+
    geom_point(aes(Long.rec,Lat.rec),color=2,alpha=0.45)+
    facet_wrap(~Labs)+
    ggtitle(Tit)
  print(p)
}
for(i in 1:length(List.tagging))
{
  fn.plot.spatial(d=List.tagging[[i]],
                  Tit=names(List.tagging)[i],
                  yrs.selected=used.released.years%>%filter(Species==SPECIES[i])%>%pull(Finyear))  
  ggsave(handl_OneDrive(paste('Analyses/Conventional tagging/F_estimation/Spatial dist_',
                              names(List.tagging)[i],"_TDGDLF.tiff",sep='')), 
         width = 10,height = 8,dpi = 300, compression = "lzw")  
}
for(i in 1:length(List.tagging.NSF))
{
  fn.plot.spatial(d=List.tagging.NSF[[i]],
                  Tit=names(List.tagging.NSF)[i],
                  yrs.selected=used.released.years%>%filter(Species==SPECIES[i])%>%pull(Finyear))
  ggsave(handl_OneDrive(paste('Analyses/Conventional tagging/F_estimation/Spatial dist_',
                              names(List.tagging)[i],"_NSF.tiff",sep='')), 
         width = 10,height = 8,dpi = 300, compression = "lzw")
}

#Select relevant years by species
for(s in 1:length(List.tagging))
{
  dis.rel.yr=used.released.years%>%filter(Species==unique(List.tagging[[s]]$Species))%>%pull(Finyear)
  List.tagging[[s]]=List.tagging[[s]]%>%
                            filter(Yr.rel%in%dis.rel.yr)
}

for(s in 1:length(List.tagging.NSF))
{
  dis.rel.yr=used.released.years%>%filter(Species==unique(List.tagging.NSF[[s]]$Species))%>%pull(Finyear)
  List.tagging.NSF[[s]]=List.tagging.NSF[[s]]%>%
    filter(Yr.rel%in%dis.rel.yr)
}

#---Get Kirkwood & Walker selectivity-------
Rango.sizes=vector('list',length(SPECIES))
names(Rango.sizes)=SPECIES
for(s in 1:length(List.tagging))
{
  Rango.sizes[[s]]=data.frame(SPECIES=SPECIES[s],
                              MIN.FL=min(List.tagging[[s]]$Rel_FL,na.rm=T),
                              MAX.FL=max(List.tagging[[s]]$CAP_FL,na.rm=T))
}
Rango.sizes=do.call(rbind,Rango.sizes)
these.FLs=seq(10*floor(min(Rango.sizes$MIN.FL)/10),ceiling(max(Rango.sizes$MAX.FL)/10)*10,1)
pred.Kirkwood.Walker=function(theta,pred.len,Mesh)
{
  Theta1=theta[1]
  Theta2=theta[2]
  Mesh=sort(Mesh)
  
  Mesh=round(Mesh*0.393701,1)  #convert cm to inches
  pred.len=pred.len *10        #length in mm
  
  d1=data.frame(Size.class=rep(pred.len,times=length(Mesh)),
                Mesh.size=rep(Mesh,each=length(pred.len)))%>%
    mutate(alpha.beta=Theta1*Mesh.size,
           beta=-0.5*((alpha.beta)-((alpha.beta*alpha.beta+4*Theta2)^0.5)),
           alpha=alpha.beta/beta,
           Rel.sel=((Size.class/(alpha*beta))^alpha)*(exp(alpha-(Size.class/beta))))%>%
    dplyr::select(-c('alpha.beta','beta','alpha'))%>%
    spread(Mesh.size,Rel.sel)
  return(d1)
}
Sel.pred=vector('list',length(SPECIES))
names(Sel.pred)=SPECIES
for(s in 1:length(SPECIES))
{
  Sel.pred[[s]]=pred.Kirkwood.Walker(theta=unlist(Selectivity.pars%>%
                                                    filter(SPECIES==SPECIES[s])%>%dplyr::select(theta1,theta2)),
                                     pred.len=these.FLs,
                                     Mesh=MESH)%>%
                        dplyr::rename(FL=Size.class)
  plot(Sel.pred[[s]]$FL,Sel.pred[[s]][,2],type='l',
       ylab='Relative selectivity',xlab="Fork length (mm)",las=1,
       cex.lab=2,cex.axis=1.5,col=2,lwd=4,main=names(Sel.pred)[s],ylim=c(0,1))
}

#---Calculate F. based on McAuley et al 2007-------
#note: incomplete... use Al Harry's instead
do.this=FALSE
if(do.this)
{
  fn.fishing.mort.McAuley=function(d,shedding,non.reporting,M,Ages.selected)  
  {
    #Fill in reported captures
    d.rec=d%>%filter(Recaptured=='Yes')%>%
      mutate(Lat.rec=ifelse(is.na(Lat.rec),Lat.rels,Lat.rec),
             Long.rec=ifelse(is.na(Long.rec),d$Long.rels,Long.rec),
             Recap.region=case_when(Lat.rec>(-33) & Lat.rec<(-21) & Long.rec<(116.5)~'West',
                                    Lat.rec>(-23) & Long.rec>(114)~'North',
                                    Lat.rec<=(-33) & Long.rec<(116.5)~'South.west',
                                    Lat.rec<=(-28) & Long.rec>=(116.5)~'South'))
    Regions=unique(d.rec$Recap.region)
    n.reg=length(Regions)   
    Months=1:12
    n.mon=length(Months)
    AgeS=sort(unique(d$pred.age))
    n.age=length(AgeS)
    Captures=array(0,c(n.age,n.reg,n.mon))
    dimnames(Captures)[[3]]=1:n.mon
    dimnames(Captures)[[2]]=Regions
    dimnames(Captures)[[1]]=paste('Age',AgeS)
    for(a in 1:n.age)
    {
      for(r in 1:n.reg)
      {
        for(m in 1:n.mon)
        {
          xx=d.rec%>%filter(pred.age==AgeS[a] & Recap.region==Regions[r] & Mn.rec==Months[m])
          if(nrow(xx)>0)
          {
            Captures[match(AgeS[a],AgeS),match(Regions[r],Regions),match(Months[m],Months)]=nrow(xx)
          }
        }
      }
    }
    
    #Calculate predicted Captures
    Pred.captures=matrix(0,nrow=n.age,ncol=n.mon)
    rownames(Pred.captures)=AgeS
    colnames(Pred.captures)=Months
    for(a in 1:n.age)
    {
      for(m in 1:n.mon)
      {
        non.rep=non.reporting%>%
          filter(Finyear==unique(d$Yr.rel))
        non.rep=non.rep[,match(Regions,names(non.reporting))]
        Pred.captures[a,m]=sum(Captures[a,,m]/(1-non.rep))
      }
    }
    Pred.captures=rowSums(Pred.captures)
    
    #Numbers of tagged individuals
    Releases=matrix(0,nrow=n.age,ncol=n.mon)
    rownames(Releases)=AgeS
    colnames(Releases)=Months
    N=Releases
    for(a in 1:n.age)
    {
      for(m in 1:n.mon)
      {
        xx=d%>%filter(pred.age==AgeS[a] & Mn.rel==Months[m])
        Releases[a,m]=nrow(xx)
      }
    }
    
    for(a in 1:n.age)
    {
      iid=which(Releases[a,]>0)[1]
      N[a,iid]=Releases[a,iid]
      iid=iid+1
      for(m in iid:n.mon)
      {
        if(m==iid) plus=0
        if(m>iid) plus=Releases[a,m-1]
        N[a,m]=(N[a,m-1]-Pred.captures[a,m-1])*exp(-((M[a]+shedding)/12))+plus
      }
    }
    
    #Predicted annual catch
    F.at.age=M*1.5
    Zeta.at.age=M+F.at.age
    N.at.age=rowSums(N)
    Pred.annual.catch=(F.at.age/Zeta.at.age)*N.at.age*(1-exp(-Zeta.at.age))
    
    epsilon=(Pred.captures-Pred.annual.catch)^2
    
  }
  
  #TDGDLF   
  #note: the low number of recaptures and the aggregation of ages yields a very low F
  #       artificially increasing the number of recaptures yields higher F
  Tagging.model=function(Ft,Recaptured.hat,eMe,tag.shedding,Tagged)
  {
    Nt=matrix(0,nrow=1,ncol=ncol(Recaptured.hat))
    iid=which(Tagged>0)[1]
    iid=iid+1
    for(t in iid:ncol(Nt))Nt[t]=(Nt[t-1]-Recaptured.hat[t-1])*exp(-(eMe+tag.shedding)/12)+Tagged[t-1]
    
    Zt=Ft+eMe  
    Pred.Recaptured=(Ft/Zt)*Nt*(1-exp(-Zt))
    
    id=which(Recaptured.hat==0)
    res=sum((Recaptured.hat[-id]-Pred.Recaptured[-id])^2)
    
    return(list(res=res,Pred.Recaptured=Pred.Recaptured,Recaptured.hat=Recaptured.hat))
    
  }
  ObjFunc=function(par)
  {
    return(Tagging.model(par,Recaptured.hat=Recaptured.hat,eMe=eMe,tag.shedding=tag.shedding,Tagged=Tagged)$res)
  }
  params=c(Ft=0.1)
  F.McAuley.2007=vector('list',length(Species.names))
  names(F.McAuley.2007)=Species.names
  for(s in 1:length(SPECIES))
  {
    unik.yrs=sort(unique(List.tagging[[s]]$Yr.rel))
    
    non.reporting=tag.non.reporting%>%
      filter(Species==SPECIES[s])%>%
      dplyr::select(-North)%>%
      select_if(~ !any(is.na(.)))
    
    F.by.year=vector('list',length(unik.yrs))
    names(F.by.year)=unik.yrs
    for(y in 1:length(unik.yrs))
    {
      dat=List.tagging[[s]]%>%
        filter(Yr.rel==unik.yrs[y])
      
      #Ages combined to estimate an overall F
      if(F.type=='Overall')
      {
        d=dat
        eMe=list.M$M.at.sel.peak_meanM[match(SPECIES[s],names(list.M$M.at.sel.peak_meanM))]
        D=non.reporting%>%filter(Finyear==unik.yrs[y])
        
        #get recapture data
        d.rec=d%>%filter(Recaptured=='Yes')%>%
          mutate(Lat.rec=ifelse(is.na(Lat.rec),Lat.rels,Lat.rec),
                 Long.rec=ifelse(is.na(Long.rec),d$Long.rels,Long.rec),
                 Recap.region=case_when(Lat.rec>(-33) & Lat.rec<(-21) & Long.rec<(116.5)~'West',
                                        Lat.rec>(-23) & Long.rec>(114)~'North',
                                        Lat.rec<=(-33) & Long.rec<(116.5)~'South.west',
                                        Lat.rec<=(-28) & Long.rec>=(116.5)~'South'))%>%
          filter(Yr.rec==unik.yrs[y])
        
        #get dimensions
        Regions=unique(d.rec$Recap.region)
        n.reg=length(Regions)   
        Months=1:12
        n.mon=length(Months)
        n.age=1
        
        #fill in total tags releases
        Tagged=matrix(0,nrow=n.age,ncol=n.mon)
        rownames(Tagged)='Combined.Age'
        colnames(Tagged)=Months
        for(m in 1:n.mon)
        {
          xx=d%>%filter(Mn.rel==Months[m])
          Tagged[n.age,m]=nrow(xx)
        }
        
        #fill in captures
        Captures=array(0,c(n.age,n.reg,n.mon))
        dimnames(Captures)[[3]]=1:n.mon
        dimnames(Captures)[[2]]=Regions
        dimnames(Captures)[[1]]='Combined.Age'
        for(a in 1:n.age)
        {
          for(r in 1:n.reg)
          {
            for(m in 1:n.mon)
            {
              xx=d.rec%>%filter(Recap.region==Regions[r] & Mn.rec==Months[m])
              if(nrow(xx)>0)
              {
                Captures[1,match(Regions[r],Regions),match(Months[m],Months)]=nrow(xx)
              }
            }
          }
        }
        
        #sum captures over regions
        Recaptured=matrix(0,nrow=n.age,ncol=n.mon)
        for(a in 1:n.age) for(m in 1:n.mon)  Recaptured[a,m]=sum(Captures[a,,m])
        
        #fill in expected captures given non-reporting
        Recaptured.hat=matrix(0,nrow=n.age,ncol=n.mon)
        rownames(Recaptured.hat)='Combined.Age'
        colnames(Recaptured.hat)=Months
        for(a in 1:n.age)
        {
          for(m in 1:n.mon)
          {
            Dm=D[,match(Regions,names(D))]
            Recaptured.hat[a,m]=sum(Captures[a,,m]/(1-Dm))
          }
        }
        
        #Estimate overall F
        nlmb <- nlminb(params, ObjFunc, gradient = NULL, hessian = TRUE)
        fit = optim(nlmb$par, ObjFunc,method='Brent',lower = 1e-6, upper = 5,hessian = TRUE)
        se <- sqrt(diag(solve(fit$hessian)))
        
        
        
        #likelihood profile
        dis=seq(0,1,by=0.001)
        plot(dis,sapply(dis,ObjFunc))
        
        
        #plot observed vs expected
        PREDS=Tagging.model(Ft=fit$par,Recaptured.hat=Recaptured.hat,eMe=eMe,tag.shedding=tag.shedding,Tagged=Tagged)
        plot(1:length(PREDS$Pred.Recaptured),PREDS$Recaptured.hat,xlab='Month',ylab='Recaptures',pch=19,
             main=paste(Species.names[match(SPECIES[s],names(Species.names))],' (year=',unik.yrs[y],
                        ', tagged=',sum(Tagged),', recaptured=',sum(Recaptured),')',sep=''))
        points(1:length(PREDS$Pred.Recaptured),PREDS$Pred.Recaptured,col=2)
        legend('topright',c('Recaptured.hat','Pred.Recaptured'),bty='n',pch=c(19,1),col=1:2)
        
        out=data.frame(Finyear=unik.yrs[y],Mean=fit$par,CV=se/fit$par)
        
        #clean log
        rm(d,D,eMe,Captures,Recaptured,Recaptured.hat,Tagged,d.rec,Regions)
      }
      
      if(F.type=='F.at.age') #not finished
      {
        #pool ages
        if(SPECIES[s]=='TK') pooled.ages=c(seq(0,12,3))
        if(SPECIES[s]=='BW') pooled.ages=1:3
        if(SPECIES[s]=='WH') pooled.ages=c(4:6)
        if(SPECIES[s]=='GM') pooled.ages=c(5:9)
        
        dat=dat%>%
          filter(pred.age<=max(pooled.ages))%>%
          mutate(pred.age.pool=cut(pred.age,breaks=pooled.ages),
                 pred.age=as.numeric(as.character(str_extract(pred.age.pool,"[0-9]+"))))
        eigs=sort(unique(dat%>%pull(pred.age)))
        M.at.age=rep(list.M$M.at.sel.peak_meanM[match(SPECIES[s],names(list.M$M.at.sel.peak_meanM))],length(eigs)) 
        
        
        for(e in 1:length(eigs))
        {
          AgeS=eigs[e]
          d=dat%>%filter(pred.age==AgeS)
          D=non.reporting%>%filter(Finyear==unik.yrs[y])
          eMe=M.at.age[e]
          
          #get recapture data
          d.rec=d%>%filter(Recaptured=='Yes')%>%
            mutate(Lat.rec=ifelse(is.na(Lat.rec),Lat.rels,Lat.rec),
                   Long.rec=ifelse(is.na(Long.rec),d$Long.rels,Long.rec),
                   Recap.region=case_when(Lat.rec>(-33) & Lat.rec<(-21) & Long.rec<(116.5)~'West',
                                          Lat.rec>(-23) & Long.rec>(114)~'North',
                                          Lat.rec<=(-33) & Long.rec<(116.5)~'South.west',
                                          Lat.rec<=(-28) & Long.rec>=(116.5)~'South'))%>%
            filter(Yr.rec==unik.yrs[y])
          #get dimensions
          Regions=unique(d.rec$Recap.region)
          n.reg=length(Regions)   
          Months=1:12
          n.mon=length(Months)
          AgeS=sort(unique(d$pred.age))
          n.age=length(AgeS)
          
          #fill in total tags releases
          Tagged=matrix(0,nrow=n.age,ncol=n.mon)
          rownames(Tagged)=AgeS
          colnames(Tagged)=Months
          for(a in 1:n.age)
          {
            for(m in 1:n.mon)
            {
              xx=d%>%filter(pred.age==AgeS[a] & Mn.rel==Months[m])
              Tagged[a,m]=nrow(xx)
            }
          }
          
          #fill in captures
          Captures=array(0,c(n.age,n.reg,n.mon))
          dimnames(Captures)[[3]]=1:n.mon
          dimnames(Captures)[[2]]=Regions
          dimnames(Captures)[[1]]=paste('Age',AgeS)
          for(a in 1:n.age)
          {
            for(r in 1:n.reg)
            {
              for(m in 1:n.mon)
              {
                xx=d.rec%>%filter(pred.age==AgeS[a] & Recap.region==Regions[r] & Mn.rec==Months[m])
                if(nrow(xx)>0)
                {
                  Captures[match(AgeS[a],AgeS),match(Regions[r],Regions),match(Months[m],Months)]=nrow(xx)
                }
              }
            }
          }
          
          #sum captures over regions
          Recaptured=matrix(0,nrow=n.age,ncol=n.mon)
          for(a in 1:n.age) for(m in 1:n.mon)  Recaptured[a,m]=sum(Captures[a,,m])
          
          #fill in expected captures given non-reporting
          Recaptured.hat=matrix(0,nrow=n.age,ncol=n.mon)
          rownames(Recaptured.hat)=AgeS
          colnames(Recaptured.hat)=Months
          for(a in 1:n.age)
          {
            for(m in 1:n.mon)
            {
              D=D[,match(Regions,names(D))]
              Recaptured.hat[a,m]=sum(Captures[a,,m]/(1-D))
            }
          }
          rm(d,D,eMe,Captures,Recaptured,Recaptured.hat,Tagged,AgeS,d.rec,Regions)
        }
      }
      
      F.by.year[[y]]=out
    }
    
    F.McAuley.2007[[s]]=do.call(rbind,F.by.year)
  }
}

#---Get F. estimated by McAuley et al 2007 for dusky and sandbar-------
Annual.Reported.F.McAuley.2007=Reported.F.McAuley.2007
for(s in 1:length(Reported.F.McAuley.2007))
{
  id=as.numeric(names(SPECIES)[match(names(Reported.F.McAuley.2007)[s],SPECIES)])
  lh=LH.data%>%filter(SPECIES==id)
  lh$to=ifelse(is.na(lh$to),-1,lh$to)
  
  #get length (in mm) from age
  AgE=Reported.F.McAuley.2007[[s]]$Age.mid.point
  if(SPECIES[s]=='BW') AgE=AgE+2 #to allow max selectivity for 0+
  Reported.F.McAuley.2007[[s]]$mean.length=10*(lh$LF_o+(lh$FL_inf-lh$LF_o)*(1-exp(-lh$K*(AgE))))
  
  
  #get selectivity
  Sel=Sel.pred[[match(names(Reported.F.McAuley.2007)[s],names(Sel.pred))]]
  
  
  Reported.F.McAuley.2007[[s]]$Selectivity=NA
  for(pp in 1:nrow(Reported.F.McAuley.2007[[s]]))
  {
    Reported.F.McAuley.2007[[s]]$Selectivity[pp]=Sel[which.min(abs(Sel$FL-Reported.F.McAuley.2007[[s]]$mean.length[pp])),2]
  }
  
  #get F
  F.approach='max'  #chose F at max selectivity
   #max F approach
  if(F.approach=='max')
  {
    Annual.Reported.F.McAuley.2007[[s]]=Reported.F.McAuley.2007[[s]]%>%
      filter(f>0)%>%
      group_by(Finyear)%>%
      summarise(Mean=max(f))%>%
      ungroup()
    
  }

  
  #F derived from F@age and selectivity
  if(F.approach=='F@age')
  {
    ierS=unique(Reported.F.McAuley.2007[[s]]$Finyear)
    for(pp in 1:length(ierS))
    {
      dd=Reported.F.McAuley.2007[[s]]%>%
                filter(Finyear==ierS[pp])
      mod <- loess(f ~ Age, dd)
      dd=dd%>%
          mutate(loessSel=predict(mod,Age),
                 F=f/Selectivity)
    }
  }
 }


#---Calculate F. based on Harry et al 2016-------  
source(handl_OneDrive('Analyses/Conventional tagging/F_estimation/Harry et al 2016 F estimation/Fishing Mortality Functions V2.R'))
MC=1e2
F.Harry=List.tagging
for(s in 1:length(List.tagging))
{
  disyiers=used.released.years%>%filter(Species==SPECIES[s])%>%pull(Finyear)
  disyiers=c(disyiers,disyiers[length(disyiers)]+1)
  d=List.tagging[[s]]%>%
    filter(Yr.rel%in%disyiers)  
  Deployed.tags=d%>%
                  dplyr::select(Yr.rel,Mn.rel)%>%
                  dplyr::rename(Year=Yr.rel,
                         Month=Mn.rel)
  Recaptured.tags=d%>%
                  filter(Recaptured=='Yes' & Rec.method=='Commercial gillnet')%>%
                  dplyr::select(Yr.rec,Mn.rec,DaysAtLarge)%>%
                  dplyr::rename(R.Year=Yr.rec,
                                R.Month=Mn.rec,
                                Days.Out=DaysAtLarge)%>%
                  mutate(R.Sector="ComNet")
  
  nn=as.numeric(names(SPECIES[s]))
  LH=LH.data%>%filter(SPECIES==nn)
  
  Reporting=tag.non.reporting%>%
                   filter(Species==SPECIES[s])%>%
                    dplyr::select(-North)%>%
                    select_if(~ !any(is.na(.)))%>%
                    gather(Region,Non.reporting,-c(Finyear,Species))%>%
                    mutate(Rep=1-Non.reporting)%>%
                    group_by(Finyear)%>%
                    summarise(Rep.rate=mean(Rep),
                              Rep.rate.sd=sd(Rep))  

  
  dd=vector('list',MC)
  f_matrix=dd
  for(xx in 1:length(dd))
  {
   
    life.history=list(M.type="Then",
                      amax=rtriangle(1,LH$Max_Age,LH$Max_Age_max,LH$Max_Age),
                      mat_a50=rtriangle(1,LH$Age_50_Mat_min,LH$Age_50_Mat_max,LH$Age_50_Mat_min),
                      k=rlnorm(n=1,log(LH$K),0.1))
    Reporting.sample=Reporting%>%
                      mutate(Rep.samp=rlnorm(nrow(Reporting),log(Rep.rate),Rep.rate.sd),
                             Rep.samp=ifelse(Rep.samp>0.99,0.99,Rep.samp))  
    
    
    dd[[xx]]=GetF(TAGS=Deployed.tags,
                   RECAPS=Recaptured.tags,
                   LHP=life.history,
                   Reporting=Reporting.sample,
                   Shedding=tag.shedding,
                   MIX=min.days.liberty)
    print(paste(SPECIES[s],xx,sep='  ---- MC iteration '))
    
    dummi=dd[[xx]]%>%
                dplyr::select(Fxt,Finyear)%>%
                group_by(Finyear)%>%
                summarise(Fxt=sum(Fxt,na.rm=T))%>%
                spread(Finyear,Fxt)
    
    f_matrix[[xx]]<-dummi
  }

  f_matrix=do.call(rbind,f_matrix)
  meltData <- melt(as.data.frame(f_matrix))
  boxplot(data=meltData, value~variable,main=SPECIES[s],ylim=c(0,1))
  
  F.Harry[[s]]=f_matrix%>%
                gather(Finyear,Mean)%>%
                group_by(Finyear)%>%
                summarise(SD=sd(Mean),
                          Mean=mean(Mean))%>%
                ungroup()%>%
                mutate(CV=SD/Mean)%>%
                dplyr::select(Finyear,Mean,CV)
}



#---Calculate F. Hoenig et al 1998 ------
#note: age-independent instantaneous rates model of Hoenig et al. (1998)    

#modify the irm_h function to fix M
irm_h_fixedM=function (relyrs = NULL, recapyrs = NULL, N = NULL, recapharv = NULL, 
                       lambda = NULL, phi = NULL, Fyr = NULL, Myr = NULL, initial = NULL, 
                       lower = c(1e-04, 1e-04), upper = c(5, 5), maxiter = 10000,M.estim=FALSE) 
{
  Fp <- rep(initial[1], length(Fyr))
  Mp <- rep(initial[2], length(Myr))
  if (!is(recapharv, "matrix")) 
    stop("recapharv is not a matrix.")
  if (is.null(relyrs) | is.null(recapyrs)) 
    stop("Missing relyrs or recapyrs.")
  if (is.null(N)) 
    stop("Ns are missing.")
  if (is.null(recapharv)) 
    stop("Missing recovery matrix for harvested fish.")
  if (is.null(lambda)) 
    stop("lambdas (reporting rates) for harvested fish) are missing.")
  if (is.null(phi)) 
    stop("hphi (initial tag survival rates) for harvested fish) are missing.")
  if (is.null(Fyr)) 
    stop("Year designations for fishing mortality estimates (Fyr) are missing.")
  if (is.null(Myr)) 
    stop("Year designations for natural mortality estimates (Myr) are missing.")
  if (Myr[1] != relyrs[1] | Myr[1] != recapyrs[1]) 
    stop("First year in Myr should be equal to the first year in relyrs and recapyrs.")
  if (Fyr[1] != relyrs[1] | Fyr[1] != recapyrs[1]) 
    stop("First year in Fyr should be equal to the first year in relyrs and recapyrs.")
  if (relyrs[1] != recapyrs[1]) 
    stop("First year in relyrs and recapyrs should be the same.")
  rec <- length(seq(recapyrs[1], recapyrs[2], 1))
  if (ncol(recapharv) != rec) 
    stop("The number of columns in recapharv does not equal the number of\n             year specified by recapyrs.")
  if (length(lambda) != rec) 
    stop("The number of values in lambda does not equal the number of\n             year specified by recapyrs.")
  if (length(phi) != rec) 
    stop("The number of values in phi does not equal the number of\n             year specified by recapyrs.")
  rel <- length(seq(relyrs[1], relyrs[2], 1))
  if (length(N) != rel) 
    stop("The number of values in N does not equal the number of\n             year specified by relyrs.")
  if (initial[1] < lower[1] | initial[1] > upper[1]) 
    stop("initial F is outside lower or upper bounds.")
  if (initial[2] < lower[2] | initial[2] > upper[2]) 
    stop("initial M is outside lower or upper bounds.")
  if(M.estim)parms <- c(Fp, Mp) else parms <- Fp
  fit.obs <- function(x) {
    F <- x[1:length(Fp)]
    if(M.estim)
    {
      M <- x[as.numeric(length(Fp) + 1):as.numeric(length(Fp) + length(Mp))]
    }else
    {
      M <- rep(initial[2], length(Myr))
    }
    
    styrR <- 1
    endyrR <- relyrs[2] - relyrs[1] + 1
    styr <- 1
    endyr <- recapyrs[2] - recapyrs[1] + 1
    tags <- NULL
    cnt <- 0
    for (t in styrR:endyrR) {
      Ntags <- 0
      for (y in as.numeric(styr + cnt):endyr) {
        Ntags <- Ntags + recapharv[t, y]
      }
      tags[t] <- Ntags
      cnt <- cnt + 1
    }
    yrvector <- seq(recapyrs[1], recapyrs[2], 1)
    Fvector <- rep(NA, length(yrvector))
    Ftemp <- Fyr
    Ftemp[length(Ftemp) + 1] <- recapyrs[2] + 1
    for (t in styr:endyr) {
      for (d in 1:length(F)) {
        if (yrvector[t] >= Ftemp[d] & yrvector[t] < Ftemp[d + 
                                                          1]) 
          Fvector[t] <- F[d]
      }
    }
    Mvector <- rep(NA, length(yrvector))
    Mtemp <- Myr
    Mtemp[length(Mtemp) + 1] <- recapyrs[2] + 1
    for (t in styr:endyr) {
      for (d in 1:length(M)) {
        if (yrvector[t] >= Mtemp[d] & yrvector[t] < Mtemp[d + 
                                                          1]) 
          Mvector[t] <- M[d]
      }
    }
    s <- array(NA, dim = c(endyrR, endyr))
    cnt <- 0
    for (t in styrR:endyrR) {
      for (y in as.numeric(styr + cnt):endyr) {
        if (t == y) {
          s[t, y] <- 1
        }
        if (t != y) {
          s[t, y] <- exp(-Fvector[y - 1] - Mvector[y - 
                                                     1])
        }
      }
      cnt <- cnt + 1
    }
    u_h <- array(NA, dim = c(endyrR, endyr))
    cnt <- 0
    for (t in styrR:endyrR) {
      for (y in as.numeric(styr + cnt):endyr) {
        u_h[t, y] <- (Fvector[y]/(Fvector[y] + Mvector[y])) * 
          (1 - exp(-Fvector[y] - Mvector[y]))
      }
      cnt <- cnt + 1
    }
    cnt <- 0
    s_prob <- array(NA, dim = c(endyrR, endyr))
    for (t in styrR:endyrR) {
      looper <- 0
      for (y in as.numeric(styr + cnt):endyr) {
        probs <- 1
        for (a in as.numeric(y - looper):y) {
          probs <- probs * s[t, a]
        }
        s_prob[t, y] <- probs
        looper <- looper + 1
      }
      cnt <- cnt + 1
    }
    exp_prob_h <- array(NA, dim = c(endyrR, endyr))
    sum_prob_h <- NULL
    cnt <- 0
    for (t in styrR:endyrR) {
      dodo <- 0
      for (y in as.numeric(styr + cnt):endyr) {
        exp_prob_h[t, y] <- lambda[y] * phi[y] * s_prob[t, 
                                                        y] * u_h[t, y]
        dodo <- dodo + exp_prob_h[t, y]
      }
      sum_prob_h[t] <- dodo
      cnt <- cnt + 1
    }
    ll_h <- array(NA, dim = c(endyrR, endyr))
    ll_ns <- NULL
    cnt <- 0
    for (t in styrR:endyrR) {
      for (y in as.numeric(styr + cnt):endyr) {
        ll_h[t, y] <- 0
        if (recapharv[t, y] != 0) {
          ll_h[t, y] <- recapharv[t, y] * log(exp_prob_h[t, 
                                                         y])
        }
      }
      cnt <- cnt + 1
    }
    for (t in styrR:endyrR) {
      ll_ns[t] <- (N[t] - tags[t]) * log(1 - sum_prob_h[t])
    }
    f_tag <- 0
    cnt <- 0
    for (t in styrR:endyrR) {
      for (y in as.numeric(styr + cnt):endyr) {
        f_tag <- f_tag + ll_h[t, y]
      }
      cnt <- cnt + 1
    }
    for (t in styrR:endyrR) {
      f_tag <- f_tag + ll_ns[t]
    }
    LL <- f_tag * -1
    LL
  }
  low <- c(rep(lower[1], length(Fyr)), rep(lower[2], length(Myr)))
  up <- c(rep(upper[1], length(Fyr)), rep(upper[2], length(Myr)))
  results <- optim(parms, fit.obs, gr = NULL, lower = low, 
                   upper = up, method = c("L-BFGS-B"), control = list(maxit = maxiter), 
                   hessian = TRUE)
  var <- diag(solve(results$hessian))
  varcov <- solve(results$hessian)
  cormat <- (cov2cor(varcov))
  if(M.estim)
  {
    corr <- list(c(paste("F", 1:length(Fp), sep = ""), 
                   paste("M", 1:length(Mp), sep = "")))
    
  }else
  {
    corr <- list(c(paste("F", 1:length(Fp), sep = "")))
    
  }
  dimnames(cormat)[1] <- corr
  dimnames(cormat)[2] <- corr
  
  Fishing.mortality=data.frame(Yr=Fyr,
                               F=results$par)%>%
                                mutate(SE=sqrt(var),
                                       CV=SE/F)
                               
  return(list(model=results,Fishing.mortality=Fishing.mortality))
}
fn.fit.Z.estim=function(d,M.fixed,Fishery,MAIN)
{
  #Create input files
  SPP=unique(d$Species)
  Non.reporting=tag.non.reporting%>%
    filter(Species==SPP)%>%
    arrange(Finyear)
  Weits=Catch.proportion.by.zone%>%filter(Species==SPP)
  tag.shed=tag.shedding
  if(Fishery=='TDGDLF')
  {
    Non.reporting$Mean.non.rep=rowWeightedMeans(as.matrix(Non.reporting[,c('South','South.west','West')]),
                                                w=c(Weits$South,Weits$South.west,Weits$West), na.rm = T) #weighted average
  }
  if(Fishery=='NSF')
  {
    Non.reporting$Mean.non.rep=rowMeans(as.matrix(Non.reporting[,c('North')]),na.rm = T) 
  }
  if(Fishery=="TDGDLF - acoustic")
  {
    dummy=unique(d$Yr.rec)
    dummy=subset(dummy,dummy!='NA')
    dummy=sort(unique(c(unique(d$Yr.rel),dummy)))
    dummy=min(dummy):max(dummy)
    Non.reporting=data.frame(Finyear=dummy,
                             Mean.non.rep=0)
    tag.shed=0
  }
  
  dis.rec.yr=Non.reporting%>%pull(Finyear)
  dis.rec.yr=subset(dis.rec.yr,dis.rec.yr!='NA')
  relyrs=range(unique(d$Yr.rel),na.rm=T)     #range of released years
  recapyrs=range(dis.rec.yr,na.rm=T)                        #range of recaptured years
  N=c(table(d$Yr.rel))                      #released numbers by year
  recapharv=d%>%   #matrix of number of recoveries (row= release year; column= recapture year; lower triangle filled with -1)
    filter(Yr.rec%in%dis.rec.yr)%>%
    group_by(Yr.rel,Yr.rec)
  recapharv=table(recapharv$Yr.rel,factor(recapharv$Yr.rec,levels=recapyrs[1]:recapyrs[2]))
  RwN=rownames(recapharv)
  ClN=colnames(recapharv)
  recapharv=matrix(recapharv,nrow=nrow(recapharv),ncol=ncol(recapharv))
  rownames(recapharv)=RwN
  colnames(recapharv)=ClN
  if(!nrow(recapharv)==length(N) | !ncol(recapharv)==length(dis.rec.yr) )
  {
    add.row=which(!names(N)%in%RwN)
    if(length(add.row)>0)
    {
      ddumy=matrix(rep(0,ncol(recapharv)),nrow=length(add.row),ncol=ncol(recapharv))
      rownames(ddumy)=names(N)[add.row]
      recapharv=rbind(recapharv,ddumy)
      recapharv=recapharv[order(rownames(recapharv)),]
    }
    
    add.col=which(!dis.rec.yr%in%ClN)
    if(length(add.col)>0)
    {
      ddumy=matrix(rep(0,nrow(recapharv)),ncol=length(add.col),nrow=nrow(recapharv))
      colnames(ddumy)=dis.rec.yr[add.col]
      recapharv=cbind(recapharv,ddumy)
      recapharv=recapharv[,order(colnames(recapharv))]
    }
  }

  recapharv[lower.tri(recapharv)]=-1

  input.list=list(relyrs=relyrs,
                  recapyrs=recapyrs,
                  N=N,
                  recapharv=recapharv,
                  lambda=1-Non.reporting$Mean.non.rep,
                  phi=rep(1-tag.shed,length(dis.rec.yr)),
                  Fyr=dis.rec.yr,
                  Myr=dis.rec.yr[1])
  
  #Estimate mortality     
  model<-with(input.list,irm_h_fixedM(relyrs = relyrs, recapyrs = recapyrs, 
                                      N = N, recapharv = recapharv,lambda = lambda,
                                      phi = phi, Fyr = Fyr, Myr = Myr, initial = c(0.1,M.fixed), 
                                      lower = c(0.0001,0.0001),upper = c(5,5), maxiter = 10000,M.estim= FALSE))
  
  with(model$Fishing.mortality,plot(Yr,F,main=MAIN,pch=19,col=3,ylim=c(0,max(F))))
  #with(model$fishing_mortality,plot(Year,F,main=MAIN)
  #with(model$natural_mortality,plot(Year,M))
  return(model)
}
list.M=list(M.overall_minM=M.overall_minM,
            M.at.sel.peak_minM=M.at.sel.peak_minM,
            M.overall_meanM=M.overall_meanM,
            M.at.sel.peak_meanM=M.at.sel.peak_meanM)
Mortality.estimation=List.tagging   
for(s in 1:length(List.tagging))
{
  dd=list.M
  par(mfcol=c(2,2))
  for(i in 1:length(dd))
  {
   dd[[i]]=fn.fit.Z.estim(d=List.tagging[[s]],
                           M.fixed=list.M[[i]][match(unique(List.tagging[[s]]$Species),names(list.M[[i]]))],
                           Fishery='TDGDLF',
                          MAIN=names(List.tagging)[s])
    legend('top',paste(names(list.M)[i],round(list.M[[i]][s],3),sep='='),bty='n')
  }
  Mortality.estimation[[s]]=dd
}

Do.NSF=FALSE  #model not converging
if(Do.NSF)
{
  Mortality.estimation.NSF=List.tagging   
  for(s in 1:length(List.tagging.NSF))
  {
    dd=list.M
    par(mfcol=c(2,2))
    for(i in 1:length(dd))
    {
      dd[[i]]=fn.fit.Z.estim(d=List.tagging.NSF[[s]],
                             M.fixed=list.M[[i]][match(unique(List.tagging.NSF[[s]]$Species),names(list.M[[i]]))],
                             Fishery='NSF',
                             MAIN=names(List.tagging.NSF)[s])
      legend('top',paste(names(list.M)[i],round(list.M[[i]][s],3),sep='='),bty='n')
    }
    Mortality.estimation.NSF[[s]]=dd
  }
}

Do.Acous=FALSE  #model not converging
if(Do.Acous)
{
  Mortality.estimation.acoustic=List.tagging.acoustic   
  for(s in 1:length(List.tagging.acoustic))
  {
    dd=list.M
    par(mfcol=c(2,2))
    for(i in 1:length(dd))
    {
      dd[[i]]=fn.fit.Z.estim(d=List.tagging.acoustic[[s]]%>%filter(!(Recaptured=='Yes'&is.na(Yr.rec))),
                             M.fixed=list.M[[i]][match(unique(List.tagging.acoustic[[s]]$Species),names(list.M[[i]]))],
                             Fishery='TDGDLF - acoustic',
                             MAIN=names(List.tagging.acoustic)[s])
      legend('top',paste(names(list.M)[i],round(list.M[[i]][s],3),sep='='),bty='n')
    }
    Mortality.estimation.acoustic[[s]]=dd
  }
}


#---Export F-------
#note: use Published estimated by McAuley et al 2007 for dusky and sandbar and
#       the calculated estimate for whiskery and gummy based on Hoenig et al
hnd.indx=handl_OneDrive("Analyses/Data_outs/")

# Published by McAuley et al 2007
for(s in 1:length(Annual.Reported.F.McAuley.2007))
{
  Nm=capitalize(Species.names[match(names(Annual.Reported.F.McAuley.2007)[s],names(Species.names))])
  
  write.csv( Annual.Reported.F.McAuley.2007[[s]]%>%mutate(CV=0.2),
             paste(hnd.indx,Nm,'/',Nm,".Fishing.mortality.TDGDLF_McAuley2007.csv",sep=""),row.names=F)
  
}

# Hoenig et al 1998 (conventional (1994:1995) and acoustic (2012:2017)
out.hoenig=c("whiskery shark", "gummy shark")
this.out=match(out.hoenig,names(Mortality.estimation))
for(s in this.out)
{
  Nm=capitalize(names(Mortality.estimation)[s])
  dis.iers=used.released.years%>%
    filter(Species==names(Species.names[match(tolower(Nm),Species.names)]))%>%
    pull(Finyear)
  OUT=Mortality.estimation[[s]]$M.at.sel.peak_meanM$Fishing.mortality%>%
            filter(Yr%in%dis.iers)%>%
            dplyr::select(Yr,F,CV)%>%
            dplyr::rename(Finyear=Yr,
                          Mean=F)
  OUT.acoust=Mortality.estimation.acoustic[[match(names(Mortality.estimation)[s],names(Mortality.estimation.acoustic))]]$M.at.sel.peak_meanM$Fishing.mortality%>%
          filter(!is.na(CV) & Yr%in%these.acous.yrs)%>%
          dplyr::select(Yr,F,CV)%>%
          dplyr::rename(Finyear=Yr,
                        Mean=F)
  write.csv(rbind(OUT,OUT.acoust),paste(hnd.indx,Nm,'/',Nm,".Fishing.mortality.TDGDLF_Hoenig.csv",sep=""),row.names=F)
}

# Harry et al 2016 approach  
# not exported as it yielded similar results as Hoenig