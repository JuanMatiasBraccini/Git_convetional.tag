############## CONVENTIONAL TAGGING DATA ANALYSIS ##############

#notes: this script maps the location of sharks tagged with conventional tags and performs a range of different analysis
#       Speed values are low because sharks are not moving in straight line so this is minimum speed

rm(list=ls(all=TRUE))   #clear log


#Control Section
              #MISSING: update all figures and tables!
# ONLINE="NO"  #if working in Argentina  
ONLINE="YES"  #Switch to YES to update Changes done by dani



setwd("C:/Matias/Analyses/Conventional tagging")  #function for settting the working directory

source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Plot.Map.R")   #function for sourcing other .R scripts

library(RODBC)  			#library for importing excel data
library(gridBase)			#for inset map
library(PBSmapping)   #for polygon
library(MASS)  			#for kernel densities
library(plotrix)  			#needed for graph legends
library(geosphere)      #for spatial statistics and trigonometry in a sphere
library(circular)     #for circular stats
library(SDMTools)   #for creating lat and long grid
library(effects)    #for plotting all effects of GLM 
library(WriteXLS)   #for creating output sheets of excel file
library(climatol)   #for rose diagram
library(bootstrap)  #for jackknife
library(lubridate)
library(modeest)    #mode
library(dplyr)
library(sjPlot)  #word tables
library(tidyverse)
library(ggplot2)
library(NISTunits)
library(fields)
library(ggridges)
library(mgcv)

#memory.limit(3900)   #set memory limit to maximum size


###### DATA SECTION ############

      #SHAPE FILE PERTH ISLANDS
PerthIs=read.table("C:/Matias/Data/Mapping/WAislandsPointsNew.txt", header=T) #function for reading txt file
Rottnest.Is=subset(PerthIs,ID%in%c("ROTT1"))
Garden.Is=subset(PerthIs,ID%in%c("ROTT3"))


    #1.2. tagging data
if(ONLINE=="NO")
{
  Tagging=read.csv("C:/Matias/Data/Tagging/Conventional_tagging/Tagging_data.csv") 
  Species.Codes=read.csv("C:/Matias/Data/Species.code.csv")    
}
if(ONLINE=="YES") source("C:/Matias/Analyses/SOURCE_SCRIPTS/Git_other/Source_conventional_data.R")


    #1.3. bathymetry data
#    bathymetry data downloaded from http://topex.ucsd.edu/cgi-bin/get_data.cgi (Topography option)
Bathymetry_120=read.table("C:/Matias/Data/Mapping/get_data112_120.cgi")
Bathymetry_138=read.table("C:/Matias/Data/Mapping/get_data120.05_138.cgi")
Bathymetry=rbind(Bathymetry_120,Bathymetry_138)





###### PARAMETERS SECTION ############
size.mat=c(130,245,100,115)   #size at maturity of TK, BW, GM, WH  (FL, cm)    
size.birth=c(42.5,75,25,25)   #size at birth of TK, BW, GM, WH  (FL, cm)
names(size.mat)=names(size.birth)=c("Sandbar shark","Dusky shark","Gummy shark","Whiskery shark")
daysAtLiberty.threshold=30  #minimum time at liberty for population dynamics
min.daysAtLiberty=2         #minimum time at liberty for descriptive movement paper

###### MANIPULATE DATA ############

#Use only sharks tagged in WA
Tagging=subset(Tagging,Long.rels<=129 | is.na(Long.rels))


#3. CREATE VARIABLES

SPECIES=c("TK","BW","GM","WH")
Species.names=c("Sandbar shark","Dusky shark","Gummy shark","Whiskery shark")
names(Species.names)=SPECIES
names(SPECIES)=Species.names

# All.Species=as.character(Species.Codes$Species[c(1:40,99:102,104:124)])
# names(All.Species)=as.character(Species.Codes$COMMON_NAME[c(1:40,99:102,104:124)])
All.Species=as.character(Species.Codes$Species)
names(All.Species)=as.character(Species.Codes$COMMON_NAME)


#Keep only relevant species
Tagging=subset(Tagging,Species%in%All.Species)

Vec.all.species=unique(Tagging$Species)

#Remove missidentified gummies
Tagging=subset(Tagging,!(Species=="GM" & Lat.rels>(-26)))


#Valid Recapture (i.e. after or equal release date)      
Tagging$Jday.rel=with(Tagging,Yr.rel+(Mn.rel/12) + (Day.rel/(12*31)))
Tagging$Jday.rec=with(Tagging,Yr.rec+(Mn.rec/12) + (Day.rec/(12*31)))
Tagging$Valid.Rec=with(Tagging,ifelse(Recaptured=="Yes" & Jday.rec>=Jday.rel,"Yes",NA))

#Validate lengths
Tagging$Rel_FL=with(Tagging,
          ifelse(Species=="TK" & Rel_FL<37,NA,
          ifelse(Species=="BW" & Rel_FL<55,NA,
          ifelse(Species=="WH" & Rel_FL<25,NA,
          ifelse(Species=="GM" & Rel_FL<25,NA,
          Rel_FL)))))

Tagging$CAP_FL=with(Tagging,
          ifelse(Species=="TK" & CAP_FL<37,NA,
          ifelse(Species=="BW" & CAP_FL<55,NA,
          ifelse(Species=="WH" & CAP_FL<25,NA,
          ifelse(Species=="GM" & CAP_FL<25,NA,
          CAP_FL)))))

#Convert all factors to character 
Tagging %>% mutate(across(where(is.factor), as.character)) -> Tagging


#Add size bins
Tagging=Tagging%>%
  mutate(Rel_FL.bin=as.character(10*floor(Rel_FL/10)),
         CAP_FL.bin=as.character(10*floor(CAP_FL/10)))

#Time at liberty
Tagging$DATE_REL=as.Date( paste( Tagging$Yr.rel,Tagging$Mn.rel , Tagging$Day.rel , sep = "-" )  , format = "%Y-%m-%d" )
Tagging$DATE_CAPTR=as.Date( paste( Tagging$Yr.rec,Tagging$Mn.rec , Tagging$Day.rec , sep = "-" )  , format = "%Y-%m-%d" )
Tagging$DaysAtLarge=as.numeric(round((Tagging$DATE_CAPTR-Tagging$DATE_REL),0))
Tagging$DaysAtLarge=with(Tagging,ifelse(DaysAtLarge<0,NA,DaysAtLarge))

Tagging=Tagging%>%
  mutate(YrsAtLarge.bin=round(DaysAtLarge/365))

#Block
Tagging=Tagging%>%
  mutate(Block=paste(trunc(abs(Lat.rels)),trunc(Long.rels)-100,sep=''))


# Table A1. -----------------------------------------------------------------------
# Summary of numbers, size and sex ratios for all shark and ray species  
Tab1.fun=function(Spec)
{
  dat=subset(Tagging,Species==Spec)
  
    #Years rel and rec
  Yrs.rel.min=min(dat$Yr.rel,na.rm=T)
  Yrs.rel.max=max(dat$Yr.rel,na.rm=T)
  Yrs.rec.min=Yrs.rec.max=NA
  if(sum(dat$Yr.rec,na.rm=T)>0)
  {
    Yrs.rec.min=min(dat$Yr.rec,na.rm=T)
    Yrs.rec.max=max(dat$Yr.rec,na.rm=T)
  }
    
  
    #Number rel and rec                
  N.rel=length(dat$Species)  
  N.rec=NA
  dat.rec=subset(dat,Recaptured=="Yes")
  if(nrow(dat.rec)>0) N.rec=nrow(dat.rec)
  
    #Recapture ratio
  Rec.ratio=NA
  if(!is.na(N.rec)) Rec.ratio=round(100*N.rec/N.rel,1)
  
    #Size rel and rec
  FL.rel.min=min(dat$Rel_FL,na.rm=T) 
  FL.rel.max=max(dat$Rel_FL,na.rm=T) 
  FL.rec.min=FL.rec.max=NA
  if(nrow(dat.rec)>0)
  {
       FL.rec.min=min(dat$CAP_FL,na.rm=T)
      FL.rec.max=max(dat$CAP_FL,na.rm=T)
  }
    
  
    #sex ratio rel and rec
  Sex.rel=table(dat$Sex)
  id.M=match("M",names(Sex.rel))
  id.F=match("F",names(Sex.rel))
  Rel.M=Sex.rel[id.M]
  Rel.F=Sex.rel[id.F]
  
  
  Rec.M=Rec.F=NA
  if(nrow(dat.rec)>0)
  {
    Sex.rel=table(dat.rec$Sex)
    id.M=match("M",names(Sex.rel))
    id.F=match("F",names(Sex.rel))
    Rec.M=Sex.rel[id.M]
    Rec.F=Sex.rel[id.F]
    
  }
  TL.measured=as.character(ifelse(Spec%in%TL.species,"TL","FL"))
  return(list(Yrs.rel.min=Yrs.rel.min,Yrs.rel.max=Yrs.rel.max,Yrs.rec.min=Yrs.rec.min,
              Yrs.rec.max=Yrs.rec.max,
              N.rel=N.rel,N.rec=N.rec,Rec.ratio=Rec.ratio,FL.rel.min=FL.rel.min,FL.rel.max=FL.rel.max,FL.rec.min=FL.rec.min,
              FL.rec.max=FL.rec.max,Rel.M=Rel.M,Rel.F=Rel.F,Rec.M=Rec.M,Rec.F=Rec.F,
              TL.measured=TL.measured))
  
}

TL.species=c('ZE','PC','FR','TN','PZ','PM','SR')     
Tabl1.list=vector('list',length(Vec.all.species))
names(Tabl1.list)=Vec.all.species
for (i in 1:length(Vec.all.species)) Tabl1.list[[i]]=Tab1.fun(Spec=Vec.all.species[i])

Tabl1.matrix=matrix(nrow=length(Vec.all.species),ncol=length(Tabl1.list[[1]]))

for (i in 1:length(Vec.all.species))Tabl1.matrix[i,]=do.call(cbind,Tabl1.list[[i]])
colnames(Tabl1.matrix)=names(Tabl1.list[[1]])
rownames(Tabl1.matrix)=Vec.all.species

Tabl1.matrix=as.data.frame(Tabl1.matrix)
Tabl1.matrix$Species=rownames(Tabl1.matrix)
Tabl1.matrix=merge(Tabl1.matrix,Species.Codes[,1:2],by="Species")
Tabl1.matrix$N.rel=as.numeric(as.character(Tabl1.matrix$N.rel))
Tabl1.matrix=Tabl1.matrix[order(-Tabl1.matrix$N.rel),]

setwd("C:/Matias/Analyses/Conventional tagging/General movement/outputs")
write.csv(Tabl1.matrix,"Paper/Table.A1.csv",row.names=F)

#word table
Tabl1.doc=Tabl1.matrix %>%
          mutate(across(where(is.factor), as.character),
                 N.rel=as.numeric(N.rel),
                 Rel.M=as.numeric(Rel.M),
                 Rel.F=as.numeric(Rel.F),
                 N.rec=as.numeric(N.rec),
                 Rec.M=as.numeric(Rec.M),
                 Rec.F=as.numeric(Rec.F),
                 FL.rel.min=as.numeric(FL.rel.min),
                 FL.rel.max=as.numeric(FL.rel.max),
                 FL.rec.min=as.numeric(FL.rec.min),
                 FL.rec.max=as.numeric(FL.rec.max))

Tabl1.doc=left_join(Tabl1.doc,Species.Codes%>%
                         dplyr::select(Species,SCIENTIFIC_NAME),
                       by='Species')%>%
              mutate(SP=paste(COMMON_NAME," (",SCIENTIFIC_NAME,")",sep=''),
                     Rel.yr=paste(Yrs.rel.min,Yrs.rel.max,sep='-'),
                     Rec.yr=paste(Yrs.rec.min,Yrs.rec.max,sep='-'),
                     N.rel.u=N.rel-Rel.M-Rel.F,
                     N.rec.u=N.rec-Rec.M-Rec.F,
                     Rel.fl=paste(round(FL.rel.min),round(FL.rel.max),sep='-'),
                     Rec.fl=paste(round(FL.rec.min),round(FL.rec.max),sep='-'),
                     SP=ifelse(TL.measured=="TL",paste(SP,'*',sep=''),SP))%>%
              dplyr::select(SP,Rel.yr,Rec.yr,Rel.M,Rel.F,N.rel.u,
                            Rec.M,Rec.F,N.rec.u,Rec.ratio,Rel.fl,Rec.fl)%>%
              rename(Species=SP,
                     Release.year=Rel.yr,
                     Recapture.year=Rec.yr,
                     Release.numbers.Male=Rel.M,
                     Release.numbers.Female=Rel.F,
                     Release.numbers.Unknown=N.rel.u,
                     Recapture.numbers.Male=Rec.M,
                     Recapture.numbers.Female=Rec.F,
                     Recapture.numbers.Unknown=N.rec.u,
                     Recapture.ratio=Rec.ratio)%>%
            replace(is.na(.),'')%>%
            mutate(Recapture.year=ifelse(Recapture.year=='NA-NA','',Recapture.year),
                   Rel.fl=ifelse(Rel.fl=='NA-NA','',Rel.fl),
                   Rec.fl=ifelse(Rec.fl=='NA-NA','',Rec.fl))%>%
            rename('Size.range(cm).Release'=Rel.fl,
                   'Size.range(cm).Recapture'=Rec.fl)
tab_df(Tabl1.doc,file="Paper/Table.A1.doc")
 


# Early exports -----------------------------------------------------------

#export data for mark recapture analysis
write.csv(Tagging,"C:/Matias/Analyses/Data_outs/Tagging_conventional.data.csv",row.names=F)

 
#export data for Sarah Jakobs analysis
do.Sarah=FALSE
if(do.Sarah)
{
  Sarah.J=subset(Tagging,Species=="GN")
  Sarah.J= Sarah.J %>% select (-c(Areas, Areas.rec,Taxa,CAES_Code,Effort,Reporting,Valid.Rec))
  write.csv(Sarah.J,"C:/Matias/Students/Sarah Jakobs/Data/Conv_tag.csv",row.names=F)
}

#export data for Colin's mpa analysis
do.Colin=FALSE
if(do.Colin)
{
  dummy=subset(Tagging,Recaptured=="Yes")
  dummy=table(dummy$Species)
  dummy=dummy[dummy>=10]
  Colin.D=subset(Tagging,Species%in%names(dummy) & Recaptured=="Yes",
                 select=c(Tag.no,Species,SCIENTIFIC_NAME,Sex,Rel_FL,
                          Lat.rels,Long.rels,Lat.rec,Long.rec,
                          Day.rel,Mn.rel,Yr.rel,Day.rec,Mn.rec,Yr.rec))
  names(Colin.D)[match("Tag.no",names(Colin.D))]="TagID"
  write.csv(Colin.D,"C:/Matias/Data/Tagging/Conventional_tagging/data_for_Colin/Colin_conv_tag.csv",row.names=F)
  
}

# Chi-square test of sex ratios -----------------------------------------------------------

X2.fun=function(Spec)
{
  dat=subset(Tagging,Species==Spec)
  dat.rec=subset(dat,Recaptured=="Yes")
  
  X2.rel=X2.rec=NA
  
  #sex ratio rel and rec
  Sex.rel=table(dat$Sex)
  id=match(c("M","F"),names(Sex.rel))
  Sex.rel=Sex.rel[id]
  if(sum(is.na(Sex.rel))>=1)Sex.rel=NULL
  if(length(Sex.rel)>1)X2.rel=chisq.test(Sex.rel)

  Sex.rec=NULL
  Tab.out=NULL
  if(nrow(dat.rec)>0)
  {
    Sex.rec=table(dat.rec$Sex)
    id=match(c("M","F"),names(Sex.rec))
    Sex.rec=Sex.rec[id]
    if(sum(is.na(Sex.rec))>=1)Sex.rec=NULL
    if(length(Sex.rec)>1)X2.rec=chisq.test(Sex.rec)
    
    if(Spec%in%SPECIES)Tab.out=data.frame(Species=dat$COMMON_NAME[1],
                       M.rel=Sex.rel[1],
                       F.rel=Sex.rel[2],
                       X2.p=round(X2.rel$p.value,2),
                       M.rec=Sex.rec[1],
                       F.rec=Sex.rec[2],
                       X2.rec.p=round(X2.rec$p.value,2))
  }
  return(list(n.rel=Sex.rel,X2.rel=X2.rel,n.rec=Sex.rec,X2.rec=X2.rec,Tab.out=Tab.out))
  
}
X2.list=vector('list',length(Vec.all.species))
names(X2.list)=Vec.all.species
for (i in 1:length(Vec.all.species)) X2.list[[i]]=X2.fun(Spec=Vec.all.species[i])

Tab.Chi.sex=do.call(rbind,sapply(X2.list,function(x)x[5]))
rownames(Tab.Chi.sex)=NULL
tab_df(Tab.Chi.sex,file="Paper/Tabl.Chi.sex.doc")

# Recapture methods  -----------------------------------------------------------
Rec.m=table(Tagging$COMMON_NAME,Tagging$CAPT_METHD)
write.csv(Rec.m,"Recapture.methods.csv")

# MAPPING  -----------------------------------------------------------
add.depth="NO"

Tagging$Lat.rec=with(Tagging,ifelse(Lat.rec==0,NA,Lat.rec))
Tagging$Long.rec=with(Tagging,ifelse(Long.rec==0,NA,Long.rec))

#Create corners for moving around land
Exmouth=cbind(113.843,-21.81416)
Shark.bay=cbind(112.921,-25.497)
Mid.point=cbind(116.425,-35.043)
Streaky.Bay=cbind(134.2052,-32.71551)
Spencer=cbind(135.6149,-35.06194)
North.West.Cape=c(114.1333,-21.8825)
Cape.Leuwin=c(115.17,-34.35)
Albany=c(117.8847,-34.7)

#bathymetry
if(add.depth=="YES")
{
  Bathymetry=Bathymetry[order(Bathymetry$V1,Bathymetry$V2),]
  xbat=sort(unique(Bathymetry$V1))
  ybat=sort(unique(Bathymetry$V2))
  reshaped=as.matrix(reshape(Bathymetry,idvar="V1", timevar="V2",v.names="V3", direction="wide"))
}

range.long=c(111,138)
range.lat=c(-37,-11)

data(worldLLhigh)

#Function just showing releases and recaptures
Map.fn=function(Spec,LABEL)
{
  dat=subset(Tagging,Species%in%Spec)
  
  plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="grey85", axes=F, xlab="", ylab="",
          border="black",bg="grey99",plt = NULL)
  points(dat$Long.rel,dat$Lat.rel,pch=21,col="grey30",cex=1.1,bg="grey60")
  points(dat$Long.rec,dat$Lat.rec,pch=21,col="black",cex=1.25,bg="white")
  
  
  #add axis
  axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = F,
       tck=-0.035,las=1,cex.axis=1.2)
  axis(side = 1, at = seq(range.long[1],range.long[2],1), labels = F,tck=-0.02) 
  axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = F,
       tck=-0.035,las=1,cex.axis=1.2)
  axis(side = 2, at = seq(range.lat[1],range.lat[2],1), labels = F,tck=-0.02)
  box()
  contour(xbat, ybat, z=reshaped[,2:ncol(reshaped)],xlim=WOz.long, ylim=WOz.lat, zlim=c(-1,-300),nlevels = 3,
          labcex=0.75,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
  legend("topright",LABEL,bty='n',cex=1.9)
}

#Function showing differences in sizes and release locations

SPECs=list("BW","GM","WH","TK")
FL.rng=list(c(100,200,300),c(50,100,150),c(50,100,150),c(50,150,250))
names(FL.rng)=unlist(SPECs)
ArEaS=list(c("closed","JANSF","WANCSF"),"WCDGDLF",c("JASDGDLF.zone1","JASDGDLF.zone2"))


LABELS=list("Dusky shark","Gummy shark","Whiskery shark","Sandbar shark")

# SPECs=list("BW","GM","WH","TK",c("LG","WW","CP","TG","NS","HZ","PN","MS","GR","GN","PZ","WB","WK"))
# LABELS=list("Dusky shark","Gummy shark","Whiskery shark","Sandbar shark","Others")
fn.inset1=function(XE,WHERE,PCHs,PiCCol,CiX,BGcol)
{
  
  ColsS="black";Y1=108
  
  line.fn=function(LONG,LAT)lines(LONG,LAT,col=ColsS,lwd=1.5)
  
  par(fig=WHERE, new = T,mgp=c(.1,.4,0),mai=c(.01,.01,.01,.01))
  plotMap(worldLLhigh, xlim=c(Y1,range.long[2]), ylim=c(-36.5,-13),col="black", axes=F, xlab="", ylab="",
          border="black",bg="white",plt = NULL)
  line.fn(c(116.5,116.5),c(-36.5,-35))
  text(123,-34.5,"JASDGDLF",col=ColsS,cex=XE,font=2)
  text(123,-35.5,"(zone 2)",col=ColsS,cex=XE,font=2)
  text(112,-34.5,"JASDGDLF",col=ColsS,cex=XE,font=2)
  text(112,-35.5,"(zone 1)",col=ColsS,cex=XE,font=2)
  
  line.fn(c(Y1,115.6),c(-33,-33))
  text(111.75,-30,"WCDGDLF",col=ColsS,cex=XE,font=2)
  
  line.fn(c(Y1,113.45),c(-26.5,-26.5))
  text(111,-23,"Closed",col=ColsS,cex=XE,font=2)
  text(111,-24,"area",col=ColsS,cex=XE,font=2)
  
  line.fn(c(114,114),c(-21.75,-17))
  text(118,-19,"WANCSF",col=ColsS,cex=XE,font=2)
  
  line.fn(c(123.75,123.75),c(-16.12,-13))
  text(126,-13.5,"JANSF",col=ColsS,cex=XE,font=2)
  box(lwd=2)
  
  text(122,-24,"Western",col="white",cex=1,font=2)
  text(122,-27,"Australia",col="white",cex=1,font=2)
  
  line.fn(c(129,129),c(-36.5,-31.7))
  text(133,-27,"South",col="white",cex=0.8,font=2)
  text(133,-28.75,"Australia",col="white",cex=0.8,font=2)
  lines(c(129,129),c(-31.7,-14.99),col="white",lwd=1.5,lty=3)
  
  #turning points
  points(North.West.Cape[1],North.West.Cape[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  points(Shark.bay[1],Shark.bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  points(Cape.Leuwin[1],Cape.Leuwin[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  points(Albany[1],Albany[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  points(Streaky.Bay[1],Streaky.Bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  points(Spencer[1],Spencer[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  
  
}
#fn.inset1(XE=0.8,WHERE=c(.375,.675,.85,.99),PCHs=21,PiCCol=1,CiX=1.25,BGcol="grey95")

fn.inset=function(WHERE)
{
  par(fig=WHERE, new = T,mgp=c(.1,.4,0),mai=c(.01,.01,.01,.01))
  plotMap(worldLLhigh, xlim=c(108,range.long[2]), ylim=c(-36.5,-13),
          col="grey20",bg="white", axes=F, xlab="", ylab="",
          border="black",plt = NULL)
  
  lines(c(129,129),c(-31.7,-14.99),col="grey95",lwd=1.5,lty=3) 
  lines(c(129,129),c(-36.5,-31.7),lwd=1.5,lty=3)
  text(122,-24,"Western",col="white",cex=.825,font=2)
  text(122,-28,"Australia",col="white",cex=.825,font=2)
  
  #text(133,-27,"South",col="white",cex=0.8,font=2)
  #text(133,-28.75,"Australia",col="white",cex=0.8,font=2)
  box(lwd=1.5)
}

fn.axis=function()
{
  axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = F,
       tck=-0.035,las=1,cex.axis=1.2)
  axis(side = 1, at = seq(range.long[1],range.long[2],1), labels = F,tck=-0.02) 
  axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = F,
       tck=-0.035,las=1,cex.axis=1.2)
  axis(side = 2, at = seq(range.lat[1],range.lat[2],1), labels = F,tck=-0.02)
  box()
}

Tagging$Lat.rels=with(Tagging,ifelse(Species=="BW" & Long.rels<114 & Lat.rels>-17.7,-21.765,Lat.rels))

blk.wht=TRUE
if(blk.wht)
{
  Le.col=c("grey65","white","black")
  BordeR=1
}else
{
   Le.col=c(rgb(1,.1,.1,alpha=.3),
            rgb(.1,1,.1,alpha=.3),
            rgb(.1,.1,1,alpha=.3))
   BordeR='transparent'
}


Map2.fn=function(Spec,BKS,y)
{  
  dat=subset(Tagging,Species%in%Spec & !is.na(Rel_FL))
  
  dat$col=cut(dat$Rel_FL,BKS)
  CUTS=levels(dat$col)
   dat$col=with(dat,ifelse(col==CUTS[1],Le.col[1],
                    ifelse(col==CUTS[2],Le.col[2],
                    ifelse(col==CUTS[3],Le.col[3],
                    NA))))
  if(Spec=="TK")dat$Long.rels=with(dat,ifelse(Long.rels<114&Lat.rels>(-18),121,Long.rels))
    
    
  
  plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="grey97", axes=F, xlab="", ylab="",
          border="black",bg="white",plt = NULL)
  points(dat$Long.rels,dat$Lat.rels,pch=21,bg=dat$col,cex=1.2,col=BordeR)
  fn.axis()
  axis(side = 2, at = seq(range.lat[1],range.lat[2],4), 
       labels = -seq(range.lat[1],range.lat[2],4),tck=-0.035,las=1,cex.axis=1.1)
  if(add=="YES")   axis(side = 1, at = seq(range.long[1],range.long[2],4),
                        labels = seq(range.long[1],range.long[2],4),tck=-0.02,cex.axis=1.2) 
  legend(107,-10,y,bty='n',cex=1.05,xjust=0) 
  if(TIT=="YES") mtext("Released",3,line=0,cex=1.2)
  
  
  plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="grey97", axes=F, xlab="", ylab="",
          border="black",bg="white",plt = NULL)
  points(dat$Long.rec,dat$Lat.rec,pch=21,col=BordeR,cex=1.2,bg=dat$col)
  fn.axis()   
  if(add=="YES")   axis(side = 1, at = seq(range.long[1],range.long[2],4),
                        labels = seq(range.long[1],range.long[2],4),tck=-0.02,cex.axis=1.1) 
  if(TIT=="YES") mtext("Recaptured",3,line=0,cex=1.2)
}
fn.lg=function(x)
{
  legend('right',x,bty='n',pch=rep(21,3),
         pt.cex=1.5,pt.bg=Le.col,title="FL (cm)")
}
add="NO"
TIT="YES"
tiff(file="Paper/Figure1_Map.tiff",width = 1000, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfrow=c(4,2),mai=c(.05,.05,.05,.05),oma=c(3,3,1,.01),mgp=c(1,.5,0))
Map2.fn("TK",BKS=c(0,100,150,270),"Sandbar shark")
fn.lg(c("<100","100-150",">150"))
TIT="NO"
Map2.fn("BW",BKS=c(0,100,200,300),"Dusky shark")
fn.lg(c("<100","100-200",">200"))
Map2.fn("GM",BKS=c(0,90,120,200),"Gummy shark")
fn.lg(c("<90","90-120",">120"))
add="YES"
Map2.fn("WH",BKS=c(0,90,120,200),"Whiskery shark")
fn.lg(c("<90","90-120",">120"))
mtext(expression("Latitude ("*~degree*S*")"),side=2,outer=T,line=1.35,font=1,las=0,cex=1.2)
mtext(expression("Longitude ("*~degree*E*")"),side=1,outer=T,line=1.75,font=1,las=0,cex=1.2)
fn.inset(WHERE=c(.2,.475,.75,.99))
dev.off()


add.Map.simple="NO"
if(add.Map.simple=="YES")
{
  tiff(file="Paper/Map2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(2,2),mai=c(.15,.15,.05,.05),oma=c(3,3,.01,.01))
  #par(mfcol=c(2,2),mai=c(.5,.6,.05,.05),oma=c(.2,.2,.01,.01))
  for(i in c(4,1:3))
  {
    Map.fn(SPECs[[i]],LABELS[[i]])
    if(i %in%c(1,3))
    {
      axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = seq(range.long[1],range.long[2],2),
           tck=-0.035,las=1,cex.axis=1.2)
    }
    if(i %in%c(1,4))
    {
      axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = seq(range.lat[1],range.lat[2],2)*-1,
           tck=-0.035,las=1,cex.axis=1.2)
    }
    if(i==2)
    {
      PCHs=25
      CiX=2.5
      BGcol="black"
      PiCCol= "white"
      legend("right",c("released","recaptured"),bty='n',pch=21,col=c("grey30","black"),
             pt.bg=c("grey60","white"),cex=1.75,pt.cex=2)
      points(North.West.Cape[1],North.West.Cape[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
      points(Shark.bay[1],Shark.bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
      points(Cape.Leuwin[1],Cape.Leuwin[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
      points(Albany[1],Albany[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
      points(Streaky.Bay[1],Streaky.Bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
      points(Spencer[1],Spencer[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
    }
    
  }
  mtext("Latitude (?S)",side=2,outer=T,line=1.5,font=1,las=0,cex=1.5)
  mtext("Longitude (?E)",side=1,outer=T,line=1.25,font=1,las=0,cex=1.5)
  
  #Add inset
  XE=0.6;ColsS="black";Y1=108
  
  line.fn=function(LONG,LAT)lines(LONG,LAT,col=ColsS,lwd=1.5)
  
  par(fig=c(.0,.7,.685,.9), new = T,mgp=c(.1,.4,0),mai=c(.01,.01,.01,.01))
  plotMap(worldLLhigh, xlim=c(Y1,range.long[2]), ylim=c(-36.5,-13),col="black", axes=F, xlab="", ylab="",
          border="black",bg="white",plt = NULL)
  line.fn(c(116.5,116.5),c(-36.5,-35))
  text(123,-34.5,"JASDGDLF",col=ColsS,cex=XE,font=2)
  text(123,-35.5,"(zone 2)",col=ColsS,cex=XE,font=2)
  text(112,-34.5,"JASDGDLF",col=ColsS,cex=XE,font=2)
  text(112,-35.5,"(zone 1)",col=ColsS,cex=XE,font=2)
  
  line.fn(c(Y1,115.6),c(-33,-33))
  text(111.75,-30,"WCDGDLF",col=ColsS,cex=XE,font=2)
  
  line.fn(c(Y1,113.45),c(-26.5,-26.5))
  text(111,-23,"Closed",col=ColsS,cex=XE,font=2)
  text(111,-24,"area",col=ColsS,cex=XE,font=2)
  
  line.fn(c(114,114),c(-21.75,-17))
  text(118,-19,"WANCSF",col=ColsS,cex=XE,font=2)
  
  line.fn(c(123.75,123.75),c(-16.12,-13))
  text(126,-13.5,"JANSF",col=ColsS,cex=XE,font=2)
  box(lwd=2)
  
  text(122,-24,"Western",col="white",cex=1,font=2)
  text(122,-27,"Australia",col="white",cex=1,font=2)
  
  line.fn(c(129,129),c(-36.5,-31.7))
  text(132.5,-34,"South",col=ColsS,cex=XE*1.15,font=2)
  text(132.5,-35.3,"Australia",col=ColsS,cex=XE*1.15,font=2)
  lines(c(129,129),c(-31.7,-14.99),col="white",lwd=1.5,lty=3)
  
  dev.off()
  
}

#Map Release and recaptures for other species
Tabl1.matrix$N.rec=as.numeric(as.character(Tabl1.matrix$N.rec))
Others.rep.map=subset(Tabl1.matrix,N.rec>=5 & !Species%in%c("TK", "BW","GM", "WH"))$Species
add="NO"
TIT="YES"
tiff(file="Paper/Figure1_Map.Other.Species.tiff",width = 1000, height = 2000,units = "px", res = 300,compression = "lzw")
par(mfrow=c(4,2),mai=c(.05,.05,.05,.05),oma=c(3,3,1,.01),mgp=c(1,.5,0))
Map2.fn("TG",BKS=c(0,100,200,300),"Tiger shark")
fn.lg(c("<100","100-200",">200"))
TIT="NO"
Map2.fn("CP",BKS=c(0,100,150,300),"Bronze whaler")
fn.lg(c("<100","100-200",">200"))
Map2.fn("WW",BKS=c(0,90,120,200),"Western wobbegong")
fn.lg(c("<90","90-120",">120"))
add="YES"
Map2.fn("LG",BKS=c(0,90,120,200),"Spinner shark")
fn.lg(c("<90","90-120",">120"))
mtext(expression("Latitude ("*~degree*S*")"),side=2,outer=T,line=1.35,font=1,las=0,cex=1.2)
mtext(expression("Longitude ("*~degree*E*")"),side=1,outer=T,line=1.75,font=1,las=0,cex=1.2)
fn.inset(WHERE=c(.2,.475,.75,.99))
dev.off()


#Fig. A1. Zones and turning points
do.Fig.A1="NO"
if(do.Fig.A1=="YES")
{
  #Shark zones
  library(rgdal)
  JA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/JA_Northern_Shark.shp", layer="JA_Northern_Shark") 
  WA_Northern_Shark=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthCoastShark_s43.shp", layer="NorthCoastShark_s43") 
  WA_Northern_Shark_2=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/NorthWestCoastShark_s43.shp", layer="NorthWestCoastShark_s43") 
  SDGDLL_zone1=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/SDGDLL_zone1.shp", layer="SDGDLL_zone1") 
  SDGDLL_zone2=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/SDGDLL_zone2.shp", layer="SDGDLL_zone2") 
  WCDGDLL=readOGR("C:/Matias/Data/Mapping/Shark_shape_files/WCDGDLL.shp", layer="WCDGDLL") 
  
  Y=range.lat
  Y[1]=-39
  
  tiff(file="Paper/MapA1.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
  par(mfrow=c(1,1),mai=c(.05,.05,.05,.05),oma=c(3,1,1,.01),mgp=c(1,.75,0))
  
  plotMap(worldLLhigh, xlim=range.long, ylim=Y,col="grey90", axes=F, xlab="", ylab="",
          border="black",bg="white",plt = NULL)  
  
  #Add zones
  plot(WA_Northern_Shark,add=T,col="grey85")
  plot(JA_Northern_Shark,ylim=Y,xlim=range.long,add=T,col="grey60")
  plot(WA_Northern_Shark_2,add=T,col="grey70",angle=0,density=seq(5,35,2))
  plot(WA_Northern_Shark_2,add=T,col="grey70",angle=90,density=seq(5,35,2))
  plot(WCDGDLL,add=T,col="grey70")
  plot(SDGDLL_zone1,add=T,col="white")
  plot(SDGDLL_zone2,add=T,col="grey55")
  box()
  
  #turning points
  PCHs=21;PiCCol=1;CiX=2;BGcol="white"
  points(North.West.Cape[1],North.West.Cape[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(North.West.Cape[1]*1.04,North.West.Cape[2]*1.01,"North West Cape",cex=1,font=2)
  
  points(Shark.bay[1],Shark.bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(Shark.bay[1]*1.03,Shark.bay[2],"Shark Bay",cex=1,font=2)
  
  points(Cape.Leuwin[1],Cape.Leuwin[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(Cape.Leuwin[1]*1.025,Cape.Leuwin[2]*0.94,"Cape Leeuwin",cex=1,font=2,srt=45)
  
  points(Albany[1],Albany[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(Albany[1]*1.0125,Albany[2]*0.955,"Albany",cex=1,font=2,srt=45)
  
  points(Streaky.Bay[1],Streaky.Bay[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(Streaky.Bay[1]*1.005,Streaky.Bay[2]*0.91,"Streaky Bay",cex=1,font=2,srt=80)
  
  points(Spencer[1],Spencer[2],pch=PCHs,col=PiCCol,cex=CiX,bg=BGcol)
  text(Spencer[1]*1.005,Spencer[2]*0.91,"Port Lincoln",cex=1,font=2,srt=80)
  
  XE=1.25;ColsS=1
  text(117,-35.6,"JASDGDLF",col=ColsS,cex=XE,font=2)
  text(123,-36.5,"(zone 2)",col=ColsS,cex=XE-0.1,font=2)
  text(113.5,-36.5,"(zone 1)",col=ColsS,cex=XE-0.1,font=2,, srt=320)
  text(112.5,-30,"WCDGDLF",col=ColsS,cex=XE,font=2, srt=90)
  text(112,-22,"Closed",col=ColsS,cex=XE,font=2, srt=75)
  
  text(118.5,-17,"WANCSF",col=ColsS,cex=XE,font=2, srt=40)
  
  text(126.5,-13,"JANSF",col=ColsS,cex=XE,font=2)
  text(132,-34,"South",col=ColsS,cex=XE,font=2)
  text(132,-35.25,"Australia",col=ColsS,cex=XE,font=2)
  
  lines(rbind(129,129),rbind(-15,-31.5),lty=2,col="black",lwd=2)
  text(122,-24,"Western",col=1,cex=2,font=2)
  text(122,-28,"Australia",col=1,cex=2,font=2)
  
  axis(side = 2, at = seq(Y[1],Y[2],4), 
       labels = -seq(Y[1],Y[2],4),tck=-0.02,las=1,cex.axis=1.35)
  axis(side = 1, at = seq(range.long[1],range.long[2],4),
       labels = seq(range.long[1],range.long[2],4),tck=-0.02,cex.axis=1.35) 
  
  
  mtext("Latitude (?S)",side=2,outer=T,line=-1,font=1,las=0,cex=2)
  mtext("Longitude (?E)",side=1,outer=T,line=1.75,font=1,las=0,cex=2)
  dev.off()
  
  
}


# Size frequency of releases  -----------------------------------------------------------
#note: dodgy recapture size info, don't use, only report size at release

SizeFreq.rel.fn=function(SPEC,YARROW)
{
  datos=subset(Tagging,Species==SPEC)
  datos1=subset(datos,Recaptured=="Yes")
  
  datos1=subset(datos1,!is.na(Rel_FL))
  datos=subset(datos,!is.na(Rel_FL) & Recaptured=="No")
  
  Rango=range(c(datos1$Rel_FL,datos$Rel_FL),na.rm=T)
  Rango[1]=floor(Rango[1]/10)*10
  Rango[2]=ceiling(Rango[2]/10)*10
  Bks=seq(Rango[1],Rango[2],10)
  Rec.FL=hist(datos1$Rel_FL,breaks=Bks,plot=F)
  NotRec.FL=hist(datos$Rel_FL,breaks=Bks,plot=F)
  xvals=barplot(rbind(NotRec.FL$counts,Rec.FL$counts),beside=T,names.arg=Rec.FL$mids,
                col=c("grey80","grey40"),cex.axis=1.25,cex.names=1.25)
  axis(1,xvals[1,]+0.5,F)
  
  legend("topright",Species.names[i],bty="n",cex=1.75)
   
  id=which(abs(NotRec.FL$mids-size.mat[i])==min(abs(NotRec.FL$mids-size.mat[i])))
  id=id[length(id)]
  HERE=0
  arrows(xvals[1,id],YARROW,xvals[1,id],0,lwd=3,length = 0.15)
  #points(xvals[1,id],HERE,pch="*",cex=3)
  box()
  
  #KS test
  Kolmo=ks.test(datos$Rel_FL,datos1$Rel_FL)
  
  #Gear used for tagging and recapturing
  Table.GR.NoRec.rel=table(datos$REL_METHD)
  Table.GR.Rec.rel=table(datos1$REL_METHD)
  Table.GR.Rec.rec=table(datos1$CAPT_METHD)
  
  return(list(Kolmo=Kolmo,Table.GR.NoRec.rel=Table.GR.NoRec.rel,
        Table.GR.Rec.rel=Table.GR.Rec.rel,Table.GR.Rec.rec=Table.GR.Rec.rec))
}

KS.result=vector("list",length(SPECIES))
names(KS.result)=Species.names
yArrow=c(200,200,35,35)
tiff(file="Paper/Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
par(mfcol=c(2,2),oma=c(3,4,1.5,.1),mar=c(1,1,1,2),mgp=c(1, 0.75, 0),las=1)
for(i in 1:length(SPECIES)) 
{
  KS.result[[i]]=SizeFreq.rel.fn(SPEC=SPECIES[i],YARROW=yArrow[i])
  if(i==2) legend("right",c("Non-recaptured","Recaptured"),fill=c('gray80','gray40'),bty="n",cex=1.25)
}
mtext("Fork length (cm)",side=1,outer=T,line=1.5,font=1,las=0,cex=2)
mtext("Frequency",side=2,outer=T,line=2.1,font=1,las=0,cex=2)
dev.off()
Tab.Kolmo=as.data.frame(do.call(rbind,sapply(KS.result,function(x)x[1])))
Tab.Kolmo=cbind(Species=names(SPECIES),Tab.Kolmo)
Tab.Kolmo=data.frame(Species=unlist(Tab.Kolmo$Species),
                      statistic=unlist(Tab.Kolmo$statistic),
                      p.value=unlist(Tab.Kolmo$p.value))
rownames(Tab.Kolmo)=NULL
tab_df(Tab.Kolmo,file="Paper/Tabl.Kolmo.Smirnov.doc")


size.at.dat=data.frame(x=size.mat,y=.015,COMMON_NAME=names(size.mat),
                       x2=size.mat,y2=0)
Tagging%>%
  mutate(Rec=ifelse(Recaptured=="Yes","Recaptured","Not recaptured"))%>%
  filter(!is.na(Rel_FL) & Species%in%c('TK','GM','BW','WH'))%>%
  ggplot(aes(x=Rel_FL, fill=COMMON_NAME)) +
  geom_density(alpha=0.4) +
  geom_segment( aes(x = x, y = y, xend = x2, yend = y2, colour = COMMON_NAME), data = size.at.dat,
                size = 1.5,arrow = arrow(length = unit(0.2, "inches")), show.legend = F)+
  facet_wrap(~ Rec, ncol = 1)+xlab("Fork length (cm)") + 
  ylab("Density") + theme_bw() +labs(fill = "Species")+ theme(legend.position = c(0.8, 0.2))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.title=element_text(size=14),
        strip.text = element_text(size=14))
ggsave("Paper/Figure2_density.tiff", width = 6,height = 10, dpi = 300,compression = "lzw")


#release depth distribution
Tagging%>%
  filter(BOTDEPTH<400 & !is.na(Rel_FL) & Species%in%c('TK','GM','BW','WH'))%>%
  left_join(size.at.dat[,c("x","COMMON_NAME")],by="COMMON_NAME")%>%
  mutate(Mature=ifelse(Rel_FL>=x,"Mature","Immature"))%>%
  ggplot(aes(x=BOTDEPTH, fill=COMMON_NAME)) +
  geom_density(alpha=0.4) +
  facet_wrap(~ Mature, ncol = 1)+ xlab("Bottom depth (m)") + 
  ylab("Density") + theme_bw() +labs(fill = "Species")+ theme(legend.position = c(0.8, 0.2))+
  theme(axis.text=element_text(size=8),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        strip.text = element_text(size=10))
ggsave("Paper/FigureX_release.depth.dist.tiff", width = 4,height = 6, dpi = 300,compression = "lzw")

fun.depth.stats=function(SPEC)
{
  Dat=Tagging%>%
    filter(Species==SPEC)%>%
    filter(!is.na(BOTDEPTH))%>%
    filter(BOTDEPTH<400 & !is.na(Rel_FL))%>%
    left_join(size.at.dat[,c("x","COMMON_NAME")],by="COMMON_NAME")%>%
    mutate(Mature=ifelse(Rel_FL>=x,"Mature","Immature"))
    
  z=Dat%>%
    group_by(Mature)%>%
    summarise(Median=median(BOTDEPTH),
              Low.95=quantile(BOTDEPTH,probs=0.025),
              Up.95=quantile(BOTDEPTH,probs=0.975),
              Min=min(BOTDEPTH),
              Max=max(BOTDEPTH))
  return(cbind(Species=Dat$COMMON_NAME[1],z))
}
STORE.z=vector('list',length=length(SPECIES))
names(STORE.z)=SPECIES
for (i in 1:length(SPECIES)) STORE.z[[i]]=fun.depth.stats(SPEC=SPECIES[i])
Combined=do.call(rbind,STORE.z)
rownames(Combined)=NULL
tab_df(Combined,file="Paper/Tabl_depth.release.doc")



# Proportion of individuals in different release condition  -----------------------------------------------------------
Prop.fn=function(SPEC)
{
  datos.All=subset(Tagging,Species==SPEC)
  datos=subset(datos.All,Recaptured=="Yes")
  
  Cond=(table(datos$CONDITION))
  Cond.all=(table(datos.All$CONDITION))
  
  Cond=Cond[match(c("1","2","3"),names(Cond))]
  Cond.all=Cond.all[match(c("1","2","3"),names(Cond.all))]
  
  Condition=data.frame(Rec=Cond,Rel=Cond.all)
  
  #Chi-square distribution
  Xsq <- chisq.test(as.table(as.matrix(Condition)))
  return(Xsq)
}
Xsq.Cond=Xsq.SexRatio=KS.result
#for(i in 1:length(SPECIES)) Xsq.Cond[[i]]=Prop.fn(SPECIES[i])


# Sex ratios of release and recapture individuals  -----------------------------------------------------------
SexRatio.fn=function(SPEC)
{
  datos.All=subset(Tagging,Species==SPEC)
  datos=subset(datos.All,Recaptured=="Yes")
  
  Cond=(table(datos$Sex))
  Cond.all=(table(datos.All$Sex))
  
  Cond=Cond[match(c("F","M"),names(Cond))]
  Cond.all=Cond.all[match(c("F","M"),names(Cond.all))]
  Condition=data.frame(Rec=Cond,Rel=Cond.all)
  
  #Chi-square distribution
  Xsq <- chisq.test(as.table(as.matrix(Condition)))
  return(list(Xsq=Xsq,ratio=Condition))  
}
#for(i in 1:length(SPECIES)) Xsq.SexRatio[[Species.names[i]]]=SexRatio.fn(SPECIES[i])


# Create spatial cells  -----------------------------------------------------------

#create 1 degree cells
cell.size=1   #1 by 1 degree cells
rango.lat=round(range(c(Tagging$Lat.rels,Tagging$Lat.rec),na.rm=T),0)
rango.long=round(range(c(Tagging$Long.rels,Tagging$Long.rec),na.rm=T),0)
Total.lat.range=seq(rango.lat[2],rango.lat[1],by=-1)
Total.long.range=seq(rango.long[1],rango.long[2],by=1)
CELLS=expand.grid(Total.lat.range,Total.long.range)
colnames(CELLS)=c("Lat.NW","Long.NW")
CELLS$name=with(CELLS,(-Lat.NW*100)+(Long.NW-100)) #as per Simpfendorfer & McAuley unpublished
CELLS$centroid.Long=CELLS$Long.+(cell.size/2)
CELLS$centroid.Lat=CELLS$Lat.-(cell.size/2)

#create 10 minute cells

#select cells within 300m depth
# note: for running this code first run Arrows function below
# Arrows(SPECIES,"all species")
# grid(length(Total.long.range),length(Total.lat.range)-1, lwd = 2)
# text(CELLS$centroid.Long,CELLS$centroid.Lat,CELLS$name,cex=0.5)
shelf.cells=c(1422:1425,1521:1524,1621:1624,1719:1722,1818:1822,1915:1921,2015:2018,
              2114:2115,2213:2214,2312:2313,2412:2413,2512:2514,2612:2614,2712:2714,
              2813:2814,2914,3014:3015,3114:3115,3214:3215,3314:3315,3414:3416,3515:3518,
              3128:3131,3225:3233,3323:3334,3419:3438,3534:3538)
CELLS=subset(CELLS,name%in%shelf.cells)  #only keep CELLS occurrying in the shelf

#add cell borders
CELLS$West=CELLS$centroid.Long-(cell.size/2)    	
CELLS$East=CELLS$centroid.Long+(cell.size/2)	
CELLS$North=CELLS$centroid.Lat+(cell.size/2)
CELLS$South=CELLS$centroid.Lat-(cell.size/2)
CELLS=CELLS[order(CELLS$name),]


#add release and recapture BLOCK 
Tagging$Rel.name=-(ceiling(Tagging$Lat.rels)) * 100 +(floor(Tagging$Long.rels)-100)
Tagging$Rec.name=-(ceiling(Tagging$Lat.rec)) * 100 +(floor(Tagging$Long.rec)-100)


#Create tagging data for pop dyn model  -----------------------------------------------------------

Tagging=subset(Tagging,!is.na(Long.rels))
Tagging=subset(Tagging,!is.na(Lat.rels))

  #remove nonsense records                                                   
Tagging.pop.din=subset(Tagging, is.na(DaysAtLarge)| DaysAtLarge>=daysAtLiberty.threshold)

Pop.din.sp=c("BW","TK","GM","WH")
Tagging.pop.din=subset(Tagging.pop.din,Species%in%Pop.din.sp)

  #SA recaptured sharks reset to zone 2 following linear trajectory
SA.rels=subset(Tagging.pop.din,Long.rels>129)
Tagging.pop.din=subset(Tagging.pop.din,Long.rels<=129)  #exclude gummy releases in SA
Tagging.pop.din$SA.rec=with(Tagging.pop.din,ifelse(Long.rec>129,"Y","N"))

store.SA.recal=vector('list',length=3)   #no sandbars recaptured in SA so no need to recalculate
fn.interpolate.SA=function(SPeC)
{
  a=subset(Tagging.pop.din,SA.rec=="Y"& Species==SPeC)
  Tgs=unique(a$Tag.no)
  n=a$DaysAtLarge
  b=with(a,gcIntermediate(cbind(Long.rels,Lat.rels),cbind(Long.rec,Lat.rec),n=n, addStartEnd=F))
  for(i in 1:length(Tgs))
  {
    s=as.data.frame(b[[i]])
    names(s)=c("Long.rec","Lat.rec")
    s$DaysAtLarge=1:nrow(s)
    s$Tag.no=Tgs[i]
    s$DATE_CAPTR=subset(a,Tag.no==Tgs[i])$DATE_REL+(1:nrow(s))
    id=which(s$Long.rec<=129)
    s=s[id[length(id)],]
    b[[i]]=s
  }
  b=do.call(rbind,b)
  return(b)
}
store.SA.recal[[1]]=fn.interpolate.SA("GM")
store.SA.recal[[2]]=fn.interpolate.SA("WH")
store.SA.recal[[3]]=fn.interpolate.SA("BW")
store.SA.recal=do.call(rbind,store.SA.recal)
names(store.SA.recal)[c(1:3,5)]=paste(names(store.SA.recal)[c(1:3,5)],".1",sep="")

Tagging.pop.din=merge(Tagging.pop.din,store.SA.recal,by="Tag.no",all.x=T)
Tagging.pop.din$DATE_CAPTR=with(Tagging.pop.din,ifelse(SA.rec=='Y',DATE_CAPTR.1,DATE_CAPTR))
Tagging.pop.din$Lat.rec=with(Tagging.pop.din,ifelse(SA.rec=='Y',Lat.rec.1,Lat.rec))
Tagging.pop.din$Long.rec=with(Tagging.pop.din,ifelse(SA.rec=='Y',Long.rec.1,Long.rec))
Tagging.pop.din$DaysAtLarge=with(Tagging.pop.din,ifelse(SA.rec=='Y',DaysAtLarge.1,DaysAtLarge))


  #add zone released
Tagging.pop.din$Rel.zone=as.character(with(Tagging.pop.din,ifelse(Long.rels>=116.5 & Lat.rels<=(-26),"Zone2",
                         ifelse(Long.rels<116.5 & Lat.rels<=(-33),"Zone1",
                         ifelse(Lat.rels>(-33) & Lat.rels<=(-26) & Long.rels<116.5,"West",
                         ifelse(Lat.rels>(-26) & Long.rels<114,"Closed",
                         ifelse(Lat.rels>(-26) & Long.rels>=114 & Long.rels<123.75,"North",
                         ifelse(Lat.rels>(-26) & Long.rels>=123.75,"Joint",NA))))))))


  #add zone recaptured
Tagging.pop.din$Rec.zone=as.character(with(Tagging.pop.din,ifelse(Long.rec>=116.5 & Lat.rec<=(-26),"Zone2",
                      ifelse(Long.rec<116.5 & Lat.rec<=(-33),"Zone1",
                      ifelse(Lat.rec>(-33) & Lat.rec<=(-26) & Long.rec<116.5,"West",
                      ifelse(Lat.rec>(-26) & Long.rec<114,"Closed",
                      ifelse(Lat.rec>(-26) & Long.rec>=114 & Long.rec<123.75,"North",
                      ifelse(Lat.rec>(-26) & Long.rec>=123.75,"Joint",NA))))))))

  #Vector analysis to visualize movement differences among sexes and size for each species
SIZE=c(240,130,110,110)

for (i in 1:length(Pop.din.sp))
{
  Size=SIZE[i]
  a=subset(Tagging.pop.din,Species==Pop.din.sp[i])
  X=range(c(a$Long.rels,a$Long.rec),na.rm=T)
  Y=range(c(a$Lat.rels,a$Lat.rec),na.rm=T)
  if(i==3) Y=range(c(a$Lat.rec,min(a$Lat.rels)),na.rm=T)
  
  fem=subset(a,Sex=="F")
  mal=subset(a,Sex=="M")
  YR.rec=sort(unique(a$Yr.rec))
  
  tiff(file=paste("C:/Matias/Analyses/Conventional tagging/General movement/outputs/Maps.for.pop.dyn/",Pop.din.sp[i],".tiff",sep=""),
       width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  
  par(mfcol=c(2,2),las=1,cex.axis=1.25,cex.lab=1.5,cex.main=2)
  #plot(X,Y,col="transparent",main=paste("Female",Pop.din.sp[i],"all years"))
  plot(X,Y,col="transparent",main="Female",ylab="Lat",xlab="Long")
  arrows(fem$Long.rels,fem$Lat.rels,fem$Long.rec,fem$Lat.rec,col="pink",lwd=1.5)
  plot(X,Y,col="transparent",main="Male",ylab="Lat",xlab="Long")
  #plot(X,Y,col="transparent",main=paste("Male",Pop.din.sp[i],"all years"))
  arrows(mal$Long.rels,mal$Lat.rels,mal$Long.rec,mal$Lat.rec,col="blue",lwd=1.5)
  
  
  b=subset(a,Rel_FL<Size)
  plot(X,Y,col="transparent",main=paste("Size >",Size, "cm"),ylab="Lat",xlab="Long")
  arrows(b$Long.rels,b$Lat.rels,b$Long.rec,b$Lat.rec,col="red",lwd=1.5)
  
  b=subset(a,Rel_FL>Size)
  plot(X,Y,col="transparent",main=paste("Size >",Size, "cm"),ylab="Lat",xlab="Long")
  arrows(b$Long.rels,b$Lat.rels,b$Long.rec,b$Lat.rec,col="darkgreen",lwd=1.5)  
#   for(x in 1:length(YR.rec))
#   {
#     fem1=subset(fem,Yr.rec==YR.rec[x])
#     mal1=subset(mal,Yr.rec==YR.rec[x])
#     par(mfcol=c(2,1))
#     plot(X,Y,col="transparent",main=paste(Pop.din.sp[i],YR.rec[x]))
#     arrows(fem1$Long.rels,fem1$Lat.rels,fem1$Long.rec,fem1$Lat.rec,col="pink",lwd=1.5)
#     plot(X,Y,col="transparent",main=paste(Pop.din.sp[i],YR.rec[x]))
#     arrows(mal1$Long.rels,mal1$Lat.rels,mal1$Long.rec,mal1$Lat.rec,col="blue",lwd=1.5)
#     
#   }
  dev.off()
}

#a) Tag groups by age class
  #growth parameters (Linf is in FL)
Gr=list(BW=c(K.f=.0367,Linf.f=374.4,to.f=-3.3,K.m=0.045,Linf.m=337,to.m=-3),
        TK=c(K.f=.040,Linf.f=244.2,to.f=-4.8,K.m=0.044,Linf.m=226,to.m=-4),
        GM=c(K.f=0.123,Linf.f=(201.9-4.6424)/1.08,to.f=-1.55,K.m=0.253,Linf.m=(138.7-4.6424)/1.08,to.m=-0.9),
        WH=c(K.f=0.369,Linf.f=120.7,to.f=-0.6,K.m=.423,Linf.m=121.5,to.m=-0.472))
Mx.age=list(BW=45,TK=33,GM=16,WH=15)

fn.group=function(DAT,K.f,Linf.f,to.f,K.m,Linf.m,to.m,BySex,MX.AGE)
{
  names(DAT)[match(c("Rel.name","Rec.name"),names(DAT))]=c("Block.rel","Block.rec")
  DAT$Number=1
  DAT$Sex=with(DAT,ifelse(Sex=="U","F",Sex))
  
  #assign age
  if(unique(DAT$Species)=="GM")DAT$Rel_FL=DAT$Rel_FL*1.0837 +4.6424   #gummy growth pars are in TL

  
  #Calculate age from length
    #note: used the inverse of Von B because there is not enough info to create age-length key for any species
  DAT$Age=with(DAT,ifelse(Sex=="F",to.f-(1/K.f)*log(1-(Rel_FL/Linf.f)),
              ifelse(Sex=="M",to.m-(1/K.m)*log(1-(Rel_FL/Linf.m)),NA))) 
    
  #fix age for length > Linf 
  fem=subset(DAT,Sex=="F");mal=subset(DAT,Sex=="M")
  max.A.f=max(fem$Age,na.rm=T);max.A.m=max(subset(mal$Age,!mal$Age=='Inf'),na.rm=T)
  DAT$Age=with(DAT,ifelse(Sex=="F" & Rel_FL>=Linf.f,runif(1,max.A.f,MX.AGE),
                        ifelse(Sex=="M" & Rel_FL>=Linf.m,runif(1,max.A.m,MX.AGE),Age)))
  DAT$Age=with(DAT,ifelse(Age<0,0,Age))
               
  #Age recapture
  DAT$Age.rec=with(DAT,ifelse(!is.na(DATE_CAPTR),Age+(DaysAtLarge/365),NA))
  
  DAT$Age=round(DAT$Age)
  DAT$Age.rec=round(DAT$Age.rec)
  
  if(BySex=="YES")
  {
    #assign tag group
    DAT$TG.blk=as.factor(with(DAT,paste(Block.rel,Yr.rel,Mn.rel,Sex,Age)))
    DAT$TG.zn=as.factor(with(DAT,paste(Rel.zone,Yr.rel,Mn.rel,Sex,Age)))
    
    #recaptures
    DAT.rec=subset(DAT,Recaptured=="Yes")
    
    #All release numbers
    BLK.rel.f=aggregate(Number~TG.blk+Block.rel+Yr.rel+Mn.rel+Sex+Age+Rel_FL,subset(DAT,Sex=="F"),sum)
    BLK.rel.m=aggregate(Number~TG.blk+Block.rel+Yr.rel+Mn.rel+Sex+Age+Rel_FL,subset(DAT,Sex=="M"),sum)
    
    Zn.rel.f=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Sex+Age+Rel_FL,subset(DAT,Sex=="F"),sum)
    Zn.rel.m=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Sex+Age+Rel_FL,subset(DAT,Sex=="M"),sum)
    
    
    #Recapture numbers
    BLK.rec.f=aggregate(Number~TG.blk+Block.rec+Yr.rec+Mn.rec,subset(DAT.rec,Sex=="F"),sum)
    BLK.rec.m=aggregate(Number~TG.blk+Block.rec+Yr.rec+Mn.rec,subset(DAT.rec,Sex=="M"),sum)
    
    Zn.rec.f=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec,subset(DAT.rec,Sex=="F"),sum)
    Zn.rec.m=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec,subset(DAT.rec,Sex=="M"),sum)
    
    return(list(BLK.rel.f=BLK.rel.f,BLK.rel.m=BLK.rel.m,Zn.rel.f=Zn.rel.f,Zn.rel.m=Zn.rel.m,
                BLK.rec.f=BLK.rec.f,BLK.rec.m=BLK.rec.m,Zn.rec.f=Zn.rec.f,Zn.rec.m=Zn.rec.m))
  }
  
  if(BySex=="NO")
  {
    #assign tag group
    DAT$TG.blk=as.factor(with(DAT,paste(Block.rel,Yr.rel,Mn.rel,Age)))
    DAT$TG.zn=as.factor(with(DAT,paste(Rel.zone,Yr.rel,Mn.rel,Age)))
    
    #recaptures
    DAT.rec=subset(DAT,Recaptured=="Yes")
    
    
    #All release numbers
    BLK.rel=aggregate(Number~TG.blk+Block.rel+Yr.rel+Mn.rel+Age+Rel_FL,DAT,sum)
    Zn.rel=aggregate(Number~TG.zn+Rel.zone+Yr.rel+Mn.rel+Age+Rel_FL,DAT,sum)
    
    #Recapture numbers
    BLK.rec=aggregate(Number~TG.blk+Block.rec+Yr.rec+Mn.rec,DAT.rec,sum)
    Zn.rec=aggregate(Number~TG.zn+Rec.zone+Yr.rec+Mn.rec,DAT.rec,sum)
     
    return(list(BLK.rel=BLK.rel,Zn.rel=Zn.rel,BLK.rec=BLK.rec,Zn.rec=Zn.rec))
  }

}

Store.group=vector('list',length(Pop.din.sp))
names(Store.group)=Pop.din.sp
Store.group.size=Store.group
for(i in 1:length(Store.group))
{
  Store.group[[i]]=fn.group(subset(Tagging.pop.din,Species==Pop.din.sp[i] & !is.na(Rel_FL)),
    K.f=Gr[[i]][1],Linf.f=Gr[[i]][2],to.f=Gr[[i]][3],K.m=Gr[[i]][4],Linf.m=Gr[[i]][5],
    to.m=Gr[[i]][6],BySex="NO",MX.AGE=Mx.age[[i]])
}


# Export Data for Population modelling------------------------------------------------------------------

  #individual based model
setwd("C:/Matias/Analyses/Data_outs")
for(i in 1:length(Pop.din.sp))
{
  a=subset(Tagging.pop.din,Species==Pop.din.sp[i] & !is.na(DaysAtLarge) & !is.na(Rec.zone),
           select=c(Tag.no,DaysAtLarge,Rel.zone,Rec.zone))
  NmS=ifelse(Pop.din.sp[i]=="BW",'Dusky shark',
      ifelse(Pop.din.sp[i]=='WH','Whiskery shark',
      ifelse(Pop.din.sp[i]=='GM','Gummy shark',
      ifelse(Pop.din.sp[i]=='TK','Sandbar shark',NA))))
  write.csv(a,paste(getwd(),'/',NmS,'/',NmS,"_Con_tag_Ind.based.mod.csv",sep=""),row.names=F)
  
  #export data for mapping
  if(Pop.din.sp[i]%in% c("GM","WH"))
  {
    HnDDl="C:/Matias/Analyses/Movement rate estimation/Joint.estim_ind.base.mod/Show Gummy and whiskery outputs/"
    a=subset(Tagging.pop.din,Species==Pop.din.sp[i] & DaysAtLarge>=30 & !is.na(Rec.zone),
             select=c(Tag.no,DaysAtLarge,Lat.rels,Long.rels,Lat.rec,Long.rec,Rel.zone,Rec.zone))
    write.csv(a,paste(HnDDl,Pop.din.sp[i],"_Raw.Conv.Tag.csv",sep=""),row.names=F)
  }
  
}

  #Other stuff
setwd("C:/Matias/Data/Tagging/Pop dyn model/Conventional")
for(i in 1:length(Store.group))
{
  a=Store.group[[i]]
  for(p in 1:length(a)) write.csv(a[[p]],paste(names(Store.group)[i],"_",names(a)[p],"_","Conv.Tag.csv",sep=""),row.names=F)
}
  #b) Tag groups by size class
#note: two size groups, one from smallest tagged to size at maturity, the other for >size at maturity
for(i in 1:length(Store.group.size))
{
  DAT=subset(Tagging.pop.din,Species==Pop.din.sp[i] & !is.na(Rel_FL))  
  DAT$size.group=with(DAT,ifelse(Rel_FL<=SIZE[i],"Juvenile",ifelse(Rel_FL>SIZE[i],"Adult",NA)))  
  DAT$Number=1
  DAT$TG.zn=as.factor(with(DAT,paste(Rel.zone,Yr.rel,size.group))) #assign tag group  
  DAT.rec=subset(DAT,Recaptured=="Yes")  #recaptures
  Smlest=subset(DAT.rec,!Rel.zone== Rec.zone) #size of smallest individual moving to different zone
  Smallest.size=min(Smlest$Rel_FL)    
  
  #All release numbers
  Zn.rel=aggregate(Number~TG.zn+Rel.zone+Yr.rel+size.group,DAT,sum)
  Zn.rel.adul=subset(Zn.rel,size.group=="Adult")
  Zn.rel.juv=subset(Zn.rel,size.group=="Juvenile")
  
   
  #Recapture numbers
  Zn.rec=aggregate(Number~TG.zn+Rec.zone+Yr.rec+size.group,DAT.rec,sum)
  Zn.rec.adul=subset(Zn.rec,size.group=="Adult")
  Zn.rec.juv=subset(Zn.rec,size.group=="Juvenile")

  #all sizes
  DAT$TG.zn=as.factor(with(DAT,paste(Rel.zone,Yr.rel))) #assign tag group 
  Zn.rel=aggregate(Number~TG.zn+Rel.zone+Yr.rel,DAT,sum)
  DAT.rec=subset(DAT,Recaptured=="Yes")  
  Zn.rec=aggregate(Number~TG.zn+Rec.zone+Yr.rec,DAT.rec,sum)
  
  Store.group.size[[i]]=list(Zn.rel=Zn.rel,Zn.rel.adul=Zn.rel.adul,Zn.rel.juv=Zn.rel.juv,
                             Zn.rec=Zn.rec,Zn.rec.adul=Zn.rec.adul,Zn.rec.juv=Zn.rec.juv,
                             Smallest.size=Smallest.size)
}

#export
for(i in 1:length(Store.group.size))
{
  a=Store.group.size[[i]]
  NmS=ifelse(names(Store.group)[i]=="BW",'Dusky shark',
      ifelse(names(Store.group)[i]=='WH','Whiskery shark',
      ifelse(names(Store.group)[i]=='GM','Gummy shark',
      ifelse(names(Store.group)[i]=='TK','Sandbar shark',NA))))
   for(p in 1:length(a)) write.csv(a[[p]],paste('C:/Matias/Analyses/Data_outs/',NmS,'/',NmS,"_Con_tag_",names(a)[p],"_","Conv.Tag_size.csv",sep=""),row.names=F)
}



# Analysis of recapture data only for paper-----------------------------------------
do.paper=TRUE
if(do.paper)
{
  setwd("C:/Matias/Analyses/Conventional tagging/General movement/outputs")
  
  #keep recaptures only with data on date and position of release and recapture
  All.rel.Tagging=Tagging
  All.rel.Tagging=subset(All.rel.Tagging,Species%in%All.Species)
  
  Tagging=subset(Tagging,Recaptured=="Yes" & !is.na(Long.rec) & !is.na(DaysAtLarge))%>%
    mutate(Over.threshold=ifelse(DaysAtLarge>=min.daysAtLiberty,"YES","NO"))
  
  Prop.less.daysAtLiberty.threshold=Tagging%>%
    group_by(COMMON_NAME,Over.threshold)%>%
    tally()%>%
    spread(Over.threshold,n,fill=0)%>%
    data.frame
  write.csv(Prop.less.daysAtLiberty.threshold,"Paper/Prop.less.daysAtLiberty.threshold.csv")
  
  
  Tagging=Tagging%>%filter(Over.threshold=="YES")
  
  #Distance moved (in metres) and bearing (in 360 degrees; 0 degrees= North) 
  #note: follow shortest path (i.e. a Great Circle) and correct for corners
  Tagging$bearing=Tagging$dist.trav_m=NA #cannot use ifelse with bearing or distcosine functions
  Pos.cols=match(c("Lat.rels","Long.rels","Lat.rec","Long.rec"),names(Tagging))  #find the position of columns Lat.rec and Long.rec 
  for(i in 1:nrow(Tagging))
  {
    pp=Tagging[i,]
    if(sum(is.na(pp[Pos.cols]))==0) Tagging$bearing[i]=bearingRhumb(cbind(pp$Long.rels,pp$Lat.rels),cbind(pp$Long.rec,pp$Lat.rec))
    if(sum(is.na(pp[Pos.cols]))==0)
    {
      Tagging$dist.trav_m[i]=distCosine(cbind(pp$Long.rels,pp$Lat.rels),cbind(pp$Long.rec,pp$Lat.rec))
      Tagging$dist.trav_m.u[i]=Tagging$dist.trav_m[i]
      
      #Corrections for across land distances
      
      Tagging$dist.trav_m[i]=with(Tagging,
                                  
                                  #Released north of Exmouth                     
                                  ifelse(pp$Lat.rels>Exmouth[2] & pp$Lat.rec<=Exmouth[2] & pp$Lat.rec>Shark.bay[2],
                                         (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Exmouth)+distCosine(Exmouth,cbind(pp$Long.rec,pp$Lat.rec))),
                                         
                                         ifelse(pp$Lat.rels>Exmouth[2] & pp$Lat.rec<Shark.bay[2],
                                                (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Exmouth)+distCosine(Exmouth,Shark.bay)+
                                                   distCosine(Shark.bay,cbind(pp$Long.rec,pp$Lat.rec))),
                                                
                                                
                                                #Releases between Exmouth and Shark bay
                                                ifelse(pp$Lat.rels<=Exmouth[2]& pp$Lat.rels>Shark.bay[2] & pp$Lat.rec>Exmouth[2],
                                                       (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Exmouth)+distCosine(Exmouth,cbind(pp$Long.rec,pp$Lat.rec))),
                                                       
                                                       ifelse(pp$Lat.rels<=Exmouth[2]& pp$Lat.rels>Shark.bay[2] & pp$Lat.rec<Shark.bay[2],
                                                              (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Shark.bay)+distCosine(Shark.bay,cbind(pp$Long.rec,pp$Lat.rec))),
                                                              
                                                              
                                                              #Releases between Shark Bay and Cape Leuwin
                                                              ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Lat.rels>Cape.Leuwin[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec > Shark.bay[2] & pp$Lat.rec<=Exmouth[2],
                                                                     (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Shark.bay)+distCosine(Shark.bay,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                     
                                                                     ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Lat.rels>Cape.Leuwin[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec>Exmouth[2],
                                                                            (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Shark.bay)+distCosine(Shark.bay,Exmouth)+distCosine(Exmouth,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                            
                                                                            ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Lat.rels>Cape.Leuwin[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec<Cape.Leuwin[2] & pp$Lat.rec>Mid.point[2] &pp$Long.rec <Mid.point[1],
                                                                                   (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                   
                                                                                   ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Lat.rels>Cape.Leuwin[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec <Shark.bay[2] & pp$Long.rec >Mid.point[1] & pp$Long.rec <Streaky.Bay[1],
                                                                                          (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Cape.Leuwin)+distCosine(Cape.Leuwin,Mid.point)+
                                                                                             distCosine(Mid.point,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                          
                                                                                          ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Lat.rels>Cape.Leuwin[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec <Shark.bay[2] & pp$Long.rec >=Streaky.Bay[1],
                                                                                                 (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Cape.Leuwin)+distCosine(Cape.Leuwin,Mid.point)+
                                                                                                    distCosine(Mid.point,Streaky.Bay)+distCosine(Streaky.Bay,Spencer)+distCosine(Spencer,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                 
                                                                                                 
                                                                                                 #Releases between Cape Leuwin and MidPoint
                                                                                                 ifelse(pp$Lat.rels<=Cape.Leuwin[2] & pp$Lat.rels>Mid.point[2] & pp$Long.rels <Mid.point[1] & pp$Lat.rec>Cape.Leuwin[2],
                                                                                                        (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                        
                                                                                                        ifelse(pp$Lat.rels<=Cape.Leuwin[2] & pp$Lat.rels>Mid.point[2] & pp$Long.rels <Mid.point[1] & pp$Long.rec>Mid.point[1] &pp$Long.rec<=Streaky.Bay[1],
                                                                                                               (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Mid.point)+distCosine(Mid.point,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                               
                                                                                                               ifelse(pp$Lat.rels<=Cape.Leuwin[2] & pp$Lat.rels>Mid.point[2] & pp$Long.rels <Mid.point[1] & pp$Long.rec>Streaky.Bay[1],
                                                                                                                      (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Mid.point)+distCosine(Mid.point,Streaky.Bay)+distCosine(Streaky.Bay,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                      
                                                                                                                      #Releases east of MidPoint and west of Streaky Bay       
                                                                                                                      ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Long.rels >=Mid.point[1] & pp$Long.rels <Streaky.Bay[1] & 
                                                                                                                               pp$Long.rec> Streaky.Bay[1] & pp$Long.rec<=Spencer[1],
                                                                                                                             (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Streaky.Bay)+distCosine(Streaky.Bay,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                             
                                                                                                                             ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Long.rels >=Mid.point[1] & pp$Long.rels <Streaky.Bay[1] & pp$Long.rec> Spencer[1],
                                                                                                                                    (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Streaky.Bay)+distCosine(Streaky.Bay,Spencer)+distCosine(Spencer,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                                    
                                                                                                                                    ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Long.rels >=Mid.point[1] & pp$Long.rels <Streaky.Bay[1] &
                                                                                                                                             pp$Long.rec< Mid.point[1] & pp$Long.rec>=Cape.Leuwin[1]& pp$Lat.rec<=Cape.Leuwin[2],
                                                                                                                                           (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Mid.point)+distCosine(Mid.point,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                                           
                                                                                                                                           ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Long.rels >=Mid.point[1] & pp$Long.rels <Streaky.Bay[1] &
                                                                                                                                                    pp$Long.rec< Mid.point[1] & pp$Lat.rec>Cape.Leuwin[2]& pp$Lat.rec<=Shark.bay[2],
                                                                                                                                                  (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+distCosine(Cape.Leuwin,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                                                  
                                                                                                                                                  ifelse(pp$Lat.rels<=Shark.bay[2] & pp$Long.rels >=Mid.point[1] & pp$Long.rels <Streaky.Bay[1] &
                                                                                                                                                           pp$Lat.rec>Exmouth[2],
                                                                                                                                                         (distCosine(cbind(pp$Long.rels,pp$Lat.rels),Mid.point)+distCosine(Mid.point,Cape.Leuwin)+distCosine(Cape.Leuwin,Shark.bay)
                                                                                                                                                          +distCosine(Shark.bay,Exmouth)+distCosine(Exmouth,cbind(pp$Long.rec,pp$Lat.rec))),
                                                                                                                                                         dist.trav_m[i]))))))))))))))))))
    }
  }
  Tagging=subset(Tagging,dist.trav_m>0)
  
  
  # Calculate distance and position for 1 year @ liberty for those @ liberty > 1 year
  Tagging=Tagging%>%
    mutate(dist.trav_m_in.1.year=ifelse(DaysAtLarge>365,
                                        365*dist.trav_m/DaysAtLarge,dist.trav_m))
  Fin.pos=as.data.frame(destPointRhumb(p=cbind(Tagging$Long.rels,Tagging$Lat.rels),
                                       b=Tagging$bearing,
                                       d=Tagging$dist.trav_m_in.1.year))%>%
    rename(Long.rec_1.year=lon,
           Lat.rec_1.year=lat)
  Tagging=cbind(Tagging,Fin.pos)%>%
    mutate(Long.rec_1.year=ifelse(DaysAtLarge>365,Long.rec_1.year,Long.rec), 
           Lat.rec_1.year=ifelse(DaysAtLarge>365,Lat.rec_1.year,Lat.rec))
  Tagging$Areas.rec_1.year=with(Tagging,
                                ifelse(Lat.rec_1.year<=(-26.5) & Lat.rec_1.year>=(-33) & Long.rec_1.year <=116,"WCDGDLF",
                                       ifelse(Lat.rec_1.year<(-33) & Lat.rec_1.year>=(-41) & Long.rec_1.year <=116.5,"JASDGDLF.zone1",
                                              ifelse(Lat.rec_1.year<(-31) & Lat.rec_1.year>=(-41) & Long.rec_1.year >116.5,"JASDGDLF.zone2",
                                                     ifelse(Long.rec_1.year>=114 & Long.rec_1.year<=123.75 & Lat.rec_1.year >=(-23),"WANCSF",
                                                            ifelse(Long.rec_1.year>123.75 & Lat.rec_1.year >=(-18),"JANSF",
                                                                   ifelse(Lat.rec_1.year>=(-26.5) & Long.rec_1.year <=114,"Closed",
                                                                          NA)))))))
  Tagging$Areas.rec_1.year=with(Tagging,
                                ifelse(Lat.rec_1.year<=(-29) & Long.rec_1.year >129,"SA",Areas.rec_1.year))
  
  Tagging=Tagging%>%
    mutate(Areas.rec_1.year=ifelse(DaysAtLarge>365,Areas.rec_1.year,Areas.rec))        
  Tagging=Tagging%>%
    mutate(Areas.rec_1.year=ifelse(is.na(Areas.rec_1.year) & Lat.rels>(-29) & Lat.rec <(-19.3),"WANCSF",
                                   ifelse(is.na(Areas.rec_1.year) & Long.rels<=(124.071) & Long.rec >=(113.15),"JASDGDLF.zone1",
                                          Areas.rec_1.year)))
  
  
  #table releases and recaptures by zone within a year at liberty
  Tagging$Areas=factor(Tagging$Areas,
                       levels=c("JANSF","WANCSF","Closed","WCDGDLF","JASDGDLF.zone1","JASDGDLF.zone2","SA"))
  Tagging$Areas.rec=factor(Tagging$Areas.rec,
                           levels=c("JANSF","WANCSF","Closed","WCDGDLF","JASDGDLF.zone1","JASDGDLF.zone2","SA"))
  Tagging$Areas.rec_1.year=factor(Tagging$Areas.rec_1.year,
                                  levels=c("JANSF","WANCSF","Closed","WCDGDLF","JASDGDLF.zone1","JASDGDLF.zone2","SA"))
  numInt=10
  col.image =rev(gray(0:(numInt-2)/(numInt-2)))
  fn.rel.rec.zn=function(DaT) 
  {
    DaT=DaT%>%
      filter(!is.na(Areas))%>%
      filter(!is.na(Areas.rec_1.year)) 
    Rel_rec=DaT%>%
      count(Areas, Areas.rec_1.year,.drop = FALSE)%>%
      spread(Areas.rec_1.year,n,fill=0)%>%
      data.frame%>%
      mutate(N=rowSums(.[-1]))%>%
      mutate(across(-c(Areas,N), ~ . / !! as.symbol("N")))%>%
      mutate_at(vars(- !! c('Areas','N')), round, 3)%>%
      mutate_if(is.numeric , replace_na, replace = 0)
    return(Rel_rec)
  }
  Rel.Rec.zn=vector('list',length(SPECIES))
  for(i in 1:length(SPECIES)) Rel.Rec.zn[[i]]=fn.rel.rec.zn(DaT=subset(Tagging,Species==SPECIES[i]))
  
  #plot table
  colfunc <- colorRampPalette(c("white","grey60", "black"))
  #colfunc <- colorRampPalette(c("white","khaki2", "darkred"))
  fn.img=function(d,MAIN)
  {
    Nms=as.character(d$Areas)
    Nms=ifelse(Nms=="JASDGDLF.zone1","Zone1",ifelse(Nms=="JASDGDLF.zone2","Zone2",Nms))
    MAT=d%>%
      dplyr::select(-c(Areas,N))%>%
      as.matrix
    BRKS=seq(0,max(MAT,na.rm = T),length.out = 20)
    image.plot(1:ncol(MAT),1:nrow(MAT),t(MAT),xaxt='n',yaxt='n',ylab="",xlab='',
               main=MAIN,col=colfunc(length(BRKS)))
    axis(1,1:ncol(MAT),Nms,las=3,cex.axis=.95)
    axis(2,1:nrow(MAT),Nms,las=1,cex.axis=.95)
    grid(nx=ncol(MAT),ny=nrow(MAT),col="black")
    box()
  }
  
  tiff(file="Paper/Figure3_Release_recap_1yr.tiff",2400,2000, units="px",res=300,compression="lzw")
  par(mfcol=c(2,2),mar=c(4.5,6,1,1),mgp=c(1,.6,0))
  for(i in 1:length(SPECIES)) fn.img(d=Rel.Rec.zn[[i]],MAIN=names(SPECIES)[i])
  mtext("Released",2,outer = T,line=-1.5,las=3,cex=1.3)
  mtext("Recaptured",1,outer = T,line=-1,cex=1.3)
  dev.off()
  
  
  fn.rel.rec.zn_size.bin=function(DaT,biN) 
  {
    DaT=DaT%>%
      filter(!is.na(Areas))%>%
      filter(!is.na(Areas.rec_1.year))%>%
      mutate(FL.bin=biN*floor(Rel_FL/biN))
    
    TAB=table(DaT$FL.bin)
    N=names(TAB[TAB>5])
    
    par(mfrow=n2mfrow(length(N)),mar=c(4.5,6,1,1),mgp=c(1,.6,0))
    for(n in 1:length(N))
    {
      Rel_rec=DaT%>%
        filter(FL.bin==N[n])
      nn=nrow(Rel_rec)
      Rel_rec=Rel_rec%>%
        count(Areas, Areas.rec_1.year,.drop = FALSE)%>%
        spread(Areas.rec_1.year,n,fill=0)%>%
        data.frame%>%
        mutate(N=rowSums(.[-1]))%>%
        mutate(across(-c(Areas,N), ~ . / !! as.symbol("N")))%>%
        mutate_at(vars(- !! c('Areas','N')), round, 3)%>%
        mutate_if(is.numeric , replace_na, replace = 0)
      fn.img(d=Rel_rec,MAIN=paste(names(SPECIES)[i]," ",N[n]," cm FL"," (n=",nn,")",sep=""))
    }
    mtext("Released",2,outer = T,line=-1.5,las=3,cex=1.3)
    mtext("Recaptured",1,outer = T,line=-1,cex=1.3)
  }
  bins=c(20,10,10,10)
  pdf("Release_recap_1yr_by.size.pdf")
  for(i in 1:length(SPECIES)) fn.rel.rec.zn_size.bin(DaT=subset(Tagging,Species==SPECIES[i]),biN=bins[i])
  dev.off()
  
  
  #plot mean direction by size bin for 1 year at liberty
  fun1=function(d)
  {
    d=d%>%
      filter(!is.na(Lat.rec) | !is.na(Long.rec))%>%
      mutate(Rel_FL.bin=factor(Rel_FL.bin,levels=sort(unique(as.numeric(Rel_FL.bin)))),
             Bearing.rad=NISTdegTOradian(bearing))
    d=subset(d,!is.na(Rel_FL.bin))
    
    
    d=subset(d,!is.na(YrsAtLarge.bin))
    d=d%>%mutate(YrsAtLarge.bin2=ifelse(YrsAtLarge.bin==0,"<1 year @ liberty",
                                        ifelse(YrsAtLarge.bin==1,"1-2 years @ liberty",
                                               ifelse(YrsAtLarge.bin==2,"2-3 years @ liberty",
                                                      ifelse(YrsAtLarge.bin>2,">3 years @ liberty",NA)))),
                 YrsAtLarge.bin2=factor(YrsAtLarge.bin2,levels=c("<1 year @ liberty",
                                                                 "1-2 years @ liberty",
                                                                 "2-3 years @ liberty",
                                                                 ">3 years @ liberty")))      
    
    Uni.area=d%>%distinct(Block,.keep_all = T)%>%dplyr::select(Block,Areas)
    
    a=d%>%
      group_by(Block,Rel_FL.bin,YrsAtLarge.bin2)%>%
      summarise(Dir=mean(bearing))%>%
      mutate(Long=100+as.numeric(substr(Block,3,4)),
             Lat=-as.numeric(substr(Block,1,2)))%>%
      data.frame
    Fin.pos=as.data.frame(destPointRhumb(p=cbind(a$Long,a$Lat),b=a$Dir,d=100*1000))%>%
      rename(Lon.final=lon,
             Lat.final=lat)
    cbind(a,Fin.pos)%>%left_join(Uni.area,by='Block')%>%
      ggplot(aes(Long,Lat))+
      geom_segment(aes(x = Long, y = Lat, xend = Lon.final, yend = Lat.final, colour = Rel_FL.bin),
                   arrow = arrow(length = unit(0.25, "cm")))+
      facet_wrap( vars(YrsAtLarge.bin2))+xlab("Longitude") + ylab("Latitude")+
      theme_bw() + labs(fill = "Release size bin (cm)")
    ggsave(paste("Mean.direction_",unique(d$COMMON_NAME),".tiff",sep=''), width = 10,height = 10, dpi = 300,compression = "lzw")
    
  }
  for(i in 1:length(SPECIES)) fun1(d=Tagging%>%filter(Species==SPECIES[i] & Recaptured=='Yes'))
  
  
  #Rate of Movement (km/day)
  Tagging$ROM=with(Tagging,ifelse(DaysAtLarge>0,(dist.trav_m/(1000*DaysAtLarge)),NA))
  Tagging$ROM=with(Tagging,ifelse(ROM<0,NA,ROM))
  
  
  #Speed (m/s)
  Tagging$speed=with(Tagging,ifelse(ROM>0,(dist.trav_m/(DaysAtLarge*60*60*24)),NA))
  
  
  #remove nonsense records                                                   
  # speed.min=0.000001  #minimum speed for straight linear movement. GET FROM ACOUSTIC TAGGING
  # speed.max=0.2  #maximum speed

  n.records.including.nonsense=nrow(Tagging)
  
  #remove records with negative speed (reporting issues) or less than threshold
  Tagging=subset(Tagging,!is.na(speed))
  n.records.without.nonsense=nrow(Tagging)
  dropped=100-round(100*n.records.without.nonsense/n.records.including.nonsense)  #percent dropped records after removing nonsense
  
  
  # Proportion of time within zones
  #note: this is used for Risk Assessment of dusky and sandbar
  do.prop.time.zn=FALSE
  if(do.prop.time.zn)
  {
    Prop.time=Tagging[,match(c("Species","FINTAGNO","Lat.rels","Long.rels","Lat.rec","Long.rec",
                               "DATE_REL","DATE_CAPTR","dist.trav_m","DaysAtLarge"),names(Tagging))]
    names(Prop.time)[match(c("FINTAGNO","Lat.rels","Long.rels","Lat.rec","Long.rec",
                             "DATE_REL","DATE_CAPTR","dist.trav_m","DaysAtLarge"),names(Prop.time))]=c("TagCode","Latitude.prev",
                                                                                                       "Longitude.prev","Latitude","Longitude","Date.local.prev","Date.local","Dist.moved.conseq.det","days.conseq.det")
    
    #Add fishing zones
    Prop.time$zone=as.character(with(Prop.time,
                                     ifelse(Longitude>=116.5 & Latitude<=(-26) & Longitude<=129,"Zone2",
                                            ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
                                                   ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
                                                          ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
                                                                 ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
                                                                        ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",
                                                                               ifelse(Longitude>129 & Latitude<=(-26),"SA",NA)))))))))
    
    Prop.time$zone=with(Prop.time,ifelse(Latitude>(-33) & 
                                           Latitude<=(-31) & Longitude>=114.8476 & Longitude<116,"Closed.metro",zone))
    
    Prop.time$zone.prev=as.character(with(Prop.time,
                                          ifelse(Longitude.prev>=116.5 & Latitude.prev<=(-26) & Longitude.prev<=129,"Zone2",
                                                 ifelse(Longitude.prev<116.5 & Latitude.prev<=(-33),"Zone1",
                                                        ifelse(Latitude.prev>(-33) & Latitude.prev<=(-26) & Longitude.prev<116.5,"West",
                                                               ifelse(Latitude.prev>(-26) & Longitude.prev<114.833,"Closed.ningaloo",
                                                                      ifelse(Latitude.prev>(-22) & Longitude.prev>=114.833 & Longitude.prev<123.75,"North",
                                                                             ifelse(Latitude.prev>(-22) & Longitude.prev>=123.75,"Joint",
                                                                                    ifelse(Longitude.prev>129 & Latitude.prev<=(-26),"SA",
                                                                                           NA)))))))))
    
    Prop.time$zone.prev=with(Prop.time,ifelse(Latitude.prev>(-33) & 
                                                Latitude.prev<=(-31) & Longitude.prev<116,"Closed.metro",zone.prev))
    
    
    Prop.time$same.zone=with(Prop.time,ifelse(!(zone==zone.prev),"N","Y"))
    
    #Total time monitored
    Prop.time$ReleaseDate=Prop.time$Date.local.prev
    shks=unique(Prop.time$TagCode)
    Total.time.monitored=Prop.time[,match(c("TagCode","days.conseq.det"),names(Prop.time))]
    names(Total.time.monitored)[2]="days.mon"
    
    #Time by zone
    Time.mon.zone=vector('list',length(shks))
    Tim.mon.zn=function(SHK)
    {  
      dat1=subset(Prop.time,TagCode==SHK)
      dat1=dat1[!duplicated(dat1$TagCode),]
      Stay.in.zone=subset(dat1,same.zone=="Y"| is.na(same.zone))
      Leave.zone=subset(dat1,same.zone=="N")
      if(nrow(Stay.in.zone)>0) Stay.in.zone$Interpolate="N"
      
      if(nrow(Leave.zone)>0)
      {
        Store=vector('list',nrow(Leave.zone))
        for(s in 1:nrow(Leave.zone))
        {
          a=Leave.zone[s,]
          n=as.numeric(a$Date.local-a$Date.local.prev)
          b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),c(Longitude,Latitude),n=n, addStartEnd=F))
          
          #Add corners
          if(a$zone.prev=="North" & a$zone%in%c("Closed.metro","Zone1","West"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$zone.prev=="North" & a$zone%in%c("Zone2","SA"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,Shark.bay)/1000
            n3=distCosine(Shark.bay,Cape.Leuwin)/1000
            n4=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b1=gcIntermediate(Exmouth,Shark.bay,n=n2, addStartEnd=F)
            b2=gcIntermediate(Shark.bay,Cape.Leuwin,n=n3, addStartEnd=F)
            b3=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
          
          if(a$zone.prev%in%c("Closed.ningaloo") & a$zone%in%c("Zone2","SA"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,Cape.Leuwin)/1000    
            n3=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=gcIntermediate(Shark.bay,Cape.Leuwin,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$zone.prev%in%c("Closed.ningaloo") & a$zone%in%c("Zone1","Closed.metro","West"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          
          if(a$zone.prev%in%c("Closed.metro","West") & a$zone%in%c("Zone2","SA"))        
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$zone.prev%in%c("Closed.metro","Zone1","West") & a$zone%in%c("Closed.ningaloo"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          
          if(a$zone.prev=="Closed.ningaloo" & a$zone%in%c("North"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Exmouth)/1000    #divide time proportional to distance
            n2=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Exmouth,n=n1, addStartEnd=F))
            b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b3)
          }
          
          if(a$zone.prev%in%c("Zone2","SA") & a$zone%in%c("Closed.ningaloo"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Shark.bay)/1000
            n3=distCosine(Shark.bay,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
            b2=with(a,gcIntermediate(Shark.bay,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b1,b2)
          }
          if(a$zone.prev%in%c("Zone2","SA") & a$zone%in%c("Closed.metro","West"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=with(a,gcIntermediate(Cape.Leuwin,c(Longitude,Latitude),n=n2, addStartEnd=F))
            b=rbind(b,b1)
          }
          if(a$zone.prev%in%c("Zone2","SA") & a$zone%in%c("North"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Cape.Leuwin)/1000    #divide time proportional to distance
            n2=distCosine(Cape.Leuwin,Shark.bay)/1000
            n3=distCosine(Shark.bay,Exmouth)/1000
            n4=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3+n4
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            n4=round(n*(n4/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Cape.Leuwin,n=n1, addStartEnd=F))
            b1=gcIntermediate(Cape.Leuwin,Shark.bay,n=n2, addStartEnd=F)
            b2=gcIntermediate(Shark.bay,Exmouth,n=n3, addStartEnd=F)
            b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n4, addStartEnd=F))
            b=rbind(b,b1,b2,b3)
          }
          
          if(a$zone.prev%in%c("Closed.metro","West") & a$zone%in%c("North"))
          {
            n=as.numeric(a$Date.local-a$Date.local.prev)
            n1=distCosine(c(a$Longitude.prev,a$Latitude.prev),Shark.bay)/1000    #divide time proportional to distance
            n2=distCosine(Shark.bay,Exmouth)/1000
            n3=distCosine(Exmouth,c(a$Longitude,a$Latitude))/1000
            TN=n1+n2+n3
            n1=round(n*(n1/TN))
            n2=round(n*(n2/TN))
            n3=round(n*(n3/TN))
            b=with(a,gcIntermediate(c(Longitude.prev,Latitude.prev),Shark.bay,n=n1, addStartEnd=F))
            b2=gcIntermediate(Shark.bay,Exmouth,n=n2, addStartEnd=F)
            b3=with(a,gcIntermediate(Exmouth,c(Longitude,Latitude),n=n3, addStartEnd=F))
            b=rbind(b,b2,b3)
          }
          
          x=a[rep(seq_len(nrow(a)), n),]
          this.na=match(c("Latitude.prev","Longitude.prev","same.zone"),names(x))
          x[,this.na]=NA
          if(!nrow(x)==nrow(b)) x=x[1:nrow(b),]
          x$Longitude=b[,1]
          x$Latitude=b[,2]
          Day.seq=seq(1,n,1)*24*60*60
          if(!length(Day.seq)==nrow(b)) Day.seq=Day.seq[1:nrow(b)]
          x$Date.local=a$Date.local.prev+Day.seq       
          
          Store[[s]]=x
        }  
        Store=do.call(rbind,Store)
        Store$zone=as.character(with(Store,
                                     ifelse(Longitude>=116.5 & Latitude<=(-26) & Longitude<=129,"Zone2",
                                            ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
                                                   ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
                                                          ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
                                                                 ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
                                                                        ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",
                                                                               ifelse(Latitude<=(-26) & Longitude>129,"SA",NA)))))))))
        Store$zone=with(Store,ifelse(Latitude>(-33) & 
                                       Latitude<=(-31) & Longitude<116,"Closed.metro",zone))
        Store$Interpolate="Y"
      }
      
      if(nrow(Leave.zone)==0) All=Stay.in.zone
      if(nrow(Leave.zone)>0) All=rbind(Store,Stay.in.zone)
      All$Dupl=with(All,paste(Date.local,zone))
      All=All[!duplicated(All$Dupl),]
      All=All[order(All$TagCode,All$Date.local),]
      
      return(All)
    }
    system.time(for (i in 1:length(shks))Time.mon.zone[[i]]=Tim.mon.zn(shks[i]))
    Time.mon.zone=do.call(rbind,Time.mon.zone)  
    
    
    #aggregate number of days within each zone
    Time.mon.zone$Days=1
    
    Species.time.zone=aggregate(Days~zone+TagCode+Species,Time.mon.zone,sum,na.rm=T)
    Species.time.zone=subset(Species.time.zone,!TagCode=="")
    aa=Species.time.zone[duplicated(Species.time.zone$TagCode),]
    aa=unique(aa$TagCode)
    
    Species.time.zone=merge(Species.time.zone,Total.time.monitored,by="TagCode",all.x=T)
    Species.time.zone$Days=with(Species.time.zone,ifelse(!TagCode%in%aa,days.mon,Days))
    
    
    Species.time.zone$prop.time.in.zn=Species.time.zone$Days/Species.time.zone$days.mon
    Species.time.zone=Species.time.zone[,-match("Days",names(Species.time.zone))]  
    
    #Show.proportions="YES"  #show proportion of days monitored by zone
    Show.proportions="NO"  #show number of days by zone
    if(Show.proportions=="NO") Species.time.zone$prop.time.in.zn=round(with(Species.time.zone,
                                                                            days.mon *prop.time.in.zn))
    
    Species.time.zone=reshape(Species.time.zone,v.names = "prop.time.in.zn", idvar = c("TagCode","Species","days.mon"),
                              timevar = "zone", direction = "wide")
    names(Species.time.zone)=c("TagCode","Species","Tot.days.mon.","Zone2","SA","Zone1",
                               "Closed.ningaloo","West","Closed.metro","North","Joint")
    
    Species.time.zone=Species.time.zone[order(Species.time.zone$Species,Species.time.zone$TagCode),
                                        match(c("Species","TagCode","Joint","North","Closed.ningaloo","Closed.metro",
                                                "West","Zone1","Zone2","SA","Tot.days.mon."),names(Species.time.zone))]  
    
    Species.time.zone[is.na(Species.time.zone)]=0
    ID.thiS=match(c("Joint","North","Closed.ningaloo","Closed.metro","West",
                    "Zone1","Zone2","SA"),names(Species.time.zone))
    Species.time.zone[,ID.thiS]=Species.time.zone[,ID.thiS]
    
    #add zone release
    Time.mon.zone$Zone.rel=Time.mon.zone$zone.prev
    dummy=subset(Time.mon.zone,select=c(TagCode,Species,Zone.rel))
    dummy=dummy[!duplicated(paste(dummy$TagCode,dummy$Species,dummy$Zone.rel)),]
    Species.time.zone=merge(Species.time.zone,dummy,by=c("TagCode","Species"),all.x=T)
    
    Species.time.zone=Species.time.zone[order(Species.time.zone$Species,Species.time.zone$Zone.rel,
                                              -Species.time.zone$Tot.days.mon.),match(c("Species","Tot.days.mon.","TagCode","Zone.rel",
                                                                                        "Joint","North","Closed.ningaloo","West" ,"Closed.metro","Zone1","Zone2","SA"),colnames(Species.time.zone))]
    
    
    
    #plot proportions
    # Zns=c("Joint","North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2","SA")
    # Zns.leg=c("Joint","North","Ningaloo","West","Metro closure","Zone1","Zone2","SA")
    # COL.prop=c("darkolivegreen1","aquamarine3","chartreuse","darkgreen","cadetblue",
    #            "deepskyblue","blue4","cornflowerblue")
    # names(COL.prop)=Zns
    
    
    #Zones (same colors as in acoustic tagging)
    Zns=c("Joint","North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2","SA")
    Zns.leg=c("Joint","WANCSF","Ningaloo","WCDGDLF","Metro closure","Zone1","Zone2","SA")
    COL.prop=c("darkolivegreen1","lightseagreen","seagreen4","lightgreen","olivedrab4","olivedrab3","mediumseagreen","cornflowerblue")
    names(COL.prop)=Zns  
    
    
    # plot(1:8,col=COL.prop,pch=19,cex=3,ylim=c(0,9))
    # text(1:8,Zns)
    # points(1:8,0:7,col=COL.prop,pch=19,cex=3)
    
    #MIN.TIM.MON=30     #minimum time monitored displayed
    #MIN.TIM.MON=365
    MIN.TIM.MON=1
    Species.time.zone_1_30=subset(Species.time.zone,Tot.days.mon.>1 & Tot.days.mon.<=30) 
    Species.time.zone_30_365=subset(Species.time.zone,Tot.days.mon.>30 & Tot.days.mon.<=365) 
    Species.time.zone_365=subset(Species.time.zone,Tot.days.mon.>365) 
    
    
    What.prop="Just zone"
    #What.prop="Zone and growth"
    
    #a. Proportion of time per zone
    if(What.prop=="Just zone")
    {
      fn.plot.prop=function(what)
      {
        a=subset(Species.time.zone,Species==what)
        id=match(Zns,names(a))
        a[,id]=a[,id]/rowSums(a[,id])
        crap=as.data.frame(as.matrix(COL.prop))
        names(crap)="Col"
        crap$Zone.rel=rownames(crap)
        a=merge(a,crap,by="Zone.rel")
        a$Sort=with(a,ifelse(Zone.rel=="North",1,
                             ifelse(Zone.rel=="Closed.ningaloo",2,
                                    ifelse(Zone.rel=="West",3,
                                           ifelse(Zone.rel=="Closed.metro",4,
                                                  ifelse(Zone.rel=="Zone1",5,
                                                         ifelse(Zone.rel=="Zone2",6,7)))))))
        a=a[order(a$Sort),]
        par(mai=c(.45,.5,.13,.6),oma = c(.15, .4, .6, 0),xpd=T,mgp=c(1,.425,0))
        r=barplot(t(a[,id]), col = COL.prop,horiz=T,beside=F,yaxt='n',cex.axis=1.1)
        legend("top",Zns.leg,bty='n',pt.cex=2,pch=15,col=COL.prop,horiz=T,inset=c(0,-.04),
               cex=0.75)
        
        
        box()
        points(rep(-0.03,length(r)),r,pch=15,cex=1.5,col=as.character(a$Col))
        #axis(2,r,a$Zone.rel,las=1,cex.axis=1.25)
        mtext("Proportion of time",1,line=1.4,cex=1.5)
        mtext("Release zone",2,line=1.5,cex=1.5)
        ll=seq(1,length(r),2)
        axis(4,r,F,tck=-0.005)
        axis(4,r[ll],a$Tot.days.mon.[ll],las=1,cex.axis=0.5,tck=-0.01)
        mtext("Days at liberty",4,line=2,cex=1.5)
      }
      hnDl="Prop.time.zone/"
      tiff(file=paste(hnDl,"dusky.tiff",sep=""),width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
      fn.plot.prop(Pop.din.sp[1])
      dev.off()
      
      tiff(file=paste(hnDl,"whiskery.tiff",sep=""),width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
      fn.plot.prop(Pop.din.sp[4])
      dev.off()
      
      tiff(file=paste(hnDl,"gummy.tiff",sep=""),width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
      fn.plot.prop(Pop.din.sp[3])
      dev.off()
      
      tiff(file=paste(hnDl,"sandbar.tiff",sep=""),width = 2400, height = 2000,units = "px", res = 300,compression = "lzw")
      fn.plot.prop(Pop.din.sp[2])
      dev.off()
      
      
    }
    
    
    #b. Proportion of time by zone with growth
    if(What.prop=="Zone and growth")
    {
      fn.plot.prop.growth=function(what,SKLER,SX,whr,CX,iidd,MAX)
      {
        a=subset(Species.time.zone,Species==what)
        id=match(Zns,names(a))
        ZONAS=colSums(a[,id])
        dropp=which(ZONAS==0)
        if(length(dropp)>0) a=a[,-match(names(dropp),names(a))]
        
        id=match(Zns,names(a))
        id=id[!is.na(id)]
        
        #add size  
        dummy=subset(Tagging,FINTAGNO%in%unique(a$TagCode), select=c(FINTAGNO,Rel_FL,Sex,Lat.rec,Long.rec))
        names(dummy)=c("TagCode","FL.rel","Sex","Latitude","Longitude")
        
        #Add fishing zones
        dummy$zone=as.character(with(dummy,
                                     ifelse(Longitude>=116.5 & Latitude<=(-26) & Longitude<=129,"Zone2",
                                            ifelse(Longitude<116.5 & Latitude<=(-33),"Zone1",
                                                   ifelse(Latitude>(-33) & Latitude<=(-26) & Longitude<116.5,"West",
                                                          ifelse(Latitude>(-26) & Longitude<114.833,"Closed.ningaloo",
                                                                 ifelse(Latitude>(-22) & Longitude>=114.833 & Longitude<123.75,"North",
                                                                        ifelse(Latitude>(-22) & Longitude>=123.75,"Joint",
                                                                               ifelse(Longitude>129 & Latitude<=(-26),"SA",NA)))))))))
        
        dummy$zone=with(dummy,ifelse(Latitude>(-33) & 
                                       Latitude<=(-31) & Longitude>=114.8476 & Longitude<116,"Closed.metro",zone))
        dummy=subset(dummy,select=c(TagCode, FL.rel, Sex, zone))
        names(dummy)[match("zone",names(dummy))]="zone.rec"
        
        #growth pars
        Gr.par=Gr[[match(what,names(Gr))]]
        K.f=Gr.par[match("K.f",names(Gr.par))];K.m=Gr.par[match("K.m",names(Gr.par))]
        Linf.f=Gr.par[match("Linf.f",names(Gr.par))];Linf.m=Gr.par[match("Linf.m",names(Gr.par))]
        to.f=Gr.par[match("to.f",names(Gr.par))];to.m=Gr.par[match("to.m",names(Gr.par))]
        
        #add size at recapture
        a=merge(a,dummy,by="TagCode")
        a=subset(a,!is.na(FL.rel))
        a$Sex=with(a,ifelse(Sex=="U","F",Sex))
        a$Age.rel=with(a,ifelse(Sex=="F",to.f-(1/K.f)*log(1-(FL.rel/Linf.f)),
                                ifelse(Sex=="M",to.m-(1/K.m)*log(1-(FL.rel/Linf.m)),NA))) 
        a=subset(a,!is.na(Age.rel))
        a$Age.rec=a$Age.rel+(a$Tot.days.mon./365)
        a$FL.rec=ifelse(a$Sex=="F",Linf.f*(1-exp(-K.f*(a$Age.rec-to.f))),
                        Linf.m*(1-exp(-K.m*(a$Age.rec-to.m))))
        a$Sort=with(a,ifelse(Zone.rel=="Joint",1,
                             ifelse(Zone.rel=="North",2,
                                    ifelse(Zone.rel=="Closed.ningaloo",3,
                                           ifelse(Zone.rel=="West",4,
                                                  ifelse(Zone.rel=="Closed.metro",5,
                                                         ifelse(Zone.rel=="Zone1",6,
                                                                ifelse(Zone.rel=="Zone2",7,8))))))))
        a=subset(a,Sex==SX)
        RL.Zns=sort(unique(a$Sort))
        These.Leg=which(names(COL.prop)%in%unique(a$zone.rec))
        
        for(tt in 1:length(RL.Zns))
        {
          bb=subset(a,Sort==RL.Zns[tt])  
          bb=bb[!duplicated(bb$TagCode),]
          bb=bb[order(bb$FL.rel),]    
          dd=bb[,match(c("Sex","Zone.rel","FL.rel","zone.rec","FL.rec"),names(bb))]
          DD=dd[,match(c("FL.rel","FL.rec"),names(dd))]
          DD=DD/MAX
          ee=bb[,id]  
          ss=unique(dd$Zone.rel)
          if(ss=="North") Srt=c("North","Closed.ningaloo","Joint","West","Closed.metro","Zone1","Zone2","SA")
          
          if(ss=="Closed.ningaloo") Srt=c("Closed.ningaloo","North","Joint","West","Closed.metro","Zone1","Zone2","SA")
          if(ss=="West") Srt=c("West","Closed.ningaloo","North","Joint","Closed.metro","Zone1","Zone2","SA")
          if(ss=="Closed.metro") Srt=c("Closed.metro","West","Closed.ningaloo","North","Joint","Zone1","Zone2","SA")
          if(ss=="Zone1") Srt=c("Zone1","Closed.metro","West","Closed.ningaloo","North","Joint","Zone2","SA")
          if(ss=="Zone2") Srt=c("Zone2","Zone1","Closed.metro","West","Closed.ningaloo","North","Joint","SA")
          if(ss=="SA") Srt=c("SA","Zone2","Zone1","Closed.metro","West","Closed.ningaloo","North","Joint")
          cols=match(Srt,colnames(ee))
          cols=cols[!is.na(cols)]
          ee=ee[,cols]    
          ee=t(ee)
          
          CLSS=COL.prop[match(rownames(ee),names(COL.prop))]  
          r=barplot(ee,col = CLSS, horiz=T,beside=F,yaxt='n',cex.axis=1.25,xlim=c(-50,max(colSums(ee))*1.05))
          box()
          points(rep(0,length(r)),r,pch=19,cex=DD[,1]*SKLER,col="red")
          points(colSums(ee),r,pch=19,cex=DD[,2]*SKLER,col="red")
          mtext(paste("Released in ",Zns.leg[RL.Zns[tt]]),2,cex=1.2,line=0.5)
          if(tt==iidd)
          {
            qq=MAX
            qqq=round(rep(qq,4)*seq(.2,1,.25))
            legend(whr,paste(qqq),bty='n',pt.cex=(qqq/qq)*SKLER,pch=19,col="red",
                   title="FL (cm)",cex=1.5)
          }   
        }
        
        #add zone legend
        par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
        plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
        legend("top",Zns.leg[These.Leg],bty='n',pt.cex=2,pch=15,
               col=COL.prop[These.Leg],xpd = TRUE,horiz=T,inset=c(0,0),cex=1.25)
        
        mtext("Days at liberty",1,line=-2,cex=1.5,outer=T)
        
      }
      
      tiff(file="Prop.time.zone/Proportion.zone.whiskery_F.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[4],SKLER=2,"F","topright",CX=1.25,iidd=4,160)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.whiskery_M.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[4],SKLER=2,"M","topright",CX=1.25,iidd=4,160)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.dusky_F.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      #par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(3, .5, .1, 0),mgp=c(1,.7,0))
      par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[1],SKLER=2,"F","bottomright",CX=1.25,iidd=2,200)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.dusky_M.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[1],SKLER=2,"M","bottomright",CX=1.25,iidd=4,200)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.gummy_F.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(2,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[3],SKLER=2,"F","bottomright",CX=1.25,iidd=2,160)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.gummy_M.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(2,1),mai=c(.1,.3,.25,.1),oma = c(4, 1, 1, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[3],SKLER=2,"M","bottomright",CX=1.25,iidd=1,160)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.sandbar_F.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(3,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 2, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[2],SKLER=3,"F","topright",CX=1.5,iidd=4,160)
      dev.off()
      
      tiff(file="Prop.time.zone/Proportion.zone.sandbar_M.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
      par(mfcol=c(3,2),mai=c(.1,.3,.2,.1),oma = c(4, 1, 2, 1),mgp=c(1,.7,0))
      fn.plot.prop.growth(Pop.din.sp[2],SKLER=3,"M","topright",CX=1.5,iidd=4,160)
      dev.off()
      
    }
    #Connectivity plot of overall movement
    
    #Image option
    Image.mig.fn=function(D,what,SP,CE,CE1,OFF)
    {
      a=subset(D,Species==what)
      if(what=="BW") a=subset(a,!Zone.rel=="West")  #remove single recapture from "WCDDLF" for dusky
      a=a[!duplicated(a$TagCode),]
      id=match(Zns,names(a))
      a[,id]=a[,id]/rowSums(a[,id])
      crap=as.data.frame(as.matrix(COL.prop))
      names(crap)="Col"
      crap$Zone.rel=rownames(crap)
      a=merge(a,crap,by="Zone.rel")
      
      a$Sort=with(a,ifelse(Zone.rel=="Joint",7,
                           ifelse(Zone.rel=="North",6,
                                  ifelse(Zone.rel=="Closed.ningaloo",5,
                                         ifelse(Zone.rel=="West",4,
                                                ifelse(Zone.rel=="Closed.metro",3,
                                                       ifelse(Zone.rel=="Zone1",2,1)))))))
      
      a$Zone.rel=as.character(a$Zone.rel)
      
      a$Zone.rel=with(a,ifelse(Zone.rel=="Closed.ningaloo","5.Ningaloo",
                               ifelse(Zone.rel=="Closed.metro","3.Metro",
                                      ifelse(Zone.rel=="Zone1","2.Zn1",
                                             ifelse(Zone.rel=="Zone2","1.Zn2",
                                                    ifelse(Zone.rel=="North","6.North",
                                                           ifelse(Zone.rel=="Joint","7.Joint",      
                                                                  "4.West")))))))
      a$Zone.rel=as.factor(a$Zone.rel)
      MATRX=aggregate(cbind(Joint,North,Closed.ningaloo,West,Closed.metro,Zone1,Zone2,SA)~Zone.rel,a,sum)
      MATRX[,2:ncol(MATRX)]=MATRX[,2:ncol(MATRX)]/rowSums(MATRX[,2:ncol(MATRX)])
      
      colnames(MATRX)[match(c("Joint","North","Closed.ningaloo","West","Closed.metro","Zone1","Zone2","SA"),
                            colnames(MATRX))]=c("7.Joint","6.North","5.Ningaloo","4.West","3.Metro","2.Zn1","1.Zn2","1.SA")
      MATRX=MATRX[,match(c("Zone.rel","1.SA","1.Zn2","2.Zn1","3.Metro","4.West","5.Ningaloo","6.North","7.Joint"),names(MATRX))]
      MATRX$Zone.rel=as.character(MATRX$Zone.rel)
      
      ALL.zns=c("1.SA","1.Zn2","2.Zn1","3.Metro","4.West","5.Ningaloo","6.North","7.Joint")
      ID=ALL.zns[which(!ALL.zns%in%MATRX$Zone.rel)]
      if(length(ID)>0)
      {
        ss=as.data.frame(matrix(nrow=length(ID),ncol=ncol(MATRX)))
        names(ss)=names(MATRX)
        ss[,]=0
        ss$Zone.rel=ID
        MATRX=rbind(MATRX,ss)
        MATRX=MATRX[order(MATRX$Zone.rel),]
      }
      MTRX=t(MATRX[,-1])  
      image(1:nrow(MTRX),1:ncol(MTRX),as.matrix(MTRX),col =couleurs,breaks=BREAKS,xaxt='n',yaxt='n',ylab="",xlab="")
      
      axis(2,1:ncol(MTRX),F,las=1)
      axis(1,1:nrow(MTRX),F)
      ZN.x=rev(Zns.leg)
      ZN.x=ifelse(ZN.x=="Metro closure","Metro",ifelse(ZN.x=="Joint","JANSF",ZN.x))
      if(what%in%c("BW", "TK")) axis(2,1:ncol(MTRX),ZN.x,las=1,cex.axis=CE)
      #if(what%in%c("TK","WH")) axis(1,1:nrow(MTRX),ZN.x,cex.axis=CE1)
      if(what%in%c("TK","WH")) text(1:ncol(MTRX),par("usr")[3]-0.2,ZN.x,cex=CE1,srt=45,pos=2,xpd=T,offset = OFF)
      box()
      mtext(SP,3)
      lines(.5:8.5,.5:8.5,lwd=1.5,col="grey60")
    }
    
    N.int=50
    do.col="N"
    colfunc <- colorRampPalette(c("navy", "cadetblue","white"))
    if(do.col=="Y")couleurs=rev(colfunc(N.int))
    if(do.col=="N")couleurs=rev(gray(seq(0,1,length=N.int)))
    BREAKS=seq(0,1,length.out=N.int+1)
    
    CeX=.8
    CeX1=.8
    
    #1 to 30 days
    tiff(file="Paper/FigA.3_1_to_30_days.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,2),mar = c(2.75, 1, 1, 1),oma=c(3,4,0,0),xpd=T,mgp=c(1,.7,0))
    Image.mig.fn(Species.time.zone_1_30,Pop.din.sp[1],"Dusky shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_1_30,Pop.din.sp[2],"Sandbar shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_1_30,Pop.din.sp[3],"Gummy shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_1_30,Pop.din.sp[4],"Whiskery shark",CeX,CeX1,0.175)
    color.legend(8.5,1,9,8,BREAKS[seq(1,length(BREAKS),5)],rect.col=couleurs,gradient="y",col=1,cex=0.7)
    mtext("Released",2,outer=T,line=2.5,cex=1.5)
    mtext("Recaptured",1,outer=T,line=1,cex=1.5)
    dev.off()
    
    
    #30 to 365 days
    tiff(file="Paper/FigA.3_30_to_365_days.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,2),mar = c(2.75, 1, 1, 1),oma=c(3,4,0,0),xpd=T,mgp=c(1,.7,0))
    Image.mig.fn(Species.time.zone_30_365,Pop.din.sp[1],"Dusky shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_30_365,Pop.din.sp[2],"Sandbar shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_30_365,Pop.din.sp[3],"Gummy shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_30_365,Pop.din.sp[4],"Whiskery shark",CeX,CeX1,0.175)
    color.legend(8.5,1,9,8,BREAKS[seq(1,length(BREAKS),5)],rect.col=couleurs,gradient="y",col=1,cex=0.7)
    mtext("Released",2,outer=T,line=2.5,cex=1.5)
    mtext("Recaptured",1,outer=T,line=1,cex=1.5)
    dev.off()
    
    
    #>365 days
    tiff(file="Paper/FigA.3_more.than.365_days.tiff",width = 2000, height = 2000,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,2),mar = c(2.75, 1, 1, 1),oma=c(3,4,0,0),xpd=T,mgp=c(1,.7,0))
    Image.mig.fn(Species.time.zone_365,Pop.din.sp[1],"Dusky shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_365,Pop.din.sp[2],"Sandbar shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_365,Pop.din.sp[3],"Gummy shark",CeX,CeX1,0.175)
    Image.mig.fn(Species.time.zone_365,Pop.din.sp[4],"Whiskery shark",CeX,CeX1,0.175)
    color.legend(8.5,1,9,8,BREAKS[seq(1,length(BREAKS),5)],rect.col=couleurs,gradient="y",col=1,cex=0.7)
    mtext("Released",2,outer=T,line=2.5,cex=1.5)
    mtext("Recaptured",1,outer=T,line=1,cex=1.5)
    dev.off()
    
    
  }
  
  #4.3. Summary statistics
  function.stats=function(SPEC)
  {
    datos=subset(Tagging,Species==SPEC)
    AREAS=unique(datos$Areas)
    
    COR.Dist.Time=Bearing.summary=STATS.COND=list()
    for(i in 1:length(AREAS))
    {
      datos1=subset(datos,Areas==AREAS[i])
      
      #.. dats at large vs distance    
      COR.Dist.Time[[AREAS[i]]]=cor(log(datos1$DaysAtLarge),log(datos1$dist.trav_m/1000))
      
      #.. bearing
      # mean(datos1$bearing, na.rm=FALSE)
      # sd(datos1$bearing, na.rm = FALSE)
      # quantile(datos1$bearing, probs = seq(0, 1, 0.25), na.rm = FALSE)
      # medianCircular(datos1$bearing, na.rm = FALSE)
      Bearing.summary[[AREAS[i]]]=summary(datos1$bearing)
      
      #.. condition effect on distance moved and time at large
      aa=sort(unique(datos1$CONDITION))
      Stats.all=list()
      for(j in 1:length(aa))
      {
        dat=subset(datos1,CONDITION==aa[j])
        stats.days=c(n=nrow(dat),summary(dat$DaysAtLarge))
        stats.dist=summary(dat$dist.trav_m/1000)
        Stats.all[[aa[j]]]=list(days=stats.days,dist=stats.dist)
      }
      STATS.COND[[AREAS[i]]]=Stats.all
    }
    return(list(Corr.Dist_time=COR.Dist.Time,Bearing.summary=Bearing.summary,Condition.summary=STATS.COND))
    
  }
  STATS=vector('list',length(SPECIES))
  names(STATS)=Species.names
  for(s in 1:length(SPECIES))
  {
    STATS[[s]]=function.stats(SPEC=SPECIES[s])
  }
  
  
  #4.4 Minimum time step for modelling movement rates
  #time spent at zones
  fn.time.zone=function(dat)
  {
    different.zone=subset(dat,!Areas==Areas.rec)
    different.zone$zone.number=with(different.zone,ifelse(Areas=="JANSF",1,
                                                          ifelse(Areas=="WANCSF",2,ifelse(Areas=="closed",3,ifelse(Areas=="WCDGDLF",4,
                                                                                                                   ifelse(Areas=="JASDGDLF.zone1",5,ifelse(Areas=="JASDGDLF.zone2",6,7)))))))
    
    
    different.zone$zone.rec.number=with(different.zone,ifelse(Areas.rec=="JANSF",1,
                                                              ifelse(Areas.rec=="WANCSF",2,ifelse(Areas.rec=="closed",3,ifelse(Areas.rec=="WCDGDLF",4,
                                                                                                                               ifelse(Areas.rec=="JASDGDLF.zone1",5,ifelse(Areas.rec=="JASDGDLF.zone2",6,7)))))))
    
    different.zone$delta.zone=abs(different.zone$zone.number-different.zone$zone.rec.number)
    
    diff=subset(different.zone,delta.zone>1)
    n.month=length(diff$DaysAtLarge[diff$DaysAtLarge<=30])
    cat(".....Number of indivdiuals moving to non-adjacent zones in less than a month=", n.month)
    
    n.yr=length(diff$DaysAtLarge[diff$DaysAtLarge<=365])
    cat("......Number of indivdiuals moving to non-adjacent zones in less than a year=", n.yr)
    
    cat("......Minimum number of days for moving to non-adjacent zones=", min(diff$DaysAtLarge))
    cat("*************************************************************************")
  }
  
  fn.time.zone(subset(Tagging,Species=="BW"))
  fn.time.zone(subset(Tagging,Species=="TK"))
  fn.time.zone(subset(Tagging,Species=="GM"))
  fn.time.zone(subset(Tagging,Species=="WH"))
  
  
  
  #5. GENERAL ANALYSES FOR MOVEMENT PATTERN PAPER
  
  #5.1 Table 2. Report all shark and ray species tagged and recaptured (numbers, distances, days)
  Numb.Sp.released=table(All.rel.Tagging$Species)
  Numb.Sp.recaptured=table(Tagging$Species)
  
  Rec.shks.rays=unique(Tagging$Species)
  fun.whatever=function(SPEC)
  {
    Dat=Tagging%>%
      filter(Species==SPEC)%>%
      mutate(dist.km=dist.trav_m/1000)
    
    dist.km=Dat%>%
      filter(!is.na(dist.km))%>%
      summarise(dist.Median=median(dist.km),
                dist.Up.95=quantile(dist.km,probs=0.975),
                dist.Max=max(dist.km))%>%
      mutate_all(funs(round(.)))
    
    Days=Dat%>%
      filter(!is.na(DaysAtLarge))%>%
      summarise(days.Median=median(DaysAtLarge),
                days.Up.95=quantile(DaysAtLarge,probs=0.975),
                days.Max=max(DaysAtLarge))%>%
      mutate_all(funs(round(.)))
    
    ROM=Dat%>%
      filter(!is.na(ROM))%>%
      summarise(ROM.Median=median(ROM),
                ROM.Up.95=quantile(ROM,probs=0.975),
                ROM.Max=max(ROM))%>%
      mutate_all(funs(round(.,3)))
    
    #max.speed.km.h=round(max(Dat$speed,na.rm=T)*3.6,3)
    
    prop.more.100=round(nrow(Dat%>%filter(dist.trav_m>=100000))/nrow(Dat),2)
    prop.more.500=round(nrow(Dat%>%filter(dist.trav_m>=500000))/nrow(Dat),2)
    
    dd=cbind(Species=Dat$COMMON_NAME[1],dist.km,Days,ROM,N=nrow(Dat),
             Prop.more.100=prop.more.100,Prop.more.500=prop.more.500)
    
    return(dd)
  }
  STORE.table1=vector('list',length=length(Rec.shks.rays))
  names(STORE.table1)=Rec.shks.rays
  for (i in 1:length(Rec.shks.rays)) STORE.table1[[i]]=fun.whatever(SPEC=Rec.shks.rays[i])
  Combined=do.call(rbind,STORE.table1)%>%
    arrange(-N)
  rownames(Combined)=NULL
  
  #add tag type and rel and rec fleet
  Dat=Tagging%>%
    filter(Species%in%Rec.shks.rays)%>%
    mutate(Tag.type=ifelse(Tag.type=='acoustic','conventional',Tag.type),
           Tag.type=ifelse(Tag.type=='conventional','Rototag',
                           ifelse(Tag.type=='conventional.dart','Dart tag',
                                  NA)),
           Species=COMMON_NAME)
  
  Tag.type=Dat%>%
    group_by(Species,Tag.type)%>%
    tally()%>%
    spread(Tag.type,n,fill='')%>%
    data.frame
  
  Rel.method=Dat%>%
    group_by(Species,Rel.method)%>%
    tally()%>%
    spread(Rel.method,n,fill='')%>%
    data.frame
  re.nm=names(Rel.method)[-1]
  Rel.method=Rel.method%>%
    rename_at(all_of(re.nm), ~ paste("Release",re.nm,sep='_'))
  
  Rec.method=Dat%>%
    group_by(Species,Rec.method)%>%
    tally()%>%
    spread(Rec.method,n,fill='')%>%
    data.frame
  re.nm=names(Rec.method)[-1]
  Rec.method=Rec.method%>%
    rename_at(all_of(re.nm), ~ paste("Recapture",re.nm,sep='_'))
  
  Add=full_join(Tag.type,Rel.method,by='Species')%>%
    full_join(Rec.method,by='Species')
  
  Combined=left_join(Combined,Add,by="Species")
  
  write.csv(Numb.Sp.released,file="Table_Numb.Sp.released.csv",row.names=F)
  write.csv(Numb.Sp.recaptured,file="Table_Numb.Sp.recaptured.csv",row.names=F)
  write.csv(Combined,file="Paper/Table1.csv",row.names=T) 
  tab_df(Combined,file="Paper/Table1.doc")
  
  
  # Observed distribution of displacement, ROM and time at liberty
  General.dist.displace=function(SPEC,maxDis)
  {
    datos=subset(Tagging,Species==SPEC)
    hist(datos$dist.trav_m/1000,breaks=breaks,col="gray",ylab="",main="",xlab="",
         xlim=c(0,maxDis),xaxt="n",cex.axis=1.5)
    #  hist(datos$dist.trav_m/1000,breaks=breaks,col="gray",ylab="",main="",xlab="",
    #       ylim=c(0,maxDisF),xlim=c(0,maxDis),xaxt="n",cex.axis=1.5)
    box()
    #  legend("topright",paste(Species.names," (n=",nrow(datos),")",sep=""),bty="n",cex=1.5)   
    
  }
  General.dist.time=function(SPEC)
  {
    datos=subset(Tagging,Species==SPEC)
    hist(datos$DaysAtLarge,breaks=breaks,col="gray",ylab="",main="",xlab="",
         xlim=c(0,breaks[length(breaks)]),xaxt="n",cex.axis=1.5)
    box()
  }
  General.speed=function(SPEC,Species.nm)
  {
    datos=subset(Tagging,Species==SPEC)
    hist(datos$ROM,breaks=breaks,col="gray",ylab="",main="",xlab="",
         xlim=c(0,breaks[length(breaks)]),xaxt="n",cex.axis=1.5)
    box()
    legend("topright",Species.nm,bty="n",cex=1.75)
  }
  General.dens.dis=function(SPECIES,LTY,COL)
  {
    datos=subset(Tagging,Species==SPECIES)
    lines(density(datos$dist.trav_m/1000,na.rm=T),lty=LTY,col=COL,lwd=2)
  }
  linetype=c(1:length(SPECIES))
  linecol=linetype
  maxDis=max(Tagging$dist.trav_m/1000,na.rm=T)
  maxDays=max(Tagging$DaysAtLarge,na.rm=T)
  maxDisFreq=c(100,300,100,100)
  maxLengFreq=c(120,80,40,40)
  maxY=c(0.0375,0.125,.075,.075)
  maxLeng=280
  BIN=20
  BIN.dist=100
  tiff(file="Paper/Observed_distribution_Dist.Days.ROM.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(length(SPECIES),3),oma=c(3.5,3.85,1.5,.5),mar=c(1,1,1,3),mgp=c(1, 0.75, 0),xpd=NA,las=2)
  #distance travelled
  breaks <- pretty(range(Tagging$dist.trav_m/1000,na.rm=T), 20)
  for(i in 1:length(SPECIES)) 
  {
    General.dist.displace(SPECIES[i],maxDis)
    axis(1,at=breaks,labels=F,tck=-0.025)
    axis(1,at=c(seq(breaks[1],breaks[length(breaks)],400)),labels=F,tck=-0.05)
    if(i==4)mtext("Displacement (km)",side=1,outer=F,line=2.5,font=1,las=0,cex=1.5)
  }
  axis(1,at=c(seq(breaks[1],breaks[length(breaks)],400)),labels=c(seq(breaks[1],breaks[length(breaks)],400))
       ,cex.axis=1.5,las=1,tck=-0.05)
  #days at liberty
  breaks <- pretty(range(Tagging$DaysAtLarge,na.rm=T), 20)
  for(i in 1:length(SPECIES)) 
  {
    General.dist.time(SPECIES[i])
    axis(1,at=breaks,labels=F,tck=-0.025)
    axis(1,at=c(seq(breaks[1],breaks[length(breaks)],1000)),labels=F,tck=-0.05)
    if(i==4)mtext("Time at liberty (days)",side=1,outer=F,line=2.5,font=1,las=0,cex=1.5)
  }
  axis(1,at=c(seq(breaks[1],breaks[length(breaks)],1000)),labels=c(seq(breaks[1],
                                                                       breaks[length(breaks)],1000)),cex.axis=1.5,las=1,tck=-0.05)
  #ROM
  breaks <- pretty(range(Tagging$ROM,na.rm=T), 20)
  for(i in 1:length(SPECIES)) 
  {
    General.speed(SPECIES[i],Species.names[i])
    axis(1,at=breaks,labels=F,tck=-0.025)
    axis(1,at=breaks[c(1,6,11,16,21)],labels=F,cex.axis=1.5,las=1,tck=-0.05)
    if(i==4)mtext("ROM (km/d)",side=1,outer=F,line=2.5,font=1,las=0,cex=1.5)
  }
  axis(1,at=breaks[c(1,6,11,16,21)],labels=breaks[c(1,6,11,16,21)],cex.axis=1.5,las=1,tck=-0.05)
  mtext("Frequency",side=2,outer=T,line=2.1,font=1,las=0,cex=1.75)
  dev.off()
  
  
  
  # Observed bearing distribution by area
  rose1 <- function (frec, fnum = 4, fint = 5, flab = 2, ang = 3 * pi/16,               
                     col = rainbow(10, 0.5, 0.92, start = 0.33, end = 0.2),
                     margen = c(0,0, .75, 0), key = TRUE, uni = "m/s", plotnew=F, ...) 
  {
    #old.par <- par(no.readonly = TRUE)
    #on.exit(par(old.par))
    if (is.matrix(frec)) 
      frec <- as.data.frame(frec)
    if (is.vector(frec)) {
      ndir <- length(frec)
      nr <- 1
    }
    else {
      ndir <- length(frec[1, ])
      nr <- nrow(frec)
    }
    fmax <- fnum * fint
    tot <- sum(frec)
    fr <- 100 * frec/tot
    key <- (nr > 1) && key
    if (key) 
      mlf <- 3
    else mlf <- 1
    par(mar = margen)#, new = plotnew)
    fx <- cos(pi/2 - (2 * pi/ndir * 0:(ndir - 1)))
    fy <- sin(pi/2 - (2 * pi/ndir * 0:(ndir - 1)))
    plot(fx, fy, xlim = c(-fmax - mlf * fint, fmax + fint), ylim = c(-fmax - 
                                                                       fint, fmax + fint), xaxt = "n", yaxt = "n", xlab = "", 
         ylab = "", bty = "n", asp = 1, type = "n", cex.main=1.5, ...)
    if (nr == 1) {
      cx <- fx * fr
      cy <- fy * fr
    }
    else {
      f <- apply(fr, 2, sum)
      cx <- fx * f
      cy <- fy * f
      for (i in nr:2) {
        f <- f - fr[i, ]
        cx <- c(cx, NA, fx * f)
        cy <- c(cy, NA, fy * f)
      }
    }
    polygon(cx, cy, col = col[nr:1])
    symbols(c(0 * 1:fnum), c(0 * 1:fnum), circles = c(fint * 
                                                        1:fnum), inches = FALSE, add = TRUE)
    segments(0 * 1:ndir, 0 * 1:ndir, fmax * fx, fmax * fy)
    fmaxi <- fmax + fint/4
    text(0, fmaxi*1.05, "N")
    text(0, -fmaxi*1.05, "S")
    text(fmaxi*1.05, 0, "E")
    text(-fmaxi*1.05, 0, "W")
    if (flab == 2) 
      for (i in 1:fnum) text(i * fint * cos(ang), i * fint * 
                               sin(ang), paste(i * fint, "%"))
    else if (flab == 1) 
      text(fmax * cos(ang), fmax * sin(ang), paste(fmax, "%"))
    if (key) {
      #legend(-fmaxi + 10 * fint, fmaxi + 2, fill = col, legend = attr(frec,"row.names"),bty='n')
      
      text(-fmaxi - 1.4 * fint, fmaxi + 0.9 * fint, uni)
    }
    invisible()
  }
  Rosa.vent.bearing1=function(AREAS,FNUM,FINT,Fig6.areas)
  {
    datos=subset(Tagging,Areas%in%AREAS)
    unicaSp=as.character(unique(datos$Species))
    unicasSp=subset(unicaSp,unicaSp%in%SPECIES)   
    unicasSp=unicasSp[match(c("TK","BW","GM","WH"),unicasSp)]
    Sp.names=Species.names[!(is.na(match(names(Species.names),unicasSp)))]
    NN=length(unicasSp)
    Bear.coord=as.data.frame(matrix(NA,nrow=NN,ncol=16))
    colnames(Bear.coord)=names.Bear.coord
    #rownames(Bear.coord)=Sp.names
    for (i in 1:NN)
    {
      datos2=subset(datos,Species==unicasSp[i])
      a=hist(datos2$bearing,plot=F,breaks=BEARINGS)
      shifted=a$counts
      add=shifted[1]+shifted[length(shifted)]
      shifted=c(add,shifted[-c(1,length(shifted))])
      Bear.coord[i,]=shifted
    }  
    if(AREAS[1]=="JASDGDLF.zone2")
    {
      
      #rosavent(Bear.coord, FNUM,FINT, ang=-3*pi/16, main=Species.names,key =T,uni="",flab=F
      #         ,col=rosa.colors)    
      rose1(Bear.coord, FNUM,FINT, ang=-3*pi/16, main=Fig6.areas,key =F,
            uni="",flab=F,col=rosa.colors)
    }else
    {
      #rosavent(Bear.coord, FNUM,FINT, ang=-3*pi/16, main=Species.names,key = F,uni="",flab=F
      #         ,col=rosa.colors)
      rose1(Bear.coord, FNUM,FINT, ang=-3*pi/16, main=Fig6.areas,key = F,
            uni="",flab=F,col=rosa.colors)
    }
  }
  CLO=c("black","grey95","grey50","grey75")
  names(CLO)=c("TK","BW","GM","WH")
  AREAS.bearing=list("Closed","WCDGDLF","JASDGDLF.zone1","JASDGDLF.zone2")
  BEARINGS=c(0,seq(0+(22.5/2),360-(22.5/2),by=360/16),360)
  names.Bear.coord=c("N","NNE","NE","ENE","E","ESE","SE","SSE", "S","SSW","SW","WSW","W","WNW","NW","NNW")
  Label.areas=list("Closed Area","WCDGDLF","Zone1 JASDGDLF","Zone2 JASDGDLF")
  Fnum=c(5,6,6,6)
  Fint=c(5,4,3.7,3.2)
  rosa.colors=CLO
  tiff(file="Paper/Figure5.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
  par(mfcol=c(2,2),oma=c(.1,.1,.1,.1),mar=c(.1,.1,.1,.2),mgp=c(1, 0.75, 0),xpd=NA,las=2)
  for(x in 1:length(AREAS.bearing))  Rosa.vent.bearing1(AREAS.bearing[[x]],Fnum[x],Fint[x],Label.areas[[x]])
  par(fig=c(.35,.625,.35,.65), new = T,mgp=c(.1,.4,0))
  plot(1:10,1:10,bty='n',pch='',xaxt='n',yaxt='n',ann=F)
  legend("center",Species.names,fill=rosa.colors,bty='n',cex=1.5)
  dev.off()
  
  
  
  # Root mean squared displacement vs time lag
  do.rmsd=FALSE
  if(do.rmsd)
  {
    #Calculate weighted mean squared displacement and root mean squared displacement (Whitehead et al 2008)
    
    #Jackkife SE
    MSD.fn=function(SPECIES)
    {
      datos=subset(Tagging,Species==SPECIES)
      
      #calculate quantities
      RMSD=SE.RMSD=NULL
      for(i in 1:length(LAG))
      {
        xdata=subset(datos,DaysAtLarge>=LAG[i])    
        n <- nrow(xdata)
        id=match("speed",names(xdata))
        
        theta <- function(x,xdata)
        { 
          displacement=LAG[i]*60*60*24*xdata[x,id]/1000 #(in km)
          squared.dis=displacement^2
          #WEIGHT=datos2$WEIGHT    #weight by effort
          #WEIGHT=datos2$DaysAtLarge #weight by days at large
          #mean.squared.dis=weighted.mean(squared.dis,WEIGHT, na.rm = T)
          mean.squared.dis=mean(squared.dis,na.rm = T)
          root.mean.squared.dis=(mean.squared.dis)^0.5
          return(root.mean.squared.dis)
        }
        results <- jackknife(1:n,theta,xdata)
        SE=results$jack.se
        Mean=mean(results$jack.values)  
        
        SE.RMSD=cbind(SE.RMSD,SE)
        RMSD=cbind(RMSD,Mean)
      }
      return(list(mean=RMSD,se=SE.RMSD))
    }
    
    #Bootstrapping SE
    # MSD.fn=function(SPECIES)
    #   {
    #     datos=subset(Tagging,Species==SPECIES)
    #     
    #     #calculate quantities
    #     RMSD=N=NULL
    #     for(i in 1:length(LAG))
    #     {
    #       datos1=subset(datos,DaysAtLarge>=LAG[i])
    #       Summary.boot=NULL
    #       for(j in 1:nboot)
    #       {
    #         datos2=datos1[sample(1:nrow(datos1), nrow(datos1), replace = TRUE),] #random sample with replacement
    # 
    #         displacement=LAG[i]*60*60*24*datos2$speed/1000 #(in km)
    #         squared.dis=displacement^2
    #         #WEIGHT=datos2$WEIGHT    #weight by effort
    #         #WEIGHT=datos2$DaysAtLarge #weight by days at large
    #         #mean.squared.dis=weighted.mean(squared.dis,WEIGHT, na.rm = T)
    #         mean.squared.dis=mean(squared.dis,na.rm = T)
    #         root.mean.squared.dis=(mean.squared.dis)^0.5
    #         Summary.boot=rbind(Summary.boot,root.mean.squared.dis)
    #       }
    #       N=cbind(N,nrow(datos1))
    #       RMSD=cbind(RMSD,Summary.boot)
    #     }
    #     MEAN=apply(RMSD,2,mean)
    #     SD=apply(RMSD,2,sd)
    #     SE=SD/(N)^.5
    #     low.CI=apply(RMSD,2,function(x){quantile(x,prob=.025)})
    #     up.CI=apply(RMSD,2,function(x){quantile(x,prob=.975)})
    #     
    #     return(list(mean=MEAN,se=SD,low.CI=low.CI,up.CI=up.CI))
    #   }
    
    
    #5. MOVEMENT DESCRIPTORS
    
    #5.1 Figure 4.
    col.species=rep(1,4)
    #col.species=c("grey50","black","grey30","grey50")
    pch.species=rep(21,4)
    bg.species=rosa.colors
    LAG=c(1,7,30,182.5,365,730,1095) #time lags (in days) 
    Time.Lag.legends=c("1 d","1 wk","1 mo","6 mo","1 yr","2 yr","3 yr")
    
    nboot=1000    #number of bootstraps
    dummy.MEAN=c(3,12,44,173,200,350,500)
    MAX.Y=1000
    tiff(file="Root mean squared displacement.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(1,1),oma=c(2,2.5,.1,.1),mar=c(2,2,.1,.1),las=2,mgp=c(2.5, 1, 0),xpd=T)
    plot(log(LAG),MSD.fn(SPECIES[1])$mean,type='p',ylab="",xlab="", yaxt='n',xaxt='n',cex.lab=1.5,
         col='transparent',ylim=c(0,MAX.Y))
    for (i in 1:length(SPECIES))
    {
      MEAN=MSD.fn(SPECIES[i])$mean
      SE=MSD.fn(SPECIES[i])$se
      #UP=MSD.fn(SPECIES[i])$up.CI
      #LOW=MSD.fn(SPECIES[i])$low.CI
      points(log(LAG),MEAN,pch=pch.species[i],col=col.species[i],cex=2,bg=bg.species[i])
      arrows(log(LAG),MEAN+SE, log(LAG), MEAN-SE, angle=90, code=3,length=.05,
             col=col.species[i],lwd=2)
      #    arrows(log(LAG),LOW, log(LAG), UP, angle=90, code=3,length=.05,
      #      col=col.species[i],lwd=2)
    }
    axis(side = 1, at = log(LAG), labels = Time.Lag.legends,tck=-0.02,las=1,cex.axis=1.2)
    axis(side = 2, at =seq(0,MAX.Y,100), labels = seq(0,MAX.Y,100),tck=-0.02,las=1,cex.axis=1.2)
    legend("topleft",Species.names,bty="n",cex=1.5,col=col.species,pch=pch.species,pt.bg=bg.species)
    mtext("RMS displacement (km)",side=2,outer=T,line=1.25,font=1,las=0,cex=1.75)
    mtext("Time lag (log scale)",side=1,outer=T,line=0.5,font=1,las=0,cex=1.75)
    text(7.07,-87,"3 yr",cex=1.21)
    dev.off()
  }  
  
  
  
  # LINEAR MODELLING FOR MOVEMENT PATTERN PAPER
  
  
  #Explore
  library(ggpubr)
  fun.explr=function(d,bin)
  {
    d=d%>%
      mutate(dist.km=dist.trav_m/1000,
             dist.km.1.yr=dist.trav_m_in.1.year/1000,
             month=factor(Mn.rel,levels=1:12),
             SP=factor(Species,levels=c("WH","GM","BW","TK")))%>%
      filter(!is.na(Rel_FL) & !is.na(Sex))%>%
      mutate(Rel_FL.bin=bin*(floor(Rel_FL/bin)),
             Rel_FL.bin=factor(Rel_FL.bin,levels=seq(min(Rel_FL.bin),max(Rel_FL.bin),bin)))
    
    a=d%>%
      ggplot(aes(x=log(dist.trav_m), fill=COMMON_NAME)) + 
      geom_histogram()+ 
      theme(legend.position="top",legend.title=element_blank())
    
    b=d%>%
      ggplot(aes(x=log(ROM), fill=COMMON_NAME)) + 
      geom_histogram(show.legend = FALSE)+ 
      theme(legend.position="top",legend.title=element_blank())
    
    p=d%>%
      ggplot(aes(x=Rel_FL, y=dist.km, colour=COMMON_NAME)) + 
      geom_point()+  geom_smooth(method = "loess") +
      theme(legend.position="top",legend.title=element_blank())
    
    p0=d%>%
      ggplot(aes(x=Rel_FL, y=ROM, colour=COMMON_NAME)) + 
      geom_point(show.legend = FALSE)+  geom_smooth(method = "loess") +
      theme(legend.position="top",legend.title=element_blank())
    
    p1=d%>%
      ggplot(aes(x=Rel_FL.bin, y=dist.km, fill=COMMON_NAME)) + 
      geom_boxplot()+ theme(legend.position="top",legend.title=element_blank())
    p2=d%>%
      ggplot(aes(x=Rel_FL.bin, y=ROM, fill=COMMON_NAME)) + 
      geom_boxplot(show.legend = FALSE)
    
    p3=d%>%
      ggplot(aes(x=month, y=dist.km, fill=COMMON_NAME)) + 
      geom_boxplot()+ theme(legend.position="top",legend.title=element_blank())
    p4=d%>%
      ggplot(aes(x=month, y=ROM, fill=COMMON_NAME)) + 
      geom_boxplot(show.legend = FALSE)
    
    p5=d%>%
      ggplot(aes(x=Areas, y=dist.km, fill=COMMON_NAME)) + 
      geom_boxplot()+ theme(legend.position="top",legend.title=element_blank())
    p6=d%>%
      ggplot(aes(x=Areas, y=ROM, fill=COMMON_NAME)) + 
      geom_boxplot(show.legend = FALSE)
    
    p7=d%>%
      ggplot(aes(x=Sex, y=dist.km, fill=COMMON_NAME)) + 
      geom_boxplot()+ theme(legend.position="top",legend.title=element_blank())
    p8=d%>%
      ggplot(aes(x=Sex, y=ROM, fill=COMMON_NAME)) + 
      geom_boxplot(show.legend = FALSE)
    
    #density distribution by month
    # ggplot(d, aes(x = ROM, y = factor(Mn.rel), fill = ..x..)) +
    #   geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
    #   scale_fill_viridis(name = "Temp. [F]", option = "C") +
    #   coord_cartesian(xlim = c(0, quantile(d$ROM,0.99))) +
    #   theme_ipsum() + labs(title = 'Month released')+
    #   xlab("ROM (km per day)")+ylab('Month')+
    #   theme(legend.position="none",panel.spacing = unit(0.1, "lines"),strip.text.x = element_text(size = 8))+
    #   facet_wrap(~COMMON_NAME)
    
    
    plot_list=list(a,b,p,p0,p1,p2,p3,p4,p5,p6,p7,p8)
    multi.page <-ggarrange(plotlist=plot_list, nrow = 2, ncol = 1)
    ggexport(multi.page, filename = "explore_linear_model.pdf")
    
  }
  fun.explr(d=Tagging%>%filter(Species%in%SPECIES),bin=10)
  
  #Fit model
  library(ciTools)
  library(multcomp)
  
  fit.mod=function(bin=10)
  {
    d=Tagging%>%
      filter(Species%in%SPECIES & Sex%in%c("F","M"))%>%
      filter(!is.na(Rel_FL))%>%
      mutate(SP=factor(Species),
             Sex=as.factor(Sex),
             CONDITION=as.factor(CONDITION),
             Log.DaysAtLarge=log(DaysAtLarge),
             Log.ROM=log(ROM),
             Log.dist=log(dist.trav_m),
             Log.Rel_FL=log(Rel_FL),
             Areas=as.factor(as.character(Areas)),
             Rel_FL.bin=bin*(floor(Rel_FL/bin)),
             Rel_FL.bin=factor(Rel_FL.bin,levels=seq(min(Rel_FL.bin),max(Rel_FL.bin),bin)))
    
    
    mod.dist<-glm(Log.dist~ Log.DaysAtLarge+SP+Log.Rel_FL+Areas+Sex+CONDITION,
                  family=gaussian(link = "identity"),data=d) 
    
    #summary(glht(mod.dist, mcp(SP="Tukey")))
    
    Sig.Distance.km.glm=anova(mod.dist,test='Chisq')	
    Dev.exp=round(Dsquared(mod.dist,adjust=F)$d3)
    
    NewData=expand.grid(Log.DaysAtLarge=log(mean(d$DaysAtLarge)),
                        SP=factor(levels(d$SP),levels=levels(d$SP)),
                        Areas=factor("JASDGDLF.zone1",levels=levels(d$Areas)),
                        Log.Rel_FL=log(seq(min(10*floor(d$Rel_FL/10)),max(d$Rel_FL),by=10)),
                        Sex=factor("F",levels=levels(d$Sex)),
                        CONDITION=factor("1",levels=levels(d$CONDITION)))
    NewData=add_ci(NewData,mod.dist,alpha=0.05,response=TRUE)
    
    NewData%>%
      ggplot(aes(x=exp(Log.Rel_FL),y=exp(pred)/1000,colour=SP))+geom_line()+
      geom_smooth(aes(ymin = exp(LCB0.025)/1000, 
                      ymax = exp(UCB0.975)/1000),
                  stat = "identity")
    
    geom_errorbar(aes(ymin=exp(LCB0.025)/1000, ymax=exp(UCB0.975)/1000), width=.1)
    
  }
  
  
  #select only species of interest
  Tagging.mySpecies=subset(Tagging,Species%in%SPECIES)
  Tagging.mySpecies=subset(Tagging.mySpecies,Sex%in%c("F","M"))
  
  
  #--Distance travelled
  #Plot fit diagnostics function
  fn.plot.diag=function(MODEL,TEXT,SPECIES)
  {
    RES=MODEL$residuals   #residuals
    Std.RES=RES/sd(RES)   #standardised residuals (res/SD(res))
    #Std.RES=rstandard(MODEL)
    
    par(mfcol=c(2,2),las=1,mar=c(3,3,2,1),oma=c(2.5,.1,.1,.1),las=1,mgp=c(2,.5,0),cex.axis=.8,cex.lab=1.1)
    qqnorm(RES,main="",ylim=c(-5,5),xlim=c(-5,5),ylab="Residuals",xlab="Quantiles of standard normal distribution")
    qqline(RES, col = 'grey40',lwd=1.5,lty=2)
    
    hist(Std.RES,xlim=c(-5,5),ylab="Frequency",xlab="Stan. residuals",main="",col="grey",breaks=50)
    box()
    
    plot(predict(MODEL),Std.RES,ylim=c(-4,12),ylab="Stan. residuals",xlab="Expected values")
    abline(0,0,lwd=1.5,lty=2,col='grey40')
    
    plot(predict(MODEL),sqrt(abs(Std.RES)),ylim=c(0,2.6),ylab="Square root of stan. residuals",xlab="Expected values")
    
    mtext(paste(TEXT,SPECIES),3,outer=T,lin=-1)
  }
  
  #Deviance explained function
  Dsquared <- function(model, adjust = F)
  {
    if(!is.null(model$deviance))
    {
      d2 <- (model$null.deviance - model$deviance) / model$null.deviance
      if (adjust)
      {
        n <- length(model$fitted.values)
        p <- length(model$coefficients)
        d2 <- 1 - ((n - 1) / (n - p)) * (1 - d2)
      }
      d3=d2*100
      d1=model$deviance
      d2=model$null.deviance
    }
    
    if(is.null(model$deviance))
    {
      d1=NA
      d2=NA
      d3=NA
    }
    return(list(d1=d1,d2=d2,d3=d3))
  }
  term.dev.explained=function(MOD)
  {
    Anova.tab=anova(MOD, test = "Chisq")
    n=2:length(Anova.tab$Deviance)
    Term.dev.exp=100*(Anova.tab$Deviance[n]/MOD$null.deviance)
    names(Term.dev.exp)=rownames(Anova.tab)[n]
    Dev.exp=sum(Term.dev.exp)
    ANOVA=as.data.frame.matrix(Anova.tab)
    Term=data.frame(Percent.dev.exp=Term.dev.exp)
    Table=ANOVA[-1,match(c("Df","Pr(>Chi)"),names(ANOVA))]
    Term=Term[match(rownames(Term), rownames(Table)),]
    Table=cbind(Table,Term)
    names(Table)[match("Term",names(Table))]="Percent.dev.exp"
    Table$term=rownames(ANOVA)[2:nrow(ANOVA)]
    Table=Table%>%select('term','Df','Pr(>Chi)','Percent.dev.exp')
    All=Table[1,]
    All[,1:ncol(All)]=NA
    All$term="model"
    All$Percent.dev.exp=round(Dev.exp,3)
    Table=rbind(Table,All)
    vars <- c(df = "Df", 'p-value' ="Pr(>Chi)")
    Table= Table %>% mutate_at(c("Pr(>Chi)","Percent.dev.exp"), round, 3) %>%
      rename(!!vars)
  }
  
  Species.size=FALSE #prelim anal showed non significant interaction
  
  #1. Displacement 
  do.displacement=TRUE
  if(do.displacement)
  {
    Dat=subset(Tagging.mySpecies,dist.trav_m>1 & Rel_FL>0 & CONDITION%in%c("1","2","3"))
    Dat$Species=as.factor(Dat$Species)
    Dat$Sex=as.factor(Dat$Sex)
    Dat$CONDITION=as.factor(Dat$CONDITION)
    Dat$Log.DaysAtLarge=log(Dat$DaysAtLarge)
    Dat$Log.Rel_FL=log(Dat$Rel_FL)
    Dat$Areas=as.character(Dat$Areas)
    Dat$Areas=as.factor(Dat$Areas)
    
    
    
    if(!isTRUE(Species.size))Distance.km.glm<-glm(log(dist.trav_m/1000)~ Log.DaysAtLarge+Species+Areas+Log.Rel_FL+Sex+CONDITION,
                                                  family=gaussian(link = "identity"),data=Dat) 
    if(Species.size) Distance.km.glm<-glm(log(dist.trav_m/1000)~ Log.DaysAtLarge+Species*Log.Rel_FL+Areas+Sex+CONDITION,
                                          family=gaussian(link = "identity"),data=Dat) 
    
    Post.hoc.comp.dist=summary(glht(Distance.km.glm, mcp(Species="Tukey")))
    
    
    #diagnostic plots (fit check)
    tiff(file="Paper/GLM.diagnostics.displacement.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    fn.plot.diag(Distance.km.glm,"Lognormal error ","Displacement km")
    dev.off()
    
    stats.Distance.km.glm=summary(Distance.km.glm)            
    
    #terms signficance    
    Sig.Distance.km.glm=anova(Distance.km.glm,test='Chisq')		
    
    BiasCor.fn=function(Median,SE) biasCorr <- exp(Median+(SE^2)/2) #function for bias corrected mean in normal space
    
    
    #combine anova and coefficients for publication
    Coefs=as.data.frame(stats.Distance.km.glm$coefficients)
    Coefs$Terms=rownames(Coefs)
    Just.coef=Coefs[,1] 
    names(Just.coef)=rownames(Coefs)
    Just.coef[c(2,11)]
    ANOVA=term.dev.explained(MOD=Distance.km.glm)
    write.csv(ANOVA,"Paper/ANOVA.displacement.csv",row.names=F)

    #Predict significant terms
    if(!isTRUE(Species.size))
    {
      NewData=expand.grid(Log.DaysAtLarge=log(mean(Dat$DaysAtLarge)),
                          Species=factor(SPECIES,levels=levels(Dat$Species)),
                          Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
                          Log.Rel_FL=log(mean(Dat$Rel_FL)),Sex=factor("F",levels=levels(Dat$Sex)),
                          CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    }
    if(Species.size)
    {
      NewData=expand.grid(Log.DaysAtLarge=log(mean(Dat$DaysAtLarge)),
                          Species=factor(SPECIES,levels=levels(Dat$Species)),
                          Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
                          Log.Rel_FL=log(seq(min(10*floor(Dat$Rel_FL/10)),max(Dat$Rel_FL),by=10)),
                          Sex=factor("F",levels=levels(Dat$Sex)),
                          CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    }
    NewData.dist=NewData
    
    NewData.dist=add_ci(NewData.dist,Distance.km.glm,alpha=0.05,response=TRUE)
    
    Pred.dist.travel=predict(Distance.km.glm,newdata=NewData,type="response",se.fit=T)
    #Pred.dist.travel$fit=BiasCor.fn(Pred.dist.travel$fit,Pred.dist.travel$se.fit)
    
    #   #FL-Species
    # FL.range.pred=80:160
    # Pred_FL_SP.dist.trav=expand.grid(Log.DaysAtLarge=log(mean(Dat$DaysAtLarge)),
    #                     Species=factor(SPECIES,levels=levels(Dat$Species)),
    #                     Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
    #                     Log.Rel_FL=log(FL.range.pred),Sex=factor("F",levels=levels(Dat$Sex)),
    #                     CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    # ss=predict(Distance.km.glm,newdata=Pred_FL_SP.dist.trav,type="response",se.fit=T)
    # Pred_FL_SP.dist.trav$pred=ss$fit
    # Pred_FL_SP.dist.trav$pred.se=ss$se.fit
  }
  
  
  #2. ROM 
  do.ROM=TRUE
  if(do.ROM)
  {
    Dat=subset(Tagging.mySpecies,speed>0 & Rel_FL>0 & CONDITION%in%c("1","2","3"))
    Dat$Species=as.factor(Dat$Species)
    Dat$Sex=as.factor(Dat$Sex)
    Dat$CONDITION=as.factor(Dat$CONDITION)
    Dat$Log.DaysAtLarge=log(Dat$DaysAtLarge)
    Dat$Log.Rel_FL=log(Dat$Rel_FL)
    Dat$Areas=as.character(Dat$Areas)
    Dat$Areas=as.factor(Dat$Areas)
    
    if(!isTRUE(Species.size)) Speed.glm<-glm(log(ROM)~ Species+Areas+Log.Rel_FL+Sex+CONDITION,
                                             family=gaussian(link = "identity"),data=Dat) 
    if(Species.size) Speed.glm<-glm(log(ROM)~ Species*Log.Rel_FL+Areas+Sex+CONDITION,
                                    family=gaussian(link = "identity"),data=Dat) 
    
    
    Post.doc.comp.ROM=summary(glht(Speed.glm, mcp(Species="Tukey")))
    
    #diagnostic plots (fit check)
    tiff(file="Paper/GLM.diagnostics.Speed.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    fn.plot.diag(Speed.glm,"Lognormal error ","ROM")
    dev.off()
    
    stats.speed.glm=summary(Speed.glm)            
    
    #terms signficance    
    Sig.speed.glm=anova(Speed.glm,test='Chisq')  	
    
    
    #combine anova and coefficients for publication
    Coefs=as.data.frame(stats.speed.glm$coefficients)
    Coefs$Terms=rownames(Coefs)
    Just.coef=Coefs[,1] 
    names(Just.coef)=rownames(Coefs)
    Just.coef[10]
    ANOVA=term.dev.explained(MOD=Speed.glm)
    write.csv(ANOVA,"Paper/ANOVA.speed.csv",row.names=F)
    
  
    #Predict significant terms
    if(!isTRUE(Species.size))
    {
      NewData=expand.grid(Species=factor(SPECIES,levels=levels(Dat$Species)),
                          Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
                          Log.Rel_FL=log(mean(Dat$Rel_FL)),Sex=factor("F",levels=levels(Dat$Sex)),
                          CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    }
    if(Species.size)
    {
      NewData=expand.grid(Species=factor(SPECIES,levels=levels(Dat$Species)),
                          Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
                          Log.Rel_FL=log(seq(min(10*floor(Dat$Rel_FL/10)),max(Dat$Rel_FL),by=10)),
                          Sex=factor("F",levels=levels(Dat$Sex)),
                          CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    }
    NewData.rom=NewData
    NewData.rom=add_ci(NewData.rom,Speed.glm,alpha=0.05,response=TRUE)
    
    Pred.speed=predict(Speed.glm,newdata=NewData,type="response",se.fit=T)
    #Pred.speed$fit=BiasCor.fn(Pred.speed$fit,Pred.speed$se.fit)
    
  }
  
  
  #3. Days at large GLM
  do.days.large=FALSE
  if(do.days.large)
  {
    Dat=subset(Tagging.mySpecies,DaysAtLarge>0 & Rel_FL>0 & CONDITION%in%c("1","2","3"))
    Dat$Species=as.factor(Dat$Species)
    Dat$Sex=as.factor(Dat$Sex)
    Dat$CONDITION=as.factor(Dat$CONDITION)
    Dat$Log.DaysAtLarge=log(Dat$DaysAtLarge)
    Dat$Log.Rel_FL=log(Dat$Rel_FL)
    Dat$Areas=as.character(Dat$Areas)
    Dat$Areas=as.factor(Dat$Areas)
    
    Days.glm.Pois<-glm(DaysAtLarge~ Species+Areas+Log.Rel_FL+Sex+CONDITION,
                       family=poisson(link = "log"),data=Dat) 
    
    Days.glm.NB<-glm.nb(DaysAtLarge~ Species*Areas+Log.Rel_FL+Sex+CONDITION,data=Dat) 
    
    
    #note: If not >> 1, then Poisson...
    OverDisp=function(MoD)
    {
      E1=resid(MoD,type='pearson')
      Dispersion=sum(E1^2)/MoD$df.resid 
      if(is.null(MoD))Dispersion=NA
      return(list(Dispersion=Dispersion))
    }
    OverDisp(Days.glm.Pois)
    OverDisp(Days.glm.NB)
    
    #diagnostic plots (fit check)
    tiff(file="Paper/Days.GLM.diagnostics.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
    fn.plot.diag(Days.glm.NB,"Negative Binomial error ","Days at large")
    dev.off()
    
    stats.days.glm=summary(Days.glm.NB)            
    
    #terms signficance    
    Sig.days.glm=anova(Days.glm.NB,test='Chisq')    
    
    
    #combine anova and coefficients for publication
    Coefs=as.data.frame(stats.days.glm$coefficients)
    Coefs$Terms=rownames(Coefs)
    Just.coef=Coefs[,1] 
    names(Just.coef)=rownames(Coefs)
    Just.coef[10]
    ANOVA=as.data.frame(Sig.days.glm)
    ANOVA=ANOVA[-1,]
    ANOVA$Terms=rownames(ANOVA)
    rownames(ANOVA)=NULL
    ANOVA=ANOVA[,match(c("Terms","Df","Deviance","Pr(>Chi)"),names(ANOVA))]
    ANOVA$"Pr(>Chi)"=round(ANOVA$"Pr(>Chi)",4)
    write.csv(ANOVA,"Paper/days.ANOVA.csv",row.names=F)
    
    #deviance explained
    Dev.exp=round(Dsquared(Days.glm.NB,adjust=F)$d3)
    write.csv(Dev.exp,"Paper/days.Dev.exp.csv",row.names=F) 
    
    
    #Predict significant terms
    #Species
    NewData=expand.grid(Species=factor(SPECIES,levels=levels(Dat$Species)),
                        Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
                        Log.Rel_FL=log(mean(Dat$Rel_FL)),Sex=factor("F",levels=levels(Dat$Sex)),
                        CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    
    
    Pred.days=predict(Days.glm.NB,newdata=NewData,type="response",se.fit=T)
    
    
    # Pred_FL_SP.days=expand.grid(Species=factor(SPECIES,levels=levels(Dat$Species)),
    #                              Areas=factor("JASDGDLF.zone1",levels=levels(Dat$Areas)),
    #                              Log.Rel_FL=log(FL.range.pred),Sex=factor("F",levels=levels(Dat$Sex)),
    #                              CONDITION=factor("1",levels=levels(Dat$CONDITION)))
    # ss=predict(Days.glm.NB,newdata=Pred_FL_SP.days,type="response",se.fit=T)
    # Pred_FL_SP.days$pred=ss$fit
    # Pred_FL_SP.days$pred.se=ss$se.fit
    
  }
  
  
  #plot predictions
  fn.bar.inter=function(what,what.se1,what.se2,CLO,YLIM,BY)
  {
    what=rbind(c(what[5:7],0),what[1:4],c(0,what[8:10]))
    what=t(what)
    
    what.se1=rbind(c(what.se1[5:7],0),what.se1[1:4],c(0,what.se1[8:10]))
    what.se1=t(what.se1)
    
    
    what.se2=rbind(c(what.se2[5:7],0),what.se2[1:4],c(0,what.se2[8:10]))
    what.se2=t(what.se2)
    
    
    
    mp <- barplot(what,beside=T, axes=FALSE, axisnames=FALSE, ylim=YLIM,col=CLO, main="", xlab="", ylab="")
    for(a in 1:ncol(mp))
    {
      arrows(mp[,a], what[,a], mp[,a] , what.se1[,a],lwd=2,angle=90,length = 0.1)
      arrows(mp[,a], what[,a], mp[,a] , what.se2[,a],lwd=2,angle=90,length = 0.1)
    }
    
    axis(2, at=seq(YLIM[1], YLIM[2], by=BY),cex.axis=1.5)
    
    axis(1, labels=c("WCDGDLF","JASDGDLF (zone1)","JASDGDLF (zone2)"),
         at = c(mean(mp[,1]),mean(mp[,2]),mean(mp[,3])),cex.axis=1.75)
    
    
    box()
  }
  
  fn.bar.inter.st=function(what,what.se1,what.se2,CLO,YLIM,BY)
  {
    what=rbind(c(what[5:7],0),what[1:4],c(0,what[8:10]))
    what=what/rowMeans(what)   #standardise to mean of 1
    what=t(what)
    
    what.se1=rbind(c(what.se1[5:7],0),what.se1[1:4],c(0,what.se1[8:10]))
    what.se1=what.se1/rowMeans(what.se1)
    what.se1=t(what.se1)
    
    
    what.se2=rbind(c(what.se2[5:7],0),what.se2[1:4],c(0,what.se2[8:10]))
    what.se2=what.se2/rowMeans(what.se2)
    what.se2=t(what.se2)
    
    
    
    mp <- barplot(what,beside=T, axes=FALSE, axisnames=FALSE, ylim=YLIM,col=CLO, main="", xlab="", ylab="")
    #   for(a in 1:ncol(mp))
    #   {
    #     arrows(mp[,a], what[,a], mp[,a] , what.se1[,a],lwd=2,angle=90,length = 0.1)
    #     arrows(mp[,a], what[,a], mp[,a] , what.se2[,a],lwd=2,angle=90,length = 0.1)
    #   }
    
    axis(2, at=seq(YLIM[1], YLIM[2], by=BY),cex.axis=1.5)
    
    axis(1, labels=c("WCDGDLF","JASDGDLF (zone1)","JASDGDLF (zone2)"),
         at = c(mean(mp[,1]),mean(mp[,2]),mean(mp[,3])),cex.axis=1.75)
    
    
    box()
  }
  
  fn.bar=function(what,what.se1,what.se2,CLO,YLIM,DO.axis,BY)
  {
    mp <- barplot(what, axes=FALSE, axisnames=FALSE, ylim=YLIM,col=CLO, main="", xlab="", ylab="")
    #arrows(mp, what, mp , what.se1,lwd=2,angle=90,length = 0.1)
    #arrows(mp, what, mp , what.se2,lwd=2,angle=90,length = 0.1)
    axis(2, at=seq(YLIM[1], YLIM[2], by=BY),cex.axis=1.5)
    axis(1, labels=c("","","",""), at = mp)
    if(DO.axis=="YES")axis(1, labels=c("Sandbar","Dusky","Gummy","Whiskery"), at = mp,cex.axis=2)
    box()
  }
  
  fn.bar.st=function(Mean,lowCI,upCI,CLO,DO.axis,BY)
  {
    
    what.se1=lowCI/mean(Mean)  #standardise to mean of 1
    what.se2=upCI/mean(Mean)
    what=Mean/mean(Mean)   
    
    YLIM=c(0,max(what.se2)*1.05)
    mp <- barplot(what, axes=FALSE, axisnames=FALSE, ylim=YLIM,col=CLO, main="", xlab="", ylab="")
    segments(c(mp), what.se1, c(mp) , what.se2,lwd=2)
    axis(2, at=seq(YLIM[1], YLIM[2], by=BY),cex.axis=1.5)
    axis(1, labels=c("","","",""), at = mp)
    if(DO.axis=="YES")axis(1, labels=c("Sandbar","Dusky","Gummy","Whiskery"), at = mp,cex.axis=1.35)
    box()
  }
  
  
  if(!isTRUE(Species.size))
  {
    tiff(file="Paper/Figure4.tiff",width = 2000, height = 2400,units = "px", res = 300,compression = "lzw")
    par(mfcol=c(2,1),mai=c(.35,.7,.06,.125),oma=c(.2,1.75,.01,.01),las=1,mgp=c(1.35,.8,0))
    NN=4
    #Species displacement
    fn.bar.st(Mean=exp(NewData.dist$pred),lowCI=exp(NewData.dist$LCB0.025),
              upCI=exp(NewData.dist$UCB0.975),CLO="grey60",DO.axis="NO",BY=.5)
    
    # fn.bar.inter.st(exp(Pred.dist.travel$fit),exp(Pred.dist.travel$fit-Pred.dist.travel$se.fit),
    #              exp(Pred.dist.travel$fit+Pred.dist.travel$se.fit),CLO,c(0,3),.5)
    mtext("Relative displacement",2,line=3,las=3,cex=1.75)
    # fn.bar.inter(exp(Pred.dist.travel$fit),exp(Pred.dist.travel$fit-Pred.dist.travel$se.fit),
    #        exp(Pred.dist.travel$fit+Pred.dist.travel$se.fit),CLO,c(0,350),50)
    # legend("topleft",c("Sandbar","Dusky","Gummy","Whiskery"),bty="n",horiz =T,fill=CLO,cex=1.8)
    #mtext("Displacement (km)",2,line=4.5,las=3,cex=1.6)
    
    #Species days at liberty
    # fn.bar.st(Pred.days$fit,Pred.days$fit-Pred.days$se.fit,
    #           Pred.days$fit+Pred.days$se.fit,"grey60",c(0,2),"NO",.5)
    # mtext("Time at liberty",2,line=4,las=3,cex=1.75)
    
    # fn.bar(Pred.days$fit,Pred.days$fit-Pred.days$se.fit,
    #        Pred.days$fit+Pred.days$se.fit,"grey60",c(0,1180),"NO",200)
    # mtext("Time at liberty (days)",2,line=4.5,las=3,cex=1.65)
    
    #Species ROM
    fn.bar.st(Mean=exp(NewData.rom$pred),lowCI=exp(NewData.rom$LCB0.025),
              upCI=exp(NewData.rom$UCB0.975),CLO="grey60",DO.axis="YES",BY=.5)
    
    # fn.bar.st(Pred.speed$fit,Pred.speed$fit-Pred.speed$se.fit,
    #           Pred.speed$fit+Pred.speed$se.fit,"grey60",c(0,2),"YES",.5)
    mtext("Relative rate of movement",2,line=3,las=3,cex=1.75)
    # fn.bar(exp(Pred.speed$fit)*3.6,exp(Pred.speed$fit-Pred.speed$se.fit)*3.6,
    #        exp(Pred.speed$fit+Pred.speed$se.fit)*3.6,"grey60",c(0,.031),"YES",.005)
    # mtext("Speed (km/h)",2,line=4.5,las=3,cex=1.65)
    dev.off()
    
  }
  
  if(Species.size)
  {
    plt.fn=function(d,fit,se)
    {
      d=cbind(d,fit,se)%>%
        mutate(low.95=fit-1.96*se,
               up.95=fit+1.96*se,
               fit.rel=fit/mean(fit),
               low.95.rel=low.95/mean(fit),
               up.95.rel=up.95/mean(fit),
               FL=exp(Log.Rel_FL))%>%
        arrange(Species)
      
      d%>%
        ggplot(aes(FL,fit.rel, colour = factor(Species)))+
        geom_point()+
        geom_errorbar(aes(ymin=low.95.rel, ymax=up.95.rel), width=2)
    }
    
    plt.fn(d=NewData.dist,fit=Pred.dist.travel$fit,se=Pred.dist.travel$se.fit)
    
    plt.fn(d=NewData.rom,fit=Pred.speed$fit,se=Pred.speed$se.fit)
    
    
    
  }
  
  #--Bearing 
  Do.Bearing=FALSE
  if(Do.Bearing)
  {
    #examine response var distribution
    rose.diag(circular(Tagging.mySpecies$bearing,units="degrees"), bins = 16, main = '')
    
    #statistical tests (see Zar 1999 page 616)
    #RAO's test
    #if Sig then distribution is not Uniform
    RAO=function(DATs) rao.spacing.test(circular(DATs$bearing,units="degrees")) 
    
    #Rayleigh's test
    #if Sig then distribution is not random
    RAYLEIGH=function(DATs) rayleigh.test(circular(DATs$bearing,units="degrees")) 
    
    Store.Circ.stats=vector('list',length=length(SPECIES))
    names(Store.Circ.stats)=SPECIES
    
    for(i in 1:length(SPECIES))
    {
      DATs=subset(Tagging.mySpecies,Species==SPECIES[i])
      dummy.store=vector('list',length(AREAS.bearing))
      names(dummy.store)=AREAS.bearing
      
      for (j in 1:length(AREAS.bearing))
      {
        dummy.store[[j]]=NULL
        if(!(i%in%c(2:4) & j==1))
        {
          if(i==4 & j%in%1:2) dummy.store[[j]]=NA
          if(!i==4 & j%in%1:2)
          {
            dat=subset(DATs,Areas%in%AREAS.bearing[[j]])
            store.rao=RAO(dat)
            store.rayleigh=RAYLEIGH(dat)
            dummy.store[[j]]=list(RaO=store.rao,RayleigH=store.rayleigh)
          }
          
          
        }
      }
      Store.Circ.stats[[i]]=dummy.store
    }
    
    #Look at Circular stats
    Store.Circ.stats
    
    
    #model        #NOTE: make sure same explanatory variables as for distance!!!ADD wEIGHT=Tagging$WEIGHT??
    
    #note: if the model falls ("Error in while (diff > tol) { : missing value where TRUE/FALSE needed"),
    #then add one covariate at a time (starting with DaysAtLarge) in stepwise manner. The lm.circular() is
    #very sensitive.
    #If failing, the use Colin's approach and do Raleigh and Rao's test for each variable.
    
    These.Covariates=Tagging.mySpecies[,match(c('bearing','Species','Sex','CONDITION','Rel_FL',
                                                'DaysAtLarge','Rel.name','Effort'),names(Tagging.mySpecies))]
    These.Covariates=These.Covariates[complete.cases(These.Covariates),]
    
    # model=lm.circular(y=circular(These.Covariates$bearing,units="degrees"), x=These.Covariates$DaysAtLarge,
    #                   init=2,type='c-l', verbose=T)
    
  }
  
}




# Not used ------------------------------------------------------------------
#1. COMPARE LATEST ACCESS DUMP WITH TAGGING

#Tab.All.releases.latest.Access=table(as.character(dat.capture$SPECIES))
# Tab.All.releases=table(as.character(Tagging$Species))
# if(!sum(Tab.All.releases.latest.Access)==sum(Tab.All.releases))print("CHECK OUT DATA")
# 
# 
# Tab.All.rec.latest.Access=table(as.character(dat.rec$SPECIES))
# Dummy=subset(Tagging,!is.na(Lat.rec))
# Tab.All.rec=table(as.character(Dummy$Species))
# if(!sum(Tab.All.rec.latest.Access)==sum(Tab.All.rec))print("CHECK OUT DATA")


#Plot data by size classes (50 cm bins)
# fun.plot.sp=function(whatDat,SPEC,MAT)
# {
#   ab.fn=function(){  abline(h=-33,lwd=2);segments(112,-26.5,114,-26.5,lwd=2)
#                      segments(116.5,-36,116.5,-33,lwd=2);segments(114,0,114,-22,lwd=2)}
#   if(whatDat=="Tagging")
#   {
#     dum=subset(Tagging,Species==SPEC)
#     HIST=hist(dum$Rel_FL,breaks=seq(0,300,10),plot=F)
#     
#     Step=50
#     dum$Bin=floor(dum$Rel_FL/Step)*Step
#     Unico=sort(unique(dum$Bin))
#     dum$Bin.col=factor(dum$Bin,labels=1:length(Unico))
#     legnd=paste(as.numeric(levels(factor(dum$Bin))),"-",as.numeric(levels(factor(dum$Bin)))+Step,sep="")
#     
#     dum$MAT=ifelse(dum$Rel_FL>=MAT,'red',"blue")
#     par(mfcol=c(1,1))
#     plot(HIST,main=paste(SPEC,"release"))
#     
#     par(mfcol=c(2,1))
#      plot(dum$Long.rels,dum$Lat.rels,main=paste(SPEC,"releases"),pch=19,col=dum$Bin.col,xlab="Long",ylab="Lat")
#     ab.fn()
#     plot(dum$Long.rec,dum$Lat.rec,main=paste(SPEC,"recaptures"),pch=19,col=dum$Bin.col,xlab="Long",ylab="Lat")
#     legend("topleft",legnd,bty='n',pch=19,col=levels(dum$Bin.col),cex=0.8)
#     ab.fn()
#   }
# 
#   
#   if(whatDat=="LatestAcces")
#   {
#     dum=subset(dat.capture,SPECIES==SPEC)
#     HIST=hist(dum$FL,breaks=seq(0,300,10),plot=F)
#     
#     Step=50
#     dum$Bin=floor(dum$FL/Step)*Step
#     Unico=sort(unique(dum$Bin))
#     dum$Bin.col=factor(dum$Bin,labels=1:length(Unico))
#     legnd=paste(as.numeric(levels(factor(dum$Bin))),"-",as.numeric(levels(factor(dum$Bin)))+Step,sep="")
#     
#     dum$MAT=ifelse(dum$FL>=MAT,'red',"blue")
#     par(mfcol=c(1,1))
#     plot(HIST,main=paste(SPEC,"release"))
#     
#     par(mfcol=c(2,1))
#     plot(dum$RELLNGDECDEG,-dum$RELLATDECDEG,main=paste(SPEC,"releases"),pch=19,col=dum$Bin.col,xlab="Long",ylab="Lat")
#     ab.fn()
#     plot(dum$RECLNGDECDEG,-dum$RECLATDECDEG,main=paste(SPEC,"recaptures"),pch=19,col=dum$Bin.col,xlab="Long",ylab="Lat")
#     legend("topleft",legnd,bty='n',pch=19,col=levels(dum$Bin.col),cex=0.8)
#     ab.fn()
#   }
# 
# }
# 
# 
# fun.plot.sp("Tagging","TK",size.mat[1])
# fun.plot.sp("LatestAcces","TK",size.mat[1])
# 
# fun.plot.sp("Tagging","BW",size.mat[2])
# fun.plot.sp("LatestAcces","BW",size.mat[2])
# 
# fun.plot.sp("Tagging","WH",size.mat[3])
# fun.plot.sp("LatestAcces","WH",size.mat[3])



#2. MANIPULATE DATA

#select relevant columns for receiver data
#Perth
# SMN_VR4G_metro=SMN_VR4G_metro[,match(c("Station No","Lat","Long"),names(SMN_VR4G_metro))]
# SMN_VR4G_metro$line="SMN.VR4"
# SMN_VR2W_metro$line="SMN.VR2"
# OTN_VR2W_metro$line="OTN"
# 
#   #Southwest
# TheseCols=c("Station No","Receiver No#","Date in","Actual Deg South","Actual Min South","Actual Deg East","Actual Min East")
# SMN_VR2W_Hamelin=SMN_VR2W_Hamelin[,match(TheseCols,names(SMN_VR2W_Hamelin))]
# SMN_VR2W_Albany=SMN_VR2W_Albany[,match(TheseCols,names(SMN_VR2W_Albany))]
# DoF_Chatham=DoF_Chatham[,match(TheseCols,names(DoF_Chatham))]
# SMN_VR2W_Hamelin$line="Hamelin"
# SMN_VR2W_Albany$line="Albany"
# DoF_Chatham$line="Chatham"
# 
# #combine receivers
# Perth=rbind(SMN_VR4G_metro,SMN_VR2W_metro,OTN_VR2W_metro)
# SouthWest=rbind(SMN_VR2W_Hamelin,SMN_VR2W_Albany,DoF_Chatham)
# SouthWest$Long=SouthWest$"Actual Deg East"+(SouthWest$"Actual Min East"/60)   #convert mins to decimals degrees
# SouthWest$Lat=-(SouthWest$"Actual Deg South"+(SouthWest$"Actual Min South"/60))


#Export data from movement rate estimationg             
#write.csv(Tagging,file="H:/Matias WA Fisheries/Analyses/Conventional tagging/Movement rate estimation/Tagging.csv",row.names =F)

# setwd("H:/Matias WA Fisheries/Analyses/Conventional tagging/General movement/outputs")
# 
# 
# 
# # Figure 1
#   #arrows or points
# # tiff(file="Figure1.points.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# par(mfcol=c(2,2),oma=c(1.5,1.5,0.001,0.001),mar=c(1,1,0.001,0.001),mgp=c(2, 0.75, 0))
# for(i in 1:length(SPECIES))
# {
#   Arrows(SPECIES[i],"Points",Species.names[i])
# # Arrows(SPECIES[i],"Arrows",Species.names[i])
# }
# mtext("Latitude (?S)",side=2,outer=T,line=-.1,font=1,las=0,cex=1.5)
# mtext("Longitude (?E)",side=1,outer=T,line=-.1,font=1,las=0,cex=1.5)
# # dev.off()
# 
# #ACA: problem with WAcoast in function plotmap!!!

# #Figure 1. Map of fishing zones
# 
# #define coordinates of plots
# South.WA.lat=c(-36,-15)
# #South.WA.lat=c(-36,-25)
# South.WA.long=c(112,138)
# #South.WA.long=c(112,129)
# S.WA.long=c(South.WA.long[2], South.WA.long[2], South.WA.long[1], South.WA.long[1])
# S.WA.lat=c(South.WA.lat[2], South.WA.lat[1], South.WA.lat[1], South.WA.lat[2])
# a=South.WA.long[1]:South.WA.long[2]
# b=seq(South.WA.lat[1],South.WA.lat[2],length.out=length(a))
# OZ.lat=c(-44.5,-11);OZ.long=c(113,155)
# PLATE=c(.01,.9,.075,.9)
# 
# tiff(file="Figure 1. Map.tiff",width = 2400, height = 2400,units = "px", res = 300, compression = "lzw")    
# 
# par(mar=c(2,2,2,2),oma=c(1,1,1,1))
# plotmap(a,b,PLATE,"dark grey",South.WA.long,South.WA.lat)
# 
# 
# text(113.5,-29.75,("West"),col="black", cex=1.5)
# text(113.5,-30.25,("coast"),col="black", cex=1.5)
# polygon(x=c(116.5,116.5,112,112),y=c(-26.5,-33,-33,-26.5),lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
# 
# text(114,-35.25,("Zone 1"),col="black", cex=1.5)
# polygon(x=c(116.5,116.5,112,112),y=c(-33,-37,-37,-33),lwd=1.5,col=rgb(.3,.3,.3,alpha=.5))
# 
# text(122,-35.25,("Zone 2"),col="black", cex=1.5)
# polygon(x=c(129,129,116.5,116.5),y=c(-30,-37,-37,-30),lwd=1.5,col=rgb(.7,.7,.7,alpha=.2))
# 
# LATT=South.WA.lat[2]:South.WA.lat[1]
# LONGG=South.WA.long[1]:South.WA.long[2]
# 
# if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=c(-37,-25),xlim=South.WA.long, zlim=c(-1,-300),
#         nlevels = 1,labcex=0.1,lty = c(1,2,3),col=c("gray20","gray20","gray20","transparent"),add=T)
# 
# 
# par(new=T,mar=c(2,2,2,2),oma=c(1,1,1,1))
# plotmap(a,b,PLATE,"dark grey",South.WA.long,c(-37,-25))
# 
# axis(side = 1, at =seq(112,129,2), labels = seq(112,129,2), tcl = .35,las=1,cex.axis=1,padj=-1.25)
# axis(side = 2, at = seq(-36,-25,2), labels = -seq(-36,-25,2),tcl = .35,las=2,cex.axis=1,hadj=.3)
# 
# text(116.73,-31.95,("Perth"),col="black", cex=1.1)
# points(115.86,-31.95,pch=19)
# text(116.73,-33.55,("Bunbury"),col="black", cex=1.1)
# points(115.6,-33.55,pch=19)
# text(117.7,-34.8,("Albany"),col="black", cex=1.1)
# points(117.8,-35,pch=19)
# text(122,-33.66,("Esperance"),col="black", cex=1.1)
# points(121.9,-33.86,pch=19)
# 
# mtext("Latitude (?S)",side=2,line=1.75,las=3,cex=1.3)
# mtext("Longitude (?E)",side=1,line=1.75,cex=1.3)
# 
# 
# 
# par(fig=c(.5,.92,.5,.92), new = T,mgp=c(.1,.4,0))
# plotMap(worldLLhigh, xlim=OZ.long,ylim=OZ.lat,plt = c(.1, 1, 0.075, 1),
#         col="dark grey",tck = 0.025, tckMinor = 0.0125, xlab="",ylab="",axes=F)
# box()
# polygon(x=S.WA.long,y=S.WA.lat,lwd=1.5,col=rgb(.1,.1,.1,alpha=.2))
# text(134,-23.5,("Australia"),col="black", cex=1.5)
# 
# 
# 
# dev.off()
# 
#   #density of points
# # for(i in 1:length(SPECIES))
# # {
#   #  pdf(file=paste(Species.names[i],".density.pdf",sep=""))    #create pdf
#   #tiff(file=paste(Species.names[i],".density.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#   #     compression = "lzw")
# #  density(SPECIES[i],Species.names[i],200)
# #  dev.off()
# }
# 
# 
# 
#   #3. EXPLORATORY ANALYSES AND SUMMARY STATISTICS
# function.explore=function(SPECIES,Species.names)
# {
#   datos=subset(Tagging,Species==SPECIES)
#   AREAS=unique(datos$Areas)
#   
#   #.. sex and size of released tagged sharks
#   tiff(file=paste(Species.names,".Sex&Size.released.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     HISTO=table(datos1$Sex,floor(datos1$Rel_FL/10)) #create histogram
#     barplot(HISTO, beside = TRUE,ylim=c(0,max(HISTO)),mgp = c(2, 0.6, 0),
#             names.arg= as.numeric(colnames(HISTO))*10,xlab="fork length (cm)",ylab="Frequency", 
#             xpd=F,axis.lty=1, axes=T,col=rev(gray(1:nrow(HISTO)/nrow(HISTO))),cex.names=1.1,las=1,cex=1.1)
#     box()
#     legend("topleft",AREAS[i],bty="n")
#     legend("topright",rownames(HISTO),fill=rev(gray(1:nrow(HISTO)/nrow(HISTO))),yjust=0, horiz=F,bty="n",cex=1.1)
#     boxplot(datos1$Rel_FL~datos1$Sex, col = "lightgray",add = F,ylab="fork length (cm)",xlab="sex")
#   }
#   dev.off()
#   
#   #.. time at large
#   tiff(file=paste(Species.names,".Time.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     boxplot(datos1$DaysAtLarge, xlab= "Time at liberty (days)")
#     legend("topleft",AREAS[i],bty="n")
#     hist(datos1$DaysAtLarge,breaks=50, xlab= "Time at liberty (days)", main="")
#     box()
#   }
#   dev.off()
#   
#   #.. distance moved
#   tiff(file=paste(Species.names,".Distance.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     boxplot(datos1$dist.trav_m/1000, xlab= "distance moved (km)")
#     legend("topleft",AREAS[i],bty="n")
#     hist(datos1$dist.trav_m/1000,breaks=50, xlab= "distance moved (km)", main="")
#     box()
#   }
#   dev.off()
#   
#   #.. distance vs time at large
#   tiff(file=paste(Species.names,".Distance VS Time.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(length(AREAS),1),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     plot(log(datos1$DaysAtLarge),log(datos1$dist.trav_m/1000), xlab= "ln(time at liberty (days))",ylab= "ln(distanced moved (km))")
#     legend("topleft",AREAS[i],bty="n")
#   }
#   dev.off()
#   
#   #     #.. speed
#   #   tiff(file=paste(Species.names,".Speed.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#   #      compression = "lzw")
#   #   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   #   for(i in 1:length(AREAS))
#   #   {
#   #     datos1=subset(datos,Areas==AREAS[i])
#   #     boxplot(datos1$speed, xlab= "speed (m/s)")
#   #     legend("topleft",AREAS[i],bty="n")
#   #     hist(datos1$speed,breaks=50, xlab= "speed (m/s)", main="")
#   #     box()
#   #   }
#   #   dev.off()
#   
#   #.. bearing
#   tiff(file=paste(Species.names,".Bearing.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(1,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     rose.diag(circular(datos1$bearing,units="degrees"), bins = 16, main = '')
#     legend("topleft",AREAS[i],bty="n")
#   }
#   dev.off()
#   
#   #.. condition effect on distance moved and time at large
#   #boxplots
#   tiff(file=paste(Species.names,".Cond_on_Dist&Time.boxplot.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     boxplot(datos1$DaysAtLarge~datos1$CONDITION, xlab="condition",ylab= "Time at liberty (days)")
#     legend("topleft",AREAS[i],bty="n")
#     boxplot(datos1$dist.trav_m/1000~datos1$CONDITION, xlab="condition",ylab= "distance moved (km)")
#     #  boxplot(datos1$speed~datos1$CONDITION, xlab="condition",ylab= "speed (m/s)")
#   }
#   dev.off()
#   
#   #histograms
#   tiff(file=paste(Species.names,".Cond_on_Dist&Time.hist.tiff",sep=""),width = 2400, height = 2400,units = "px", res = 300,
#        compression = "lzw")
#   par(mfcol=c(2,length(AREAS)),oma=c(0.8,.8,0.8,.8),mar=c(3,3,1,1),mgp=c(1.5, 0.5, 0))
#   for(i in 1:length(AREAS))
#   {
#     datos1=subset(datos,Areas==AREAS[i])
#     HISTO=table(datos1$CONDITION,floor(datos1$DaysAtLarge/10)) #create histogram
#     #COLORES=rev(gray(1:nrow(HISTO)/nrow(HISTO)))
#     COLORES=c(2,3,4,"violet")
#     barplot(HISTO, beside = TRUE,ylim=c(0,max(HISTO)),mgp = c(2, 0.6, 0),
#             names.arg= as.numeric(colnames(HISTO))*10,xlab=" Time at liberty (days)",ylab="Frequency", 
#             xpd=F,axis.lty=1, axes=T,col=COLORES,cex.names=1.1,las=1,cex=1.1)
#     box()
#     legend("topleft",AREAS[i],bty="n")
#     legend("topright",rownames(HISTO),fill=COLORES,yjust=0, horiz=F,bty="n",cex=1.1)
#     
#     HISTO=table(datos1$CONDITION,floor((datos1$dist.trav_m/1000)/100)) #create histogram
#     barplot(HISTO, beside = TRUE,ylim=c(0,max(HISTO)),mgp = c(2, 0.6, 0),
#             names.arg= as.numeric(colnames(HISTO))*100,xlab=" distance moved (km)",ylab="Frequency", 
#             xpd=F,axis.lty=1, axes=T,col=COLORES,cex.names=1.1,las=1,cex=1.1)
#     box()
#   }
#   dev.off()
# }
#   #apply exploratory function
# #  STATS=list()
# #  for(s in 1:length(SPECIES))
# #    {
# #      function.explore(SPECIES[s],Species.names[s])
# #      STATS[[Species.names[s]]]=function.stats(SPECIES[s],Species.names[s])
# #    }
# 
# 
# 
#   #4. GENERAL ANALYSES FOR MOVEMENT PATTERN PAPER
# General.numbers=function(SPECIES)
# {
#   datos=subset(All.rel.Tagging,Species==SPECIES)
#   datos1=subset(Tagging,Species==SPECIES)
#   
#   n.rel=sum(!is.na(datos$Lat.rels))
#   n.rec=sum(!is.na(datos$Lat.rec))
#   n.rec1=sum(!is.na(datos1$Lat.rec))
#   return(cbind(n.rel=n.rel,n.rec=n.rec,n.rec.used=n.rec1))
#   
# } 
# #4.1 All releases and recaptures
# ALL=NULL
# for(s in 1:length(SPECIES))
#   {
#     ALL=rbind(ALL,General.numbers(SPECIES[s]))
#   }
# ALL=as.data.frame(ALL)
# ALL$prop.used=round(100*ALL$n.rec.used/ALL$n.rec,0)
# ALL$Species=SPECIES
# write.csv(ALL,file="ALL.rel.rec.csv")
# 
# 
# #4.2 Table 1.
#5.1 Table 2. Summary of tags and recaptures

# columns: SPECIES, YEAR OF TAGGING, NUMBER TAGGED (males in brackets) by year of recaptured, 
#NUMBER RECAPTURED (males in brackets) by year of recaptured, final row with proportion recapture by year
# Table1=function(SPECIES)
# {
#   datos=subset(All.rel.Tagging,Species==SPECIES)
#   datos1=subset(Tagging,Species==SPECIES)
#   datos1=subset(datos1,!is.na(Yr.rec))
#   theseyears=sort(unique(datos$Yr.rel))
#   recyears=sort(unique(datos1$Yr.rec))
#   TABLE=list()
#   for (i in 1:length(theseyears))
#   {
#     datos2=subset(datos,Yr.rel==theseyears[i])
#     datos3=subset(datos1,Yr.rel==theseyears[i])
#     #     n.tagged=table(datos1$Sex,useNA='ifany')
#     #     n.females=n.tagged[match("F",names(n.tagged))]
#     #     n.males=n.tagged[match("M",names(n.tagged))]
#     #     n=n.females+n.males
#     n=sum(!is.na(datos2$Lat.rels))
#     n.yrs.rec=table(datos3$Yr.rec,useNA='ifany')
#     recaptures=as.matrix(t(n.yrs.rec))
#     yrs=colnames(recaptures)
#     if(length(yrs)>0) colnames(recaptures)=paste("rec.",yrs,sep="")
#     TABLE[[i]]=cbind("Year of tagging"=theseyears[i],"No. shark tagged"=n,recaptures)
#   }
#   return(TABLE)
# }
# for(s in 1:length(SPECIES))
#   {
#     write.csv(do.call(cbind,Table1(SPECIES[s])),file=paste("Table2.",Species.names[s],".csv",sep=""))
# 
#   }
# 
# 
# 
# #4.3 Figure 2.
# 
# 
# 
# 
# # tiff(file="Figure2.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# # par(mfcol=c(length(SPECIES),2),oma=c(3,3,1.5,1.5),mar=c(1,1,1,2),mgp=c(1, 0.75, 0),xpd=NA,las=2)
# # for(i in 1:length(SPECIES)) 
# #   {
# #     SizeFreq.rel.fn(SPECIES[i],maxLengFreq[i],maxLeng,BIN,Species.names[i],size.mat[i])
# #     axis(1,at=c(seq(0,maxLeng,BIN)),labels=F,tck=-0.02)
# #     axis(1,at=c(seq(0,maxLeng,BIN*5)),labels=F,tck=-0.04)
# #     if(i==1) text(maxLeng/2.3,maxLengFreq[i]+13,"Release",cex=1.5)
# #   }
# # axis(1,at=c(seq(0,maxLeng,BIN*5)),labels=c(seq(0,maxLeng,BIN*5)),cex.axis=1.5,las=1,tck=-0.04)
# # text(maxLeng/2,-13,"Fork length (cm)",cex=2)
# # 
# # 
# # for(i in 1:length(SPECIES)) 
# #   {
# #     SizeFreq.rec.fn(SPECIES[i],maxLengFreq[i],maxLeng,BIN,Species.names[i])
# #     axis(1,at=c(seq(0,maxLeng,BIN)),labels=F,tck=-0.02)
# #     axis(1,at=c(seq(0,maxLeng,BIN*5)),labels=F,tck=-0.04)
# #     if(i==1) text(maxLeng/2.3,maxLengFreq[i]+13,"Recapture",cex=1.5)
# #   }
# # axis(1,at=c(seq(0,maxLeng,BIN*5)),labels=c(seq(0,maxLeng,BIN*5)),cex.axis=1.5,las=1,tck=-0.04)
# # text(maxLeng/2,-13,"Fork length (cm)",cex=2)
# # mtext("Numbers",side=2,outer=T,line=1.5,font=1,las=0,cex=1.5)
# # dev.off()
# 
# 
# #test of effect of tagging on size distribution
# KOLMOGOV1=KOLMOGOV2=list()
# for(i in 1:length(SPECIES)) 
#   {
# #    KOLMOGOV1[[Species.names[i]]]=KolmoTest.rel.rec.fn(SPECIES[i])  
#     KOLMOGOV2[[Species.names[i]]]=KolmoTest.Nonrec.rec.fn(SPECIES[i])
# #    write.csv(do.call(cbind,KOLMOGOV1[[i]]),file=paste("KOLMOGOV.rel.rec.",Species.names[i],".csv",sep=""))
#     write.csv(do.call(cbind,KOLMOGOV2[[i]]),file=paste("KOLMOGOV.rec.nonrec.",Species.names[i],".csv",sep=""))
#   }
# 
# 
# 
# 
# 
# #Table of gears
# tableGear.rec=table(Tagging$Species,Tagging$CAPT_METHD,useNA='ifany')
# tableGear.Area.rec=table(Tagging$Areas,Tagging$CAPT_METHD,Tagging$Species,useNA='ifany')
# 
# 
# 
# 
# 
# 
# #Explore bearing at different lags
# 
# COL=c(2,3,4)
# 
# t=1;z=3;i=1
# plot(xdata$Long.rels,xdata$Lat.rels, ylim=c(-37,-15),xlim=c(112,136),col="transparent")
# for (i in 1:length(SPECIES[1:3]))
# {
#   datos=subset(Tagging,Species==SPECIES[i])
#   datos=subset(datos,Areas==AREAS[z])
#   xdata=subset(datos,DaysAtLarge>=LAG[t])
#   arrows(xdata$Long.rels,xdata$Lat.rels,xdata$Long.rec,xdata$Lat.rec,col=COL[i])
# }
# legend('topright',SPECIES[1:3],bty='n',lty=1,col=COL)
# legend('topleft',paste("LAG=",LAG[t],sep=""),bty='n')
# 
# 
# # 
# # par(mfcol=c(2,2),oma=c(.1,.1,.1,.1),mar=c(.1,.1,.1,.2),mgp=c(1, 0.75, 0),xpd=NA,las=2)
# # FNUM=c(5,5,4,4.1)
# # FINT=c(5,3.8,3.5,3.6)
# # for(i in 1:length(SPECIES))
# # {
# #   Rosa.vent.bearing(SPECIES[i],FNUM[i],FINT[i],Species.names[i])
# # }
# 
# 
# #x.text=rep(-1,length(SPECIES))
# #y.text=c(1.35,1.25,1.25,1.28)
# # tiff(file="Figure6.tiff",width = 2400, height = 2400,units = "px", res = 300,compression = "lzw")
# # par(mfcol=c(length(SPECIES),length(AREAS)),oma=c(3,3,1,.1),mar=c(1,1,.01,.1),mgp=c(1, 0.75, 0),xpd=NA,las=2)
# # for (j in 1:length(AREAS))
# # {
# #   for(i in 1:length(SPECIES)) 
# #     {
# #       General.dist.bearing(SPECIES[i],AREAS[j],16,Species.names[i])
# #       if(j==1) text(x.text[i],y.text[i],Species.names[i],cex=1.5)
# #     }
# # }
# # dev.off()
# 
# 
# 
# 
# 
# 
# # 6.LINEAR MODELLING
# 
# #6. 1 Displacement
# for (i in 1:length(SPECIES[1:3]))
# {
#   GLM=fn.GLM.distance(SPECIES[i])
#   write.table(GLM$ANOVA,file=paste("ANOVA.displacement.",Species.names[i],".csv",sep=""),
#               row.names=F,sep=",")
#   write.table(GLM$Coefs,file=paste("Coefs.displacement.",Species.names[i],".csv",sep=""),
#               row.names=T,sep=",")
# }
# 
# 
# 
# 
# 
# #NOT USED
# # plotlat=c(-35.5,-21.5)
# # plotlong=c(112.7,119.2)
# # 
# # Perth.city=c(115.866,-31.95)
# # Rotnest=c(115.50,-32.02)
# # 
# # long.Ningaloo=c(113,115.5);lat.Ningaloo=c(-23.2,-21.5)
# # long.Perth=c(114.5,116.75);lat.Perth=c(-32.75,-31.25)
# # long.SW=c(114.3,118.8);lat.SW=c(-35.385,-34.15)
# # 
# # WOz.long=c(113,122)
# # WOz.lat=c(-36,-19)
# # 
# # edgeX.Ningaloo=c(113.25,114.25,114.25,113.25);edgeY.Ningaloo=c(-21.5,-21.5,-23.5,-23.5)
# # edgeX.Perth=c(115,116,116,115);edgeY.Perth=c(-31,-31,-33,-33)
# # edgeX.SW=c(114.5,119,119,114.5);edgeY.SW=c(-33.7,-33.7,-35.5,-35.5)
# # 
# # 
# # #4.1.1 create mapping functions for showing all acoustic receivers
# Oz <- function()
# {
#   plotMap(worldLLhigh, xlim=c(110,155), ylim=c(-44.5,-11),col="white", axes=F, xlab="", ylab="",
#           border="black",bg="light grey")
#   text(133,-25,("Australia"),col="black", cex=1)
#   points(153.2199,-11.53294,col="light grey",pch=21,bg="light grey")
#   polygon(x=c(112,138,138,112),y=c(-37,-37,-15,-15),lwd=1)
# }
# 
# WestOz <- function()
# {
#   par(mar=c(0.02,0.01,0.02,0.03),mgp=c(2.5, 0.75, 0))
#   plotMap(worldLLhigh, xlim=WOz.long, ylim=WOz.lat,col="light grey", axes=F, xlab="Longitude (?E)", ylab="Latitude (?S)",
#           cex.lab=1.5)
#   text(118,-26,("Western"),col="black", cex=1.75)
#   text(118,-27,("Australia"),col="black", cex=1.75)   
#   polygon(x=edgeX.Ningaloo,y=edgeY.Ningaloo,lwd=2)
#   polygon(x=edgeX.Perth,y=edgeY.Perth,lwd=2)
#   polygon(x=edgeX.SW,y=edgeY.SW,lwd=2)
#   axis(side = 1, at = seq(WOz.long[1]+1,WOz.long[2],2), labels = seq(WOz.long[1]+1,WOz.long[2],2),
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 1, at = seq(WOz.long[1],WOz.long[2],1), labels = F,tck=-0.02) 
#   axis(side = 2, at = seq(WOz.lat[1],WOz.lat[2],2), labels = seq(WOz.lat[1],WOz.lat[2],2)*-1,
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 2, at = seq(WOz.lat[1],WOz.lat[2],1), labels = F,tck=-0.02)
#   box()
#   if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],xlim=WOz.long, ylim=WOz.lat, zlim=c(-1,-300),nlevels = 3,
#           labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
#   
# }
# 
# Ningaloo.array <- function()
# {
#   plotMap(worldLLhigh, xlim=long.Ningaloo, ylim=lat.Ningaloo,col="light grey", axes=F, xlab="", ylab="")
#   text(114.18,-21.99,("Exmouth"),col="black", cex=1.5)
#   points(114.125,-21.9,col="black", pch=19,cex=1.25)
#   points(Ningaloo$Recovery_longitude,Ningaloo$Recovery_latitude,cex=1,pch=19,col="gray42")    
#   box()
#   contour(xbat, ybat, reshaped[,2:ncol(reshaped)],xlim=long.Ningaloo, ylim=lat.Ningaloo, zlim=c(-1,-300),nlevels = 3,
#           labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
# }
# 
# Perth.array <- function()
# {
#   plotMap(worldLLhigh, xlim=long.Perth, ylim=lat.Perth,col="light grey", axes=F, xlab="", ylab="")
#   text(116.1,-32,("Perth"),col="black", cex=1.5)
#   points(115.83,-32,col="black", pch=19,cex=1.25)
#   points(Perth$Long,Perth$Lat,cex=1,pch=19,col="gray42")
#   polygon(x=Rottnest.Is$Longitude,y=Rottnest.Is$Latitude,col="dark grey")  #add missing islands
#   polygon(x=Garden.Is$Longitude,y=Garden.Is$Latitude,col="dark grey")
#   box()
#   contour(xbat, ybat, reshaped[,2:ncol(reshaped)],xlim=long.Perth, ylim=lat.Perth, zlim=c(-1,-300),nlevels = 3,
#           labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
# }
# 
# SW.array <- function()
# {
#   plotMap(worldLLhigh, xlim=long.SW, ylim=lat.SW,col="light grey", axes=F, xlab="", ylab="")
#   text(115.645,-34.32,("Augusta"),col="black", cex=1.5)
#   points(115.155,-34.32,col="black", pch=19,cex=1.25)
#   text(117.8,-34.89,("Albany"),col="black", cex=1.5)
#   points(117.85,-35.025,col="black", pch=19,cex=1.25)
#   points(SouthWest$Long,SouthWest$Lat,cex=1,pch=19,col="gray42")
#   box()
#   contour(xbat, ybat, reshaped[,2:ncol(reshaped)],xlim=long.SW, ylim=lat.SW, zlim=c(-1,-300),nlevels = 3,
#           labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
# }
# 
# 
# 
# #Figure 1. movement patterns paper. Map
# 
# #4.1.2 create mapping functions for showing tag recaptures
# overall.range.lat=floor(range(c(Tagging$Lat.rels,Tagging$Lat.rec),na.rm=T))
# overall.range.long=floor(range(c(Tagging$Long.rels,Tagging$Long.rec),na.rm=T))
# Arrows=function(SPECIES,PlotWhat,legend)
# {
#   datos=subset(Tagging,Species%in%SPECIES)
#   range.lat=floor(range(c(datos$Lat.rels,datos$Lat.rec),na.rm=T))
#   range.lat[2]=range.lat[2]+1
#   range.long=floor(range(c(datos$Long.rels,datos$Long.rec),na.rm=T))
#   if(range.lat[2]-range.lat[1]>2*(range.long[2]-range.long[1]))
#   {
#     range.long[2]=floor(range.long[2]*1.05)
#     range.long[1]=floor(range.long[1]*0.99)
#   }
#   range.long[1]=112
#   range.long[2]=range.long[2]+1
#   range.long=overall.range.long
#   range.lat=overall.range.lat
#   
#   plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="light grey", axes=F, xlab="",
#           ylab="",cex.lab=1.5)
#   #Points
#   if(PlotWhat=="Points")
#   {
#     unicasAreas2=unique(datos$Areas)
#     unicasAreas=unicasAreas2[match(c("JASDGDLF.zone2","JASDGDLF.zone1","WCDGDLF","closed",
#                                      "WANCSF"),unicasAreas2)]   #order areas
#     unicosColors=ColorAreas[match(unicasAreas2,unique(Tagging$Areas))]
#     unicasLetras=LetrasAreas[match(names(unicosColors),names(LetrasAreas))]
#     for(j in 1:length(unicasAreas))
#     {
#       datos2=subset(datos,Areas==unicasAreas[j])
#       with(datos2,points(Long.rels,Lat.rels,col="black",pch=24,bg=unicosColors[j]))
#       with(datos2,points(Long.rec,Lat.rec,col="black",pch=21,bg=unicosColors[j]))        
#     }
#     if(SPECIES=="TK")
#     {
#       legend("topright",c("release","recapture"),pch=c(24,21),cex=1.,bty="n")
#       #legend("center",unicasAreas,pch=c(15),cex=1,bty="n",col=unicosColors)
#       legend("center",unicasLetras,pch=c(15),cex=1,bty="n",col=unicosColors)
#       
#       #scale bar
#       segments(range.long[2]-5,-25.5,range.long[2]-5+4.45,-25.5,lwd=2)
#       text(range.long[2]-5+2,-26.35,"500 km",cex=1.25)
#     }
#   }
#   #Arrows
#   if(PlotWhat=="Arrows")
#   {
#     with(datos,arrows(Long.rels,Lat.rels,Long.rec,Lat.rec,length = 0.075,lty=1,col="gray20")) 
#   }
#   
#   axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = seq(range.long[1],range.long[2],2),
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 1, at = seq(range.long[1],range.long[2],1), labels = F,tck=-0.02) 
#   axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = seq(range.lat[1],range.lat[2],2)*-1,
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 2, at = seq(range.lat[1],range.lat[2],1), labels = F,tck=-0.02)
#   box()
#   #  legend("topleft",legend,bty="n",cex=1.)
#   text(117.1006,-15.99248,legend,cex=1.)
#   #add bathymetry
#   if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=range.lat,xlim=range.long, zlim=c(-1,-300),nlevels = 3,
#             labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
#     
#   segments(112.3,-26,113.2,-26,lwd=3,col="grey20",lty=1)
#   segments(114.58,-33,115.73,-33,lwd=3,col="grey20",lty=1)
#   segments(116.5,-35.39,116.5,-34.92,lwd=3,col="grey20",lty=1)
#   segments(129,-33.33,129,-31.74,lwd=3,col="grey20",lty=1)
#   segments(118.48,-18,122.33,-18,lwd=3,col="grey20",lty=1)
#   segments(123.75,-16.14,123.75,-14.97,lwd=3,col="grey20",lty=1)
#   
#   if(SPECIES=="GM")
#   {
#     text(122.73,-35.75,"A",cex=1.5)
#     text(114,-35.25,"B",cex=1.5)
#     text(113.35,-30,"C",cex=1.5)
#     text(113.35,-20,"D",cex=1.5)
#     text(121.5,-16,"E",cex=1.5)
#     text(124.9,-16,"F",cex=1.5)
#     vp <- baseViewports()
#     pushViewport(vp$inner,vp$figure,vp$plot)
#     pushViewport(viewport(x=0.4,y=0.7,width=.35,height=.35,just=c("left","top")))
#     par(fig=gridFIG(),new=T)  
#     Oz()
#   }
#   
# }
# 
# density=function(SPECIES,legend,N.GRIDS)
# {
#   datos=subset(Tagging,Species==SPECIES)
#   range.lat=floor(range(c(datos$Lat.rels,datos$Lat.rec),na.rm=T))
#   range.lat[2]=range.lat[2]+1
#   range.long=floor(range(c(datos$Long.rels,datos$Long.rec),na.rm=T))
#   if(range.lat[2]-range.lat[1]>2*(range.long[2]-range.long[1]))
#   {
#     range.long[2]=floor(range.long[2]*1.05)
#     range.long[1]=floor(range.long[1]*0.99)
#   }
#   range.long[1]=112
#   par(oma=c(0.1,0.1,0.1,0.1),mar=c(0.2,0.2,0.2,0.2),mgp=c(2.5, 1.1, 0))
#   plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="light grey", axes=F, xlab="Longitude (?E)",
#           ylab="Latitude (?S)",cex.lab=1.5)
#   
#   #add density
#   n.grids=N.GRIDS
#   Density <- kde2d(datos$Long.rec, datos$Lat.rec, h=c(1.1,1.1),n = n.grids, lims = c(range.long,range.lat))
#   Cols.Range=5:n.grids/n.grids
#   #col.image=rev(gray(Cols.Range)) #colors for image
#   col.image=rev(heat.colors(length(Cols.Range), alpha = 1))
#   ImageBreaks=seq(range(Density$z)[1],range(Density$z)[2],length.out=(length(col.image)+1))
#   colLeg=rep("transparent",length(ImageBreaks))  
#   colLeg[c(1,seq(length(colLeg)/3,length(colLeg),length.out=3))]="black"
#   image(Density,xlab="",ylab="",main = "",ylim=range.lat, xlim=range.long,zlim = c(0, max(Density$z)),
#         col =col.image,breaks=ImageBreaks,add=T)
#   par(new=T)    										#to keep all graphs
#   plotMap(worldLLhigh, xlim=range.long, ylim=range.lat,col="light grey", axes=F, xlab="Longitude (?E)",
#           ylab="Latitude (?S)",cex.lab=1.5)
#   axis(side = 1, at = seq(range.long[1],range.long[2],2), labels = seq(range.long[1],range.long[2],2),
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 1, at = seq(range.long[1],range.long[2],1), labels = F,tck=-0.02) 
#   axis(side = 2, at = seq(range.lat[1],range.lat[2],2), labels = seq(range.lat[1],range.lat[2],2)*-1,
#        tck=-0.035,las=1,cex.axis=1.2)
#   axis(side = 2, at = seq(range.lat[1],range.lat[2],1), labels = F,tck=-0.02)
#   box()
#   legend("topleft",legend,bty="n")
#   #add bathymetry
#   if(add.depth=="YES") contour(xbat, ybat, reshaped[,2:ncol(reshaped)],ylim=range.lat,xlim=range.long, zlim=c(-1,-300),nlevels = 3,
#           labcex=0.5,lty = c(1,2,3),col=c("gray10","gray30","gray60","transparent"),add=T)
#   #add releases and recaptures
#   # points(datos$Long.rels, datos$Lat.rels)
#   # points(datos$Long.rec,datos$Lat.rec,col=2)
#   color.legend(range.long[2]-1,range.lat[2],range.long[2],range.lat[2]-10,legend=round(ImageBreaks,3),
#                rect.col=col.image,gradient="y",bg="white",cex=.75,col=colLeg)
# }
# 
# #Kolmo-Smirnov test of distributions recaptures and releases
# KolmoTest.rel.rec.fn=function(SPEC)
# {
#   datos=subset(Tagging,Species==SPEC)
#   datos$CAP_FL=ifelse(datos$CAP_FL<datos$Rel_FL,NA,datos$CAP_FL)
#   
# }
# 
# 
# SizeFreq.rel.fn=function(SPEC,maxDisFreq,maxLeng,BIN,Species.names,size.mat)
# {
#   datos=subset(Tagging,Species==SPEC)
#   hist(datos$Rel_FL,breaks=seq(0,maxLeng,BIN),col="gray",ylab="",main="",xlab="",ylim=c(0,maxDisFreq),
#        xlim=c(0,maxLeng),xaxt="n",cex.axis=1.5)
#   box()
#   legend("topright",paste(Species.names," (n=",nrow(datos),")",sep=""),bty="n",cex=1.5)
#   arrows(size.mat,maxDisFreq*.2,size.mat,1,lwd=1,length=0.1)
# }
# 
# SizeFreq.rec.fn=function(SPEC,maxDisFreq,maxLeng,BIN,Species.names)
# {
#   datos=subset(Tagging,Species==SPEC)
#   datos$CAP_FL=ifelse(datos$CAP_FL<datos$Rel_FL,NA,datos$CAP_FL)
#   DATA=as.numeric(as.character(datos$CAP_FL))
#   DATA=subset(DATA,DATA<999)
#   hist(DATA,breaks=seq(0,maxLeng,BIN),col="gray",ylab="",main="",xlab="",
#        ylim=c(0,maxDisFreq),xlim=c(0,maxLeng),xaxt="n",cex.axis=1.5)
#   box()
# }
# 
# 
# 
# 
# 
# #Kolmo-Smirnov test of distributions recaptures and non recaptures
# KolmoTest.Nonrec.rec.fn=function(SPEC)
# {
#   datos.All=subset(All.rel.Tagging,Species==SPEC)
#   datos.All$CAP_FL=ifelse(datos.All$CAP_FL<datos.All$Rel_FL,NA,datos.All$CAP_FL)
#   
#   datos1=subset(datos.All,is.na(Yr.rec))
#   datos1=subset(datos1,is.na(Long.rec))
#   datos1=subset(datos1,!is.na(Rel_FL))
#   datos1=subset(datos1,!is.na(Rel_FL))
#   
#   datos2=subset(datos.All,!(is.na(Yr.rec)))
#   datos2=subset(datos2,!is.na(Rel_FL))
#   
#   Kolmo=ks.test(datos1$Rel_FL,datos2$Rel_FL)
#   return(Kolmo)  
# }
# 
# #KolmoTest.Nonrec.rec.fn=function(SPECIES)
# # {
# #    datos.All=subset(All.rel.Tagging,Species==SPECIES)
# #    datos=subset(Tagging,Species==SPECIES)
# #    datos.All=datos.All[,match(c("Tag no","Rel_FL"),names(datos.All))]
# #    datos=datos[,match(c("Tag no","Rel_FL"),names(datos))]
# #    ThisMatch=match(datos$"Tag no",datos.All$"Tag no")
# #    data.rec=datos.All[ThisMatch,]
# #    data.Not.rec=datos.All[-ThisMatch,]
# #    Kolmo=ks.test(data.rec$Rel_FL,data.Not.rec$Rel_FL)
# #    return(Kolmo)  
# # }



