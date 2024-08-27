# Mark recapture growth estimation using Alex's WAFishBiology package
#note: assuming vonB growth
#      due to gear selectivity, adding random 'recaptures': neonates for gummy and whiskery
#                                                           large individuals for dusky and sandbar

#MISSING: consider setting gummy and whiskery Lo to 0 age rather than age 1.

# library(devtools)
# devtools::install_github("SAlexHesp/L3AssessRPackage", build_vignettes=TRUE, force=TRUE)
# devtools::install_github("SAlexHesp/WAFishBiologyRPackage", build_vignettes=TRUE, force=TRUE)

library(tidyverse)
library(WAFishBiology)

#
#---Bring in tagging data-------
if(!exists('handl_OneDrive')) source('C:/Users/myb/OneDrive - Department of Primary Industries and Regional Development/Matias/Analyses/SOURCE_SCRIPTS/Git_other/handl_OneDrive.R')
source(handl_OneDrive("Analyses/SOURCE_SCRIPTS/Git_other/Source_conventional_data.R"))
LH.data=read.csv(handl_OneDrive('Data/Life history parameters/Life_History.csv'))
min.days.liberty=30 #30

SPECIES=c("TK","BW","WH","GM")
names(SPECIES)=c(18007,18003,17003,17001)
#1. Create useful objects
LH.data=LH.data%>%filter(SPECIES %in%c(18007,18003,17003,17001))

#keep only individuals tagged with conventional tags and released in conditions 1 or 2
dat=Tagging%>%
  mutate(DATE_REL=as.Date( paste( Yr.rel,Mn.rel , Day.rel , sep = "-" ), format = "%Y-%m-%d" ),
         DATE_CAPTR=as.Date( paste( Yr.rec,Mn.rec , Day.rec , sep = "-" ), format = "%Y-%m-%d" ),
         DaysAtLarge=as.numeric(round((DATE_CAPTR-DATE_REL),0)),
         DaysAtLarge=ifelse(DaysAtLarge<0,NA,DaysAtLarge))%>%
  filter(Yr.rel>=1993 & CONDITION%in%c(1,2))%>%
  filter(Tag.type=='conventional' & Species%in%SPECIES)%>%
  filter(Recaptured=='Yes')%>%
  dplyr::select(Species,Sex,Rel_FL,CAP_FL,DATE_REL, DATE_CAPTR, DaysAtLarge)%>%
  filter(DaysAtLarge>min.days.liberty)%>%
  mutate(deltaFL=CAP_FL-Rel_FL)%>%
  filter(deltaFL>0)

dat%>%
  ggplot(aes(Rel_FL,CAP_FL,color=DaysAtLarge))+
  geom_point()+
  facet_wrap(~Species,scales='free')

# ---------Estimate growth----------------------------------------------------------
#GrowthCrvChoice  1 = double logistic
#                 2 = Gaussian function  params=log(c(0.1, 80, 40, 10))
#                 3 = von Bertalanffy growth curve Linf, K, SD
#                 4 = Gompertz growth curve

GrowthCrvChoice=3 

Fit.to.random=TRUE  #add dummy data to improve biological realism. either size at birth for GM & WH or large individuals for BW & TK
n.dummy=20    #number of random individuals added

Model.fit=vector('list',length(SPECIES))
names(Model.fit)=SPECIES
for(s in 1:length(SPECIES))
{
  dd=dat%>%filter(Species==SPECIES[s])
  hh=as.numeric(names(SPECIES)[s])
  lh=LH.data%>%filter(SPECIES==hh)
  
  if(Fit.to.random)
  {
    if(SPECIES[s]%in%c('TK','BW')) 
    {
      if(SPECIES[s]=='TK')
      {
        AgE=c(25,20,21)
        age.length.old=data.frame(Age=AgE,
                                  FL=c(150,150,148))%>%
          mutate(Rel_FL=rnorm(length(AgE),lh$LF_o,1),
                 CAP_FL=FL,
                 DaysAtLarge=365*Age)
      }
      if(SPECIES[s]=='BW')
      {
        AgE=c(30,32,25)
        age.length.old=data.frame(Age=AgE,
                                  FL=c(260,275,250))%>%
          mutate(Rel_FL=rnorm(length(AgE),lh$LF_o,1),
                 CAP_FL=FL,
                 DaysAtLarge=365*Age)
      }
      dummy=dd[1:nrow(age.length.old),]%>%
        mutate(Rel_FL=age.length.old$Rel_FL,
               CAP_FL=age.length.old$CAP_FL,
               DaysAtLarge=age.length.old$DaysAtLarge)
      
    }
     
    # if(SPECIES[s]%in%c('TK','BW'))     #aca: don't do random sample, replace by known size at age for large individuals
    # {
    #   dummy=dd[1:n.dummy,]%>%
    #     mutate(Rel_FL=runif(n.dummy,lh$FL_inf*.7,lh$FL_inf*.8),
    #            CAP_FL=Rel_FL+runif(n.dummy,5,10),
    #            DaysAtLarge=round(runif(n.dummy,365,2*365)))
    # }
    if(SPECIES[s]%in%c('GM','WH'))    #review issue of age at birth =0    ACA
    {
      dummy=dd[1:n.dummy,]%>%
        mutate(Rel_FL=rep(0,n.dummy),
               CAP_FL=rnorm(n.dummy,lh$LF_o,1),
               DaysAtLarge=rep(365,n.dummy))
    }

    dd=rbind(dd,dummy)
    
  }
  

  #MaxAge=lh$Max_Age
  #MaxLen = with(lh,(Max.TL-b_FL.to.TL)/a_FL.to.TL)
  params=log(c(lh$FL_inf, lh$K,10))
  nstep=50
  Obs_delta_t=dd$DaysAtLarge
  Obs_Initlen=dd$Rel_FL
  Obs_Finlen=dd$CAP_FL
  nobs=length(dd$DaysAtLarge)
  Model.fit[[s]]=GetTaggingGrowthModelResults(params,nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs)
  
}

Model.fit[[s]]$ParamEst


#check fit
PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, Model.fit[[s]])
Model.fit$TK$ParamEst
Model.fit$TK$convergence

Model.fit$BW$ParamEst
Model.fit$WH$ParamEst
Model.fit$GM$ParamEst





# From Alex ---------------------------------------------------------------
do.from.Alex=FALSE
if(do.from.Alex)
{
  # simulate tag-recapture data, based on von Bertalanffy growth curve
  set.seed(123)
  nstep = 50 # number of steps for numerical integration
  MaxLen = 240
  GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
  vbLinf = 250
  vbK = 0.3
  StandDev = 10
  params = log(c(vbLinf, vbK, StandDev))
  nobs = 200
  CalculationStage = 1
  res=SimulateTagRecaptureData(GrowthCrvChoice, nstep, nobs, MaxLen, params)
  Obs_delta_t=res$Obs_delta_t
  Obs_Initlen=res$Obs_Initlen
  Obs_Finlen=res$Obs_Finlen
  plot(Obs_Initlen,Obs_Finlen-Obs_Initlen)
  # range(Obs_Initlen)
  # range(Obs_Finlen)
  
  # fit model
  vbLinf = 280
  vbK = 0.2
  StandDev = 10
  params = log(c(vbLinf, vbK, StandDev))
  MaxAge = 30
  FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs)
  FittedRes$convergence
  FittedRes$ParamEst
  
  PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, FittedRes)
  
  # specified 'true' values used to simulate tagging data
  # vbLinf = 250
  # vbK = 0.3
  # StandDev = 10
  
  # estimated values
  # > FittedRes$ParamEst
  # Estimate lw_95%CL up_95%CL
  # vb_Linf    249.03   240.79   257.56
  # vb_K         0.31     0.28     0.33
  # StandDev     9.59     8.70    10.58
  
  
  
  # *********************************
  # fit to real data - sandbar sharks
  # *********************************
  
  dat = read.csv("Sandbar.csv",header=T)
  head(dat)
  
  Obs_delta_t=dat$Obs_delta_t
  Obs_Initlen=dat$Obs_Initlen
  Obs_Finlen=dat$Obs_Finlen
  plot(Obs_Initlen,Obs_Finlen)
  nobs = length(Obs_Initlen)
  
  # fit model
  vbLinf = 200
  vbK = 0.2
  StandDev = 10
  params = log(c(vbLinf, vbK, StandDev))
  MaxAge = 30
  FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs)
  FittedRes$convergence
  FittedRes$ParamEst
  PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, FittedRes)
  
  
  # **************************************************************************************
  # explore growth results - simulate age and length data with estimated growth parameters
  # with specified CV
  # **************************************************************************************
  
  library(L3Assess)
  set.seed(123)
  SampleSize=1000
  MaxAge = 30
  TimeStep = 1 # model timestep (e.g. 1 = annual, 1/12 = monthly)
  NatMort = 0.1
  FishMort = 0.1
  MaxLen = 300
  LenInc = 2
  midpt = seq(0,MaxLen - LenInc, LenInc) + (LenInc/2)
  MLL=NA # (minimum legal length) # retention set to 1 for all lengths if MLL set to NA, otherwise retention is knife-edged at MLL
  SelectivityType=2 # 1=selectivity inputted as vector, 2=asymptotic logistic selectivity curve
  L50 = 30 # selectivity
  L95 = 50 # selectivity
  SelectivityVec = NA # selectivity vector
  DiscMort = 0
  GrowthCurveType = 1 # 1 = von Bert, 2 = Schnute
  Linf = 144.96
  vbK = 0.17
  CVSizeAtAge = 0.08
  GrowthParams = c(Linf, vbK)
  RefnceAges = NA
  Res=SimLenAndAgeFreqData(SampleSize, MaxAge, TimeStep, NatMort, FishMort, MaxLen, LenInc, MLL, SelectivityType,
                           L50, L95, SelectivityVec, DiscMort, GrowthCurveType, GrowthParams, RefnceAges, CVSizeAtAge)
  
  plot(Res$ObsDecAge, Res$ObsLenClMidPt, ylim=c(0,200))
  abline(h=144.96)
  
  
  
  # *******************************
  # fit to real data - gummy sharks
  # *******************************
  
  dat = read.csv("GM.csv",header=T)
  head(dat)
  
  Obs_delta_t=dat$Obs_delta_t
  Obs_Initlen=dat$Obs_Initlen
  Obs_Finlen=dat$Obs_Finlen
  
  
  FitToSizeAtBirth = TRUE
  if (FitToSizeAtBirth) { # perhaps, more elegant approach would be use distribution for size at birth in likelihood function
    nRandSamples = 20
    rand_FinalLen = rnorm(nRandSamples,27,1)
    hist(rand_FinalLen)
    Obs_delta_t = c(Obs_delta_t,rep(365,nRandSamples))
    Obs_Initlen = c(Obs_Initlen,rep(0,nRandSamples))
    Obs_Finlen = c(Obs_Finlen,rand_FinalLen)  
  }
  
  plot(Obs_Initlen,Obs_Finlen, ylim=c(0,150))
  nobs = length(Obs_Initlen)
  
  # fit model
  set.seed(123)
  nstep = 50 # number of steps for numerical integration
  MaxLen = 240
  GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
  vbLinf = 200
  vbK = 0.2
  StandDev = 10
  params = log(c(vbLinf, vbK, StandDev))
  MaxAge = 30
  FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs)
  FittedRes$convergence
  FittedRes$ParamEst
  PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, FittedRes)
  
  # not specifying size at birth
  # > FittedRes$ParamEst
  # Estimate lw_95%CL up_95%CL
  # vb_Linf    124.41   119.75   129.26
  # vb_K         1.05     0.71     1.56
  # StandDev     6.81     5.59     8.30
  
  # specifying size at birth 
  # > FittedRes$ParamEst
  # Estimate lw_95%CL up_95%CL
  # vb_Linf    147.69   137.42   158.74
  # vb_K         0.21     0.18     0.23
  # StandDev     6.46     5.62     7.43
  
  # **********************************
  # fit to real data - whiskery sharks
  # **********************************
  
  dat = read.csv("WH.csv",header=T)
  head(dat)
  
  Obs_delta_t=dat$Obs_delta_t
  Obs_Initlen=dat$Obs_Initlen
  Obs_Finlen=dat$Obs_Finlen
  plot(Obs_Initlen,Obs_Finlen)
  nobs = length(Obs_Initlen)
  
  # fit model
  set.seed(123)
  nstep = 50 # number of steps for numerical integration
  MaxLen = 240
  GrowthCrvChoice = 3 # 1 = double logistic, 2 = Gaussian function, 3 = von Bertalanffy growth curve, 4 = Gompertz growth curve
  
  vbLinf = 150
  vbK = 0.1
  StandDev = 10
  params = log(c(vbLinf, vbK, StandDev))
  MaxAge = 30
  FittedRes=GetTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs)
  FittedRes$convergence
  FittedRes$ParamEst
  PlotFittedTaggingGrowthModelResults(params, nstep, Obs_delta_t, Obs_Initlen, Obs_Finlen, MaxAge, nobs, FittedRes)
  
}

