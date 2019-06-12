i=3
datos=subset(Tagging,Species%in%SPECIES[i])
latran=c(-31,-33)
lonran=c(116,118)
nonsense=subset(datos,Lat.rec< latran[1] & Lat.rec > latran[2] & Long.rec > lonran[1] & Long.rec < lonran[2])

Tag.nonsense.rec=c(476,1151,3021)
these.nonsense.rec=match(Tag.nonsense.rec,Tagging$FINTAGNO)
na.these.cols=match(c("Lat.rec","Long.rec"),names(Tagging))
Tagging[these.nonsense.rec,na.these.cols]=NA  #set to NA the recapture location  