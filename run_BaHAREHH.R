library("rstan")
rstan_options(javascript = FALSE)
options(mc.cores = parallel::detectCores())
setwd("/proj/wildman/BaHAREHH/")

latname="Vulpes vulpes"
startyear=2003
endyear=2021
usenew=TRUE
datsavestr<-paste0('/proj/wildman/BaHAREHH/StanInputs',latname,'New',usenew,startyear,'to',endyear,'.rdata')
load(datsavestr)

############################
#Run stan
Areafit <- stan(file = "BaHAREHH.stan" , data = Data,chains=nchains, iter = iterations, warmup = warmup,thin=thin,init = Inits,control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),include = TRUE, pars = c("a","phi","logQ","loga","logmu","beta"))



####################
#Posterior prediction
set.seed(1)
draws=(iterations-warmup)*nchains
NrOfYears=Data$NrOfYears
NrOfLans=Data$NrOfLans
NrOfNon0LanLans=Data$NrOfNon0LanLans
NrOfNon0LanAndNotFullycoveredKretsar=Data$NrOfNon0LanAndNotFullycoveredKretsar
KretsListLanSeqID=Data$KretsListLanSeqID
KretsYearSeq=Data$KretsYearSeq
NonZeroLanLansIndex=Data$NonZeroLanLansIndex
WhichisNon0LanKrets=Data$WhichisNon0LanKrets
NrOfKretsar=Data$NrOfKretsar
M=Data$M
FullyCovered=Data$FullyCovered
KretsListA=Data$KretsListA
WhichisNon0LanLan=Data$WhichisNon0LanLan
Ahat=Data$Ahat
WhichisNon0LanAndNotFullycoveredKrets=Data$WhichisNon0LanAndNotFullycoveredKrets
NonFullyCoveredIndex=Data$NonFullyCoveredIndex


sumN=array(dim=NrOfYears)
lgammaaphi =array(dim=NrOfNon0LanLans , NrOfYears)


Khat=array( NrOfNon0LanAndNotFullycoveredKretsar )
Khattot=array( NrOfNon0LanAndNotFullycoveredKretsar )
Khattotyear=array( NrOfYears )
nuhat=array( NrOfNon0LanAndNotFullycoveredKretsar )


shape=array( NrOfNon0LanAndNotFullycoveredKretsar )

ext<-rstan::extract(Areafit, pars="a")

a=array(ext$a[,,],dim=c(draws,NrOfLans,NrOfYears))
aMat=mapply(function(i, j) a[,i, j], KretsListLanSeqID, KretsYearSeq)

ext<-rstan::extract(Areafit, pars="phi")
phi=array(ext$phi[,,],dim=c(draws,NrOfNon0LanLans,NrOfYears))
phiMat=mapply(function(i, j) phi[,i, j], NonZeroLanLansIndex[KretsListLanSeqID[WhichisNon0LanKrets]], KretsYearSeq[WhichisNon0LanKrets])
phiMatenlarged=array(dim=c(draws,NrOfKretsar))
phiMatenlarged[,WhichisNon0LanKrets]=phiMat

ext<-rstan::extract(Areafit, pars="logQ")
logQMat=matrix(unname(unlist(ext)),nrow=draws)


Qmat=exp(logQMat)
Mmat=matrix(rep(M,draws),nrow=draws,byrow=TRUE)
Nmat=Mmat
Nmat[,!FullyCovered]=Nmat[,!FullyCovered]+Qmat
KretsListAMat=matrix(rep(KretsListA,draws),nrow=draws,byrow=TRUE)
mMat=KretsListAMat/Nmat
logNmat=log(Nmat)
Qmatenlarged=array(dim=c(draws,NrOfKretsar))
Qmatenlarged[,]=0
Qmatenlarged[,!FullyCovered]=Qmat

lgammaaphi=lgamma(array(a[, WhichisNon0LanLan, ] ,dim=c(draws,NrOfNon0LanLans,NrOfYears))+phi)
lgammaaphiMat=mapply(function(i, j) lgammaaphi[,i, j], NonZeroLanLansIndex[KretsListLanSeqID[WhichisNon0LanKrets]], KretsYearSeq[WhichisNon0LanKrets])
lgammaaphiMatenlarged=array(dim=c(draws,NrOfKretsar))
lgammaaphiMatenlarged[,WhichisNon0LanKrets]=lgammaaphiMat


lgammaa=lgamma(a)
lgammaaMat=mapply(function(i, j) lgammaa[,i, j], NonZeroLanLansIndex[KretsListLanSeqID[WhichisNon0LanKrets]], KretsYearSeq[WhichisNon0LanKrets])
lgammaaMatenlarged=array(dim=c(draws,NrOfKretsar))
lgammaaMatenlarged[,WhichisNon0LanKrets]=lgammaaMat


ext<-rstan::extract(Areafit, pars="loga")

loga=array(ext$loga[,,],dim=c(draws,NrOfLans,NrOfYears))
logaMat=mapply(function(i, j) loga[,i, j], KretsListLanSeqID, KretsYearSeq)


ext<-rstan::extract(Areafit, pars="logmu")
logmuMat=matrix(unname(unlist(ext)),nrow=draws)
logmuMatenlarged=array(dim=c(draws,NrOfKretsar))
logmuMatenlarged[,WhichisNon0LanKrets]=logmuMat


ext<-rstan::extract(Areafit, pars="beta")
betaMat=array(ext$beta[,,],dim=c(draws,NrOfNon0LanLans,NrOfYears))
betaMat=mapply(function(i, j) betaMat[,i, j], NonZeroLanLansIndex[KretsListLanSeqID[WhichisNon0LanKrets]], KretsYearSeq[WhichisNon0LanKrets])
betaMatenlarged=array(dim=c(draws,NrOfKretsar))
betaMatenlarged[,WhichisNon0LanKrets]=betaMat


logAhat=log(Ahat)
logbMat=logaMat[,WhichisNon0LanAndNotFullycoveredKrets]+
  logQMat[,NonFullyCoveredIndex[WhichisNon0LanAndNotFullycoveredKrets]]-
  matrix(rep(logAhat[WhichisNon0LanAndNotFullycoveredKrets],draws),nrow=draws,byrow=TRUE)
bMat=exp(logbMat)

lognuhat=aMat[,WhichisNon0LanAndNotFullycoveredKrets] * logbMat - 
  aMat[,WhichisNon0LanAndNotFullycoveredKrets] * (1/phiMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets] * (logmuMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]+log(mMat[,WhichisNon0LanAndNotFullycoveredKrets])) + logNmat[,WhichisNon0LanAndNotFullycoveredKrets]) + 
  (aMat[,WhichisNon0LanAndNotFullycoveredKrets] + 
     phiMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]) * 
  (-logbMat+1/phiMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]*
     (logmuMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]+log(mMat[,WhichisNon0LanAndNotFullycoveredKrets])) + 
     logNmat[,WhichisNon0LanAndNotFullycoveredKrets]) +
  lgammaaphiMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets] - lgammaaMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets] 



nuhat=exp(lognuhat)
shape=nuhat*Qmatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]*betaMatenlarged[,WhichisNon0LanAndNotFullycoveredKrets]

Khat=matrix(rnbinom(n=NrOfNon0LanAndNotFullycoveredKretsar*draws,mu=nuhat*Qmatenlarged[,WhichisNon0LanAndNotFullycoveredKrets],size=shape),nrow = draws);
Khat[nuhat==0]=0
Khat[shape==0]=0