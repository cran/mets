## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  #dev="png",
  comment = "#>"  
)
fullVignette <- Sys.getenv("_R_FULL_VIGNETTE_") %in% c("1","TRUE")
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 runb <- 0
 options(warn=-1)
 set.seed(1000) # to control output in simulatins for p-values below.

n <- 200; k.boot <- 10; 

dat <- kumarsimRCT(n,rho1=0.5,rho2=0.5,rct=2,censpar=c(0,0,0,0),
          beta = c(-0.67, 0.59, 0.55, 0.25, 0.98, 0.18, 0.45, 0.31),
    treatmodel = c(-0.18, 0.56, 0.56, 0.54),restrict=1)
dfactor(dat) <- dnr.f~dnr
dfactor(dat) <- gp.f~gp
drename(dat) <- ttt24~"ttt24*"
dat$id <- 1:n
dat$ftime <- 1

## -----------------------------------------------------------------------------
weightmodel <- fit <- glm(gp.f~dnr.f+preauto+ttt24,data=dat,family=binomial)
wdata <- medweight(fit,data=dat)

## -----------------------------------------------------------------------------
aaMss2 <- binreg(Event(time,status)~gp+dnr+preauto+ttt24+cluster(id),data=dat,time=50,cause=2)
summary(aaMss2)

## ----label=firstmodel---------------------------------------------------------
### binomial regression ###########################################################
aaMss <- binreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
		time=50,weights=wdata$weights,cause=2)
summary(aaMss)

ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

## ----label=multiplemodels-----------------------------------------------------
### lin-ying model ################################################################
aaMss <- aalenMets(Surv(time/100,status==2)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
		   weights=wdata$weights)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### cox model ###############################################################################
aaMss <- phreg(Surv(time,status==2)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
	       weights=wdata$weights)
summary(aaMss)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### Fine-Gray #############################################################3
aaMss <- cifreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
		weights=wdata$weights,propodds=NULL,cause=2)
summary(aaMss)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### logit model  #############################################################3
aaMss <- cifreg(Event(time,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
		weights=wdata$weights,cause=2)
summary(aaMss)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

### binomial outcome  ############################
aaMss <- binreg(Event(ftime,status)~dnr.f0+dnr.f1+preauto+ttt24+cluster(id),data=wdata,
		time=50,weights=wdata$weights,cens.weights=1,cause=2)
summary(aaMss)
ll <- mediatorSurv(aaMss,fit,data=dat,wdata=wdata)
summary(ll)
if (runb>0) { bll <- BootmediatorSurv(aaMss,fit,data=dat,k.boot=k.boot); summary(bll)}

## ----label=multinom, cache=TRUE, eval=fullVignette----------------------------
# library(mets)
# data(tTRACE)
# dcut(tTRACE) <- ~.
# 
# weightmodel <- fit <- mlogit(wmicat.4 ~agecat.4+vf+chf,data=tTRACE,family=binomial)
# wdata <- medweight(fit,data=tTRACE)
# 
# aaMss <- binreg(Event(time,status)~agecat.40+ agecat.41+ vf+chf+cluster(id),data=wdata,
# 		time=7,weights=wdata$weights,cause=9)
# summary(aaMss)
# MultMed <- mediatorSurv(aaMss,fit,data=tTRACE,wdata=wdata)
# summary(MultMed)

## ----results="hide", echo=FALSE-----------------------------------------------
## To save time building the vignettes on CRAN, we cache time consuming computations
if (fullVignette) {
  MultMed[c('iid','iid.w','iid.surv')] <- NULL
  saveRDS(MultMed, "data/MultMed.rds")
} else {
  MultMed <- readRDS("data/MultMed.rds")
}

## -----------------------------------------------------------------------------
summary(MultMed)

## -----------------------------------------------------------------------------
sessionInfo()

