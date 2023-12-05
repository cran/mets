## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets) 
set.seed(100)

n <- 200
ddf <- mets:::gsim(n,covs=1,null=0,cens=1,ce=1,betac=c(0.3,1))
true <- apply(ddf$TTt<2,2,mean)
true
datat <- ddf$datat
## set-random response on data, only relevant after status==2 
response <- rbinom(n,1,0.5)
datat$response <- as.factor(response[datat$id]*datat$Count2)
datat$A000 <- as.factor(1)
datat$A111 <- as.factor(1)

bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.f~A0.f,
		augmentR1=~X11+X12+TR, augmentR0=~X01+X02,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f))
bb

estimate(coef=bb$riskG$riskG01[,1],vcov=crossprod(bb$riskG.iid$riskG01))

estimate(coef=bb$riskG$riskG01[,1],vcov=crossprod(bb$riskG.iid$riskG01),f=function(p) c(p[1]/p[2],p[3]/p[4]))
estimate(coef=bb$riskG$riskG01[,1],vcov=crossprod(bb$riskG.iid$riskG01),f=function(p) c(p[1]-p[2],p[3]-p[4]))


## -----------------------------------------------------------------------------
## 2 levels for each response , fixed weights 
datat$response.f <- as.factor(datat$response)
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),
		estpr=c(0,0),pi0=0.5,pi1=0.5)
bb

## 2 levels for each response ,  estimated treat probabilities
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),estpr=c(1,1))
bb


## 2 and 3 levels for each response , fixed weights 
datat$A1.23.f <- as.numeric(datat$A1.f)
dtable(datat,~A1.23.f+response)
datat <- dtransform(datat,A1.23.f=2+rbinom(nrow(datat),1,0.5),
		    Count2==1 & A1.23.f==2 & response==0)
dtable(datat,~A1.23.f+response)
datat$A1.23.f <- as.factor(datat$A1.23.f)
dtable(datat,~A1.23.f+response|Count2==1)
###
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.23.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),
		estpr=c(1,0),pi1=c(0.3,0.5))
bb

## 2 and 3 levels for each response , estimated 
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.23.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),estpr=c(1,1))
bb

## 2 and 1 level for each response 
datat$A1.21.f <- as.numeric(datat$A1.f)
dtable(datat,~A1.21.f+response|Count2==1)
datat <- dtransform(datat,A1.21.f=1,Count2==1 & response==1)
dtable(datat,~A1.21.f+response|Count2==1)
datat$A1.21.f <- as.factor(datat$A1.21.f)
dtable(datat,~A1.21.f+response|Count2==1)
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.21.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),estpr=c(1,1))
bb

## known weights 
bb <- binregTSR(Event(entry,time,status)~+1+cluster(id),datat,time=2,cause=c(1),response.code=2,
		treat.model0=A0.f~+1, treat.model1=A1.21.f~A0.f*response.f,
		augmentR0=~X01+X02, augmentR1=~X11+X12,
		augmentC=~X01+X02+A11t+A12t+X11+X12+TR, cens.model=~strata(A0.f),estpr=c(1,0),pi1=c(0.5,1))
bb

## -----------------------------------------------------------------------------
sessionInfo()

