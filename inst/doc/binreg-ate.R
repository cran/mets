## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets) 
set.seed(100)
###
n <- 400
kumar <- kumarsim(n,depcens=1)
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]
dtable(kumar,~cause)
dfactor(kumar) <- gp.f~gp
kumar$id <- 1:n
kumar$idc <- sample(100,n,TRUE)
kumar$ids <- sample(n,n)
kumar$id2 <- sample(n,n)
kumar2 <- kumar[order(kumar$id2),]
kumar$int <- interaction(kumar$gp,kumar$dnr)
kumar2$int <- interaction(kumar2$gp,kumar2$dnr)
clust <- 0

b2 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
	treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(b2)

b5 <- binregATE(Event(time,cause)~int+preauto+ttt24,kumar,cause=2,
		treat.model=int~preauto+ttt24,cens.code=0,time=60)
summary(b5)

## -----------------------------------------------------------------------------
ib2 <- logitIPCWATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
	treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(ib2)

ib5 <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,
	 time=60,cens.model=~strata(gp,dnr))
summary(ib5)

ibs <- logitIPCW(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,cens.code=0,time=60)
summary(ibs)

check <- 0
if (check==1) {
require(riskRegression)
e.wglm <- wglm( regressor.event=~gp.f+dnr+preauto+ttt24, formula.censor = Surv(time,cause==0)~+1, times = 60, data = kumar, product.limit=TRUE,cause=2)
summary(e.wglm)$coef
estimate(ibs)

es.wglm <- wglm( regressor.event=~gp.f+dnr+preauto+ttt24, 
formula.censor = Surv(time,cause==0)~strata(gp,dnr), times = 60, 
data = kumar, product.limit=TRUE,cause=2)
summary(es.wglm)$coef
estimate(ib5)
}

## -----------------------------------------------------------------------------
kumar$cause2 <- 1*(kumar$cause==2)

b3 <- logitATE(cause2~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)

###library(targeted)
###b3a <- ate(cause2~gp.f|dnr+preauto+ttt24| dnr+preauto+ttt24,kumar,family=binomial)
###summary(b3a)

## calculate also relative risk
estimate(coef=b3$riskDR,vcov=b3$var.riskDR,f=function(p) p[1]/p[2])

## -----------------------------------------------------------------------------
b3 <- normalATE(time~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)

## -----------------------------------------------------------------------------
sessionInfo()

