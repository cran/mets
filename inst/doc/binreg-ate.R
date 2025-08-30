## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
###install.packages("mets")

library(mets) 
set.seed(100)
###
n <- 400
kumar <- kumarsim(n,depcens=1)
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]
dtable(kumar,~cause)
dfactor(kumar) <- gp.f~gp
kumar$id <- 1:400
kumar$idc <- sample(100,400,TRUE)
kumar$ids <- sample(400,400)
kumar$id2 <- sample(400,400)
kumar2 <- kumar[order(kumar$id2),]
kumar$int <- interaction(kumar$gp,kumar$dnr)
kumar2$int <- interaction(kumar2$gp,kumar2$dnr)
clust <- 0

b2 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar2,cause=2,
	treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(b2)

b5 <- binregATE(Event(time,cause)~int+preauto+ttt24,kumar,cause=2,
		treat.model=int~preauto+ttt24,cens.code=0,time=60)
summary(b5)

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

