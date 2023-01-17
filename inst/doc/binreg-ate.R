## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets) 
set.seed(100)

n <- 400
kumar <- kumarsim(n,depcens=1)
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]
dtable(kumar,~cause)
dfactor(kumar) <- gp.f~gp

### censoring model must be adapted to size of data
###c2 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
###	treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr,preauto,ttt24))
###summary(c2)

c2 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
	treat.model=gp.f~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr))
summary(c2)


c1 <- binregATE(Event(time,cause)~gp.f+dnr+preauto+ttt24,kumar,cause=2,
		treat.model=gp.f~dnr+preauto+ttt24,time=60)
summary(c1)

kumar$int <- interaction(kumar$gp,kumar$dnr)

b5 <- binregATE(Event(time,cause)~int+preauto+ttt24,kumar,cause=2,
		treat.model=int~preauto+ttt24,cens.code=0,time=60)
summary(b5)

## -----------------------------------------------------------------------------

kumar$cause2 <- 1*(kumar$cause==2)

b3 <- logitATE(cause2~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)

## -----------------------------------------------------------------------------
b3 <- normalATE(time~gp.f+dnr+preauto+ttt24,kumar,treat.model=gp.f~dnr+preauto+ttt24)
summary(b3)

## -----------------------------------------------------------------------------
sessionInfo()

