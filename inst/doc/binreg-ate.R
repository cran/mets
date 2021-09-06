## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets) 

n <- 3000
kumar <- kumarsim(n,depcens=1)
kumar$cause <- kumar$status
kumar$ttt24 <- kumar[,6]
dtable(kumar,~cause)

c2 <- binregATE(Event(time,cause)~gp+dnr+preauto+ttt24,kumar,cause=2,
	treat.model=gp~dnr+preauto+ttt24,time=40,cens.model=~strata(gp,dnr,preauto,ttt24))
summary(c2)

c1 <- binregATE(Event(time,cause)~gp+dnr+preauto+ttt24,kumar,cause=2,
		treat.model=gp~dnr+preauto+ttt24,time=60)
summary(c1)

## -----------------------------------------------------------------------------

kumar$cause2 <- 1*(kumar$cause==2)

b3 <- logitATE(cause2~gp+dnr+preauto+ttt24,kumar,treat.model=gp~dnr+preauto+ttt24)
summary(b3)

## quite similar to binregATE
b2 <- logitATE(Event(time,cause)~gp+dnr+preauto+ttt24,kumar,cause=2,time=60,treat.model=gp~dnr+preauto+ttt24)
summary(b2)


## -----------------------------------------------------------------------------
sessionInfo()

