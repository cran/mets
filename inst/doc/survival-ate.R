## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
set.seed(100)

data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
dfactor(bmt) <- tcell~tcell
bmt$event <- (bmt$cause!=0)*1

fg1 <- cifreg(Event(time,cause)~tcell+platelet+age,bmt,cause=1,
	      cox.prep=TRUE,propodds=NULL)
summary(survivalG(fg1,bmt,50))

fg2 <- cifreg(Event(time,cause)~tcell+platelet+age,bmt,cause=2,
	      cox.prep=TRUE,propodds=NULL)
summary(survivalG(fg2,bmt,50))

ss <- phreg(Surv(time,event)~tcell+platelet+age,bmt)
summary(survivalG(ss,bmt,50))


## -----------------------------------------------------------------------------

br1 <- binregATE(Event(time,cause)~tcell+platelet+age,bmt,cause=1,
		 time=40,treat.model=tcell~platelet+age)
summary(br1)

sr1 <- binregATE(Event(time,event)~tcell+platelet+age,bmt,cause=1,
		 time=40, treat.model=tcell~platelet+age)
summary(sr1)

## -----------------------------------------------------------------------------
sessionInfo()

