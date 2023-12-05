## ----include = FALSE----------------------------------------------------------
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

## survival situation
sr1 <- binregATE(Event(time,event)~tcell+platelet+age,bmt,cause=1,
		 time=40, treat.model=tcell~platelet+age)
summary(sr1)

## relative risk effect 
estimate(coef=sr1$riskDR,vcov=sr1$var.riskDR,f=function(p) p[2]/p[1],null=1)

## competing risks 
br1 <- binregATE(Event(time,cause)~tcell+platelet+age,bmt,cause=1,
		 time=40,treat.model=tcell~platelet+age)
summary(br1)


## -----------------------------------------------------------------------------
br1 <- binreg(Event(time,cause)~tcell+platelet+age,bmt,cause=1,time=40)
Gbr1 <- binregG(br1,data=bmt)
summary(Gbr1)

## contrasting average age to +2-sd age, Avalues
Gbr2 <- binregG(br1,data=bmt,varname="age",Avalues=c(0,2))
summary(Gbr2)

## -----------------------------------------------------------------------------
sessionInfo()

