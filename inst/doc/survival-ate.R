## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets)
set.seed(100)

data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
dfactor(bmt) <- tcell~tcell
bmt$event <- (bmt$cause!=0)*1

fg1 <- cifregFG(Event(time,cause)~tcell+platelet+age,bmt,cause=1)
summary(survivalG(fg1,bmt,time=50))

fg2 <- cifregFG(Event(time,cause)~tcell+platelet+age,bmt,cause=2)
summary(survivalG(fg2,bmt,time=50))

cif1time <- survivalGtime(fg1,bmt)
plot(cif1time,type="risk"); 

## -----------------------------------------------------------------------------
ss <- phreg(Surv(time,event)~tcell+platelet+age,bmt)
sss <- survivalG(ss,bmt,time=50)
summary(sss)

Gtime <- survivalGtime(ss,bmt)
plot(Gtime)

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

