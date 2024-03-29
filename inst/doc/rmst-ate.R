## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
set.seed(100)

     data(bmt); bmt$time <- bmt$time+runif(nrow(bmt))*0.001
     # E( min(T;t) | X ) = exp( a+b X) with IPCW estimation 
     out <- resmeanIPCW(Event(time,cause!=0)~tcell+platelet+age,bmt,
                     time=50,cens.model=~strata(platelet),model="exp")
     summary(out)
     
      ### same as Kaplan-Meier for full censoring model 
     bmt$int <- with(bmt,strata(tcell,platelet))
     out <- resmeanIPCW(Event(time,cause!=0)~-1+int,bmt,time=30,
                                  cens.model=~strata(platelet,tcell),model="lin")
     estimate(out)
     out1 <- phreg(Surv(time,cause!=0)~strata(tcell,platelet),data=bmt)
     rm1 <- resmean.phreg(out1,times=30)
     summary(rm1)
     
     ## competing risks years-lost for cause 1  
     out <- resmeanIPCW(Event(time,cause)~-1+int,bmt,time=30,cause=1,
                                 cens.model=~strata(platelet,tcell),model="lin")
     estimate(out)
     ## same as integrated cumulative incidence 
     rmc1 <- cif.yearslost(Event(time,cause)~strata(tcell,platelet),data=bmt,times=30,cause=1)
     summary(rmc1)

     ## plotting the years lost for different horizon's and the two causes 
     par(mfrow=c(1,3))
     plot(rm1,years.lost=TRUE,se=1)
     ## cause refers to column of cumhaz for the different causes
     plot(rmc1,cause=1,se=1)
     plot(rmc1,cause=2,se=1)

## -----------------------------------------------------------------------------

 dfactor(bmt) <- tcell~tcell
 bmt$event <- (bmt$cause!=0)*1
 out <- resmeanATE(Event(time,event)~tcell+platelet,data=bmt,time=40,treat.model=tcell~platelet)
 summary(out)
 
 out1 <- resmeanATE(Event(time,cause)~tcell+platelet,data=bmt,cause=1,outcome="rmst-cause",
		    time=40,treat.model=tcell~platelet)
 summary(out1)

 out2 <- resmeanATE(Event(time,cause)~tcell+platelet,data=bmt,cause=2,outcome="rmst-cause",
		    time=40,treat.model=tcell~platelet)
 summary(out2)

## -----------------------------------------------------------------------------
sessionInfo()

