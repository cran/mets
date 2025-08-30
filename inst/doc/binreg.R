## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in random noise just below.
 data(bmt)
 bmt$time <- bmt$time+runif(nrow(bmt))*0.01

 # logistic regresion with IPCW binomial regression 
 out <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50)
 summary(out)

## -----------------------------------------------------------------------------
 predict(out,data.frame(tcell=c(0,1),platelet=c(1,1)),se=TRUE)

## -----------------------------------------------------------------------------
 outs <- binreg(Event(time,cause)~tcell+platelet,bmt,time=50,cens.model=~strata(tcell,platelet))
 summary(outs)

## -----------------------------------------------------------------------------
 outs <- binreg(Event(time,cause)~tcell,bmt,time=50,cens.model=~strata(tcell))
 summary(outs)

## -----------------------------------------------------------------------------
ps <-  predict(outs,data.frame(tcell=c(0,1)),se=TRUE)
ps
sum( c(1,-1) * ps[,1])

## -----------------------------------------------------------------------------
dd <- data.frame(tcell=c(0,1))
p <- predict(outs,dd)

riskdifratio <- function(p,contrast=c(1,-1)) {
   outs$coef <- p
   p <- predict(outs,dd)[,1]
   pd <- sum(contrast*p)
   r1 <- p[1]/p[2]
   r2 <- p[2]/p[1]
   return(c(pd,r1,r2))
}
     
estimate(outs,f=riskdifratio,dd,null=c(0,1,1))

## -----------------------------------------------------------------------------
run <- 0
if (run==1) {
library(prodlim)
pl <- prodlim(Hist(time,cause)~tcell,bmt)
spl <- summary(pl,times=50,asMatrix=TRUE)
spl
}

## -----------------------------------------------------------------------------
 data(bmt)
 dcut(bmt,breaks=2) <- ~age 
 out1<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
			  strata(platelet,agecat.2),data=bmt,cause=1,time=40)
 summary(out1)

 out2<-BinAugmentCifstrata(Event(time,cause)~platelet+agecat.2+
     strata(platelet,agecat.2)+strataC(platelet),data=bmt,cause=1,time=40)
 summary(out2)

## -----------------------------------------------------------------------------
sessionInfo()

