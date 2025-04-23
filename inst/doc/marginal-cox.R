## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  ##dev="png",
  dpi=50,
  fig.width=7.15, fig.height=5.5,
  out.width="600px",
  fig.retina=1,
  comment = "#>"  
)
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in simulatins for p-values below.
 n <- 1000
 k <- 5
 theta <- 2
 data <- simClaytonOakes(n,k,theta,0.3,3)

## -----------------------------------------------------------------------------
   out <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
   summary(out)
   # robust standard errors attached to output
   rob <- robust.phreg(out)

## -----------------------------------------------------------------------------
   # making iid decomposition of regression parameters
   betaiid <- IC(out)
   head(betaiid)
   # robust standard errors
   crossprod(betaiid/NROW(betaiid))^.5
   # same as 

## -----------------------------------------------------------------------------
  plot(rob,se=TRUE,robust=TRUE,col=3)

## -----------------------------------------------------------------------------
  pp <-  predict(out,data[1:20,],se=TRUE,robust=TRUE)
  plot(pp,se=TRUE,whichx=1:10)

## -----------------------------------------------------------------------------
tt <- twostageMLE(out,data=data)
summary(tt)

## -----------------------------------------------------------------------------
gout <- gof(out)
gout

## -----------------------------------------------------------------------------
  plot(gout)

## -----------------------------------------------------------------------------
 out <- phreg(Surv(time,status)~x+strata(cluster),data=data)
 summary(out)

## -----------------------------------------------------------------------------
sessionInfo()

