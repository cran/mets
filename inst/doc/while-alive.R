## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## ----while-alive--------------------------------------------------------------
data(hfactioncpx12)

dtable(hfactioncpx12,~status)
dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,death.code=2)
summary(dd)

dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
		   death.code=2,trans=.333)
summary(dd,type="log")


## ----while-alive-reg----------------------------------------------------------
hfactioncpx12$age <- rnorm(741)[hfactioncpx12$id]
hfactioncpx12$sex <- rbinom(741,1,0.5)[hfactioncpx12$id]

dd <- WA_reg(Event(entry,time,status)~treatment+age+sex+cluster(id),hfactioncpx12,time=2,death.code=2)
summary(dd)

## ----wa-composite-------------------------------------------------------------
hfactioncpx12$marks <- runif(nrow(hfactioncpx12))

##ddmg <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
##cause=1:2,death.code=2,marks=hfactioncpx12$marks)
##summary(ddmg)

ddm <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
cause=1:2,death.code=2,marks=hfactioncpx12$status)


## -----------------------------------------------------------------------------
sessionInfo()

