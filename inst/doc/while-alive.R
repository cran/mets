## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
data(hfactioncpx12)

dtable(hfactioncpx12,~status)
dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,death.code=2)
summary(dd)

dd <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
		   death.code=2,trans=.333)
summary(dd,type="log")


## -----------------------------------------------------------------------------
hfactioncpx12$marks <- runif(nrow(hfactioncpx12))

##ddmg <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
##cause=1:2,death.code=2,marks=hfactioncpx12$marks)
##summary(ddmg)

ddm <- WA_recurrent(Event(entry,time,status)~treatment+cluster(id),hfactioncpx12,time=2,
cause=1:2,death.code=2,marks=hfactioncpx12$status)


## -----------------------------------------------------------------------------
sessionInfo()

