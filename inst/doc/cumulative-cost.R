## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  #dev="svg",
  dpi=50,
  fig.width=7, fig.height=5.5,
  out.width="600px",
  fig.retina=1,
  comment = "#>"
)
fullVignette <- Sys.getenv("_R_FULL_VIGNETTE_") %in% c("1","TRUE")
library(mets)

## -----------------------------------------------------------------------------
library(mets)
data(hfactioncpx12)
hf <- hfactioncpx12
hf$severity <- abs((5+rnorm(741)*2))[hf$id]

proc_design <- mets:::proc_design
## marginal mean using formula  
outNZ <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id)
			 +marks(severity),hf,cause=1,death.code=2)
plot(outNZ,se=TRUE)
summary(outNZ,times=3) 

outN <- recurrentMarginal(Event(entry,time,status)~strata(treatment)+cluster(id),data=hf,
		       cause=1,death.code=2)
plot(outN,se=TRUE,add=TRUE)
summary(outN,times=3) 

## -----------------------------------------------------------------------------
outNZ3 <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id)+marks(severity),data=hf,
		  cause=1,death.code=2,time=3,cens.model=~strata(treatment),model="lin")
summary(outNZ3)
head(iid(outNZ3))

outN3 <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id),data=hf,cause=1,death.code=2,time=3,
		 cens.model=~strata(treatment),model="lin")
summary(outN3)
head(iid(outN3))

## -----------------------------------------------------------------------------
propNZ <- recreg(Event(entry,time,status)~treatment+marks(severity)+cluster(id),data=hf,cause=1,death.code=2)
summary(propNZ) 
plot(propNZ,main="Baselines")
     
GL <- recreg(Event(entry,time,status)~treatment+cluster(id),hf,cause=1,death.code=2)
summary(GL)
plot(GL,add=TRUE,col=2)

## -----------------------------------------------------------------------------
ooNZ <- prob.exceed.recurrent(Event(entry,time,status)~strata(treatment)+cluster(id)+marks(severity),data=hf,
			      cause=1,death.code=2,exceed=c(1,5,10,20))
plot(ooNZ,strata=1)
plot(ooNZ,strata=2,add=TRUE)
summary(ooNZ,times=3)

## -----------------------------------------------------------------------------
hf$lost5 <- 5-hf$time

RecLost <- recregIPCW(Event(entry,time,status)~-1+treatment+cluster(id)+marks(lost5),data=hf,
		   cause=1,death.code=2,time=5,cens.model=~strata(treatment),model="lin")
summary(RecLost)
head(iid(RecLost))

## -----------------------------------------------------------------------------
sessionInfo()

