% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim-pc-hazard.R
\name{sim.cif}
\alias{sim.cif}
\alias{sim.cifs}
\alias{subdist}
\alias{pre.cifs}
\alias{sim.cifsRestrict}
\alias{simsubdist}
\alias{invsubdist}
\title{Simulation of output from Cumulative incidence regression model}
\usage{
sim.cif(cif,n,data=NULL,Z=NULL,drawZ=TRUE,cens=NULL,rrc=NULL,cumstart=c(0,0),...)
}
\arguments{
\item{cif}{output form prop.odds.subdist or ccr (cmprsk), can also call invsubdist with 
with cumulative and linear predictor}

\item{n}{number of simulations.}

\item{data}{to extract covariates for simulations (draws from observed
covariates).}

\item{Z}{to use these covariates for simulation rather than drawing new ones.}

\item{drawZ}{to random sample from Z or not}

\item{cens}{specifies censoring model, if "is.matrix" then uses cumulative
hazard given, if "is.scalar" then uses rate for exponential, and if not
given then takes average rate of in simulated data from cox model.}

\item{rrc}{possible vector of relative risk for cox-type censoring.}

\item{cumstart}{to start cumulatives at time 0 in 0.}

\item{...}{arguments for invsubdist}
}
\description{
Simulates data that looks like fit from fitted cumulative incidence model
}
\examples{
data(bmt)

scif <-  cifreg(Event(time,cause)~tcell+platelet+age,data=bmt,cause=1,prop=NULL)
summary(scif)  
plot(scif)
################################################################
#  simulating several causes with specific cumulatives 
################################################################

cif1 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=1,prop=NULL)
cif2 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=2,prop=NULL)
# dd <- sim.cifsRestrict(list(cif1,cif2),200,data=bmt)
dd <- sim.cifs(list(cif1,cif2),200,data=bmt)
scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=1)
scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=2)
   
par(mfrow=c(1,2))   
plot(cif1); plot(scif1,add=TRUE,col=2)
plot(cif2); plot(scif2,add=TRUE,col=2)
}
\author{
Thomas Scheike
}
\keyword{survival}
