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
set.seed(1000) # to control output in simulations for p-values below.

## -----------------------------------------------------------------------------
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- CPH_HPN_CRBSI$mechanical

rr <- sim_recurrent(400,base1,death.cumhaz=dr)
rr$x <- rnorm(nrow(rr)) 
rr$strata <- floor((rr$id-0.01)/100)
dlist(rr,.~id| id %in% c(1,7,9))

## -----------------------------------------------------------------------------
#  to fit non-parametric models with just a baseline 
xr <- phreg(Surv(entry,time,status)~cluster(id),data=rr)
xdr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
par(mfrow=c(1,3))
plot(xdr,se=TRUE)
title(main="death")
plot(xr,se=TRUE)
# robust standard errors 
rxr <-   robust_phreg(xr,fixbeta=1)
plot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)

# marginal mean of expected number of recurrent events 
out <- recurrent_marginal(Event(entry,time,status)~cluster(id),data=rr,cause=1,death.code=2)
plot(out,se=TRUE,ylab="marginal mean",col=2)

## -----------------------------------------------------------------------------
summary(out,times=c(1000,2000))

## -----------------------------------------------------------------------------
xr <- phreg(Surv(entry,time,status)~strata(strata)+cluster(id),data=rr)
xdr <- phreg(Surv(entry,time,death)~strata(strata)+cluster(id),data=rr)
par(mfrow=c(1,3))
plot(xdr,se=TRUE)
title(main="death")
plot(xr,se=TRUE)
rxr <-   robust_phreg(xr,fixbeta=1)
plot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)

out <- recurrent_marginal(Event(entry,time,status)~strata(strata)+cluster(id),
			 data=rr,cause=1,death.code=2)
plot(out,se=TRUE,ylab="marginal mean",col=1:2)

## -----------------------------------------------------------------------------
test_logrankRecurrent(out)

dd <- test_marginalMean(Event(entry,time,status)~strata(strata)+cluster(id),
			 data=rr,cause=1,death.code=2)
dd
summary(dd)
dd$RAUCl

## -----------------------------------------------------------------------------
# cox case
xr <- phreg(Surv(entry,time,status)~x+cluster(id),data=rr)
xdr <- phreg(Surv(entry,time,death)~x+cluster(id),data=rr)
par(mfrow=c(1,3))
plot(xdr,se=TRUE)
title(main="death")
plot(xr,se=TRUE)
rxr <- robust_phreg(xr)
plot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)

out <- recurrentMarginalPhreg(xr,xdr)
plot(out,se=TRUE,ylab="marginal mean",col=1:2)

#### predictions witout se's 
###outX <- recmarg(xr,dr,Xr=1,Xd=1)
###plot(outX,add=TRUE,col=3)

## -----------------------------------------------------------------------------
rr <- sim_recurrent_list(100,list(base1,base1,base4),death.cumhaz=list(dr,base4),cens=3/5000,dependence=0)
dtable(rr,~status+death,level=2)
mets:::showfitsimList(rr,list(base1,base1,base4),list(dr,base4))

## -----------------------------------------------------------------------------
rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,cens=3/5000,dependence=4,var.z=1)
rr <-  count_history(rr)

rr <- transform(rr,statusD=status)
rr <- dtransform(rr,statusD=3,death==1)
dtable(rr,~statusD+status+death,level=2,response=1)

##xr <- phreg(Surv(start,stop,status==1)~cluster(id),data=rr)
##dr <- phreg(Surv(start,stop,death)~cluster(id),data=rr)
# marginal mean of expected number of recurrent events 
out <- recurrent_marginal(Event(start,stop,statusD)~cluster(id),data=rr,cause=1,death.code=3)

times <- 500*(1:10)
recEFF1 <- recurrent_marginalAIPCW(Event(start,stop,statusD)~cluster(id),data=rr,times=times,cens.code=0,
				   death.code=3,cause=1,augment.model=~Nt)
with( recEFF1, cbind(times,muP,semuP,muPAt,semuPAt,semuPAt/semuP))

times <- 500*(1:10)
###recEFF14 <- recurrent_marginalAIPCW(Event(start,stop,statusD)~cluster(id),data=rr,times=times,cens.code=0,
###death.code=3,cause=1,augment.model=~Nt+Nt2+expNt+NtexpNt)
###with(recEFF14,cbind(times,muP,semuP,muPAt,semuPAt,semuPAt/semuP))

recEFF14 <- recurrent_marginalAIPCW(Event(start,stop,statusD)~cluster(id),data=rr,times=times,cens.code=0,
death.code=3,cause=1,augment.model=~Nt+I(Nt^2)+I(exp(-Nt))+ I( Nt*exp(-Nt)))
with(recEFF14,cbind(times,muP,semuP,muPAt,semuPAt,semuPAt/semuP))

plot(out,se=TRUE,ylab="marginal mean",col=2)
k <- 1
for (t in times) {
	ci1 <- c(recEFF1$muPAt[k]-1.96*recEFF1$semuPAt[k],
  	         recEFF1$muPAt[k]+1.96*recEFF1$semuPAt[k])
	ci2 <- c(recEFF1$muP[k]-1.96*recEFF1$semuP[k],
  	         recEFF1$muP[k]+1.96*recEFF1$semuP[k])
	lines(rep(t,2)-2,ci2,col=2,lty=1,lwd=2)
	lines(rep(t,2)+2,ci1,col=1,lty=1,lwd=2)
	k <- k+1
}
legend("bottomright",c("Eff-pred"),lty=1,col=c(1,3))

## -----------------------------------------------------------------------------
n <- 200
X <- matrix(rbinom(n*2,1,0.5),n,2)
colnames(X) <- paste("X",1:2,sep="")
###
r1 <- exp( X %*% c(0.3,-0.3))
rd <- exp( X %*% c(0.3,-0.3))
rc <- exp( X %*% c(0,0))
fz <- NULL
rr <- mets:::sim_GLcox(n,base1,dr,var.z=0,r1=r1,rd=rd,rc=rc,fz,model="twostage",cens=3/5000) 
rr <- cbind(rr,X[rr$id+1,])

dtable(rr,~statusD+status+death,level=2,response=1)

times <- seq(500,5000,by=500)
recEFF1x <- recurrent_marginalAIPCW(Event(start,stop,statusD)~cluster(id),data=rr,times=times,
				   cens.code=0,death.code=3,cause=1,augment.model=~X1+X2)
with(recEFF1x, cbind(muP,muPA,muPAt,semuP,semuPA,semuPAt,semuPAt/semuP))

out <- recurrent_marginal(Event(start,stop,statusD)~cluster(id),data=rr,cause=1,death.code=3)
summary(out,times=times)

## -----------------------------------------------------------------------------
n <- 100
X <- matrix(rbinom(n*2,1,0.5),n,2)
colnames(X) <- paste("X",1:2,sep="")
###
r1 <- exp( X %*% c(0.3,-0.3))
rd <- exp( X %*% c(0.3,-0.3))
rc <- exp( X %*% c(0,0))
fz <- NULL
rr <- mets:::sim_GLcox(n,base1,dr,var.z=1,r1=r1,rd=rd,rc=rc,fz,cens=1/5000,type=2) 
rr <- cbind(rr,X[rr$id+1,])

out  <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,cens.code=0)
outs <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,cens.code=0,
		cens.model=~strata(X1,X2))
summary(out)$coef
summary(outs)$coef

## checking baseline
par(mfrow=c(1,1))
plot(out)
plot(outs,add=TRUE,col=2)
lines(scalecumhaz(base1,1),col=3,lwd=2)

## -----------------------------------------------------------------------------
outipcw  <- recregIPCW(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,
		cens.code=0,times=2000)
outipcws <- recregIPCW(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,
	    cens.code=0,times=2000,cens.model=~strata(X1,X2))
summary(outipcw)$coef
summary(outipcws)$coef

## -----------------------------------------------------------------------------
out  <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=c(1,3),
		death.code=3,cens.code=0)
summary(out)$coef

## -----------------------------------------------------------------------------
rr$binf <- rbinom(nrow(rr),1,0.5) 
rr$statusDC <- rr$statusD
rr <- dtransform(rr,statusDC=4, statusD==3 & binf==0)
rr$weight <- 1
rr <- dtransform(rr,weight=2,statusDC==3)

outC  <- recreg(Event(start,stop,statusDC)~X1+X2+cluster(id),data=rr,cause=c(1,3),
		 death.code=c(3,4),cens.code=0)
summary(outC)$coef

outCW  <- recreg(Event(start,stop,statusDC)~X1+X2+cluster(id),data=rr,cause=c(1,3),
		  death.code=c(3,4),cens.code=0,wcomp=c(1,2))
summary(outCW)$coef

plot(out,ylab="Mean composite")
plot(outC,col=2,add=TRUE)
plot(outCW,col=3,add=TRUE)

## -----------------------------------------------------------------------------
out  <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,cens.code=0)
summary(out)
baseiid <- iidBaseline(out,time=3000)
GLprediid(baseiid,rr[1:5,])

## -----------------------------------------------------------------------------
outA  <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,cause=1,death.code=3,
		 cens.code=0,augment.model=~Nt+X1+X2)
summary(outA)$coef

## -----------------------------------------------------------------------------
library(mets)
rho1 <- 1; rho2 <- 0.5
rate <- c(1,1)
tt <- seq(0, 6, by = 0.01)
base1 <- cbind(tt,rho1 * (1 - exp(-tt/rate[1])))
drcumhaz <- cbind(tt,rho2 * (1 - exp(-tt/rate[2])))
base13 <- cpred(base1,c(1,3))[,2]

dats <- mets:::sim_GLRA(100,base1,drcumhaz,varz=1)
datsA <- dats[[1]]
dtable(datsA,~statusD)
datsRA <- dats[[2]]
dtable(datsRA,~statusD)


## -----------------------------------------------------------------------------
   ## handling admin censoring with IPCW 
    outRR  <- recreg(Event(start,stop,statusD)~Z1+Z2+cluster(id), datsA, cause = 1, 
		       cens.code = c(7),death.code=3)
    estimate(outRR)

    ## Full Adm-censoring statusA, timeA
    out0  <- recreg(Event(start,stop,statusD) ~ Z1 + Z2+cluster(id) , datsA, cause = 1, 
		       cens.code = 9, death.code=3, adm.cens.time=datsA$censorA)
    estimate(out0)

## -----------------------------------------------------------------------------
    ## std Right censoring  on combined censoring time 
    outR  <- recreg(Event(start,stop,statusD)~Z1 + Z2+cluster(id) , datsRA, cause = 1, 
		       cens.code = c(0,7),death.code=3)
    estimate(outR)

    ## Combined R-IPCW + Adm-censoring status, time
    out1  <- recreg(Event(start,stop,statusD) ~ Z1 + Z2+cluster(id) , datsRA, cause = 1, 
		       cens.code = 0, death.code=3, adm.cens.time=datsRA$censorA)
    estimate(out1)

    ## Combined R-IPCW + Adm-censoring status + censoring modelling time
    out1c  <- recreg(Event(start,stop,statusD) ~ Z1 + Z2+cluster(id) , datsRA, cause = 1, 
		     cens.code = 0, death.code=3, adm.cens.time=datsRA$censorA,cens.model=~strata(Z1,Z2))
     estimate(out1c)

## -----------------------------------------------------------------------------
predR0 <- predict(outR,data.frame(Z1=0:1,Z2=0))
plot(predR0,ylim=c(0,1.5))

predR <- predict(outR,data.frame(Z1=0:1,Z2=0),se=1)
plot(predR,se=1,add=TRUE)

pred1c <- predict(out1c,data.frame(Z1=0:1,Z2=0),se=1)
plot(pred1c,se=1,add=TRUE)

## -----------------------------------------------------------------------------
km07 <- km(Event(start,stop,statusD %in% c(0,7))~strata(Z1,Z2),datsRA)
km0 <- km(Event(start,stop,statusD==0)~strata(Z1,Z2),datsRA)
plot(km07,col=1)
plot(km0,add=TRUE,col=2,lwd=2)

## -----------------------------------------------------------------------------
set.seed(100)
dr <- CPH_HPN_CRBSI$terminal
base1 <- CPH_HPN_CRBSI$crbsi 
n <- 200
X <- matrix(rbinom(n*2,1,0.5),n,2)
colnames(X) <- paste("X",1:2,sep="")
###
r1 <- exp( X %*% c(0.3,-0.3))
rd <- exp( X %*% c(0.3,-0.3))
rc <- exp( X %*% c(0,0))
fz <- NULL
## type=3 is cox-cox and type=2 is Ghosh-Lin/Cox model 
rr <- mets:::sim_GLcox(n,base1,dr,var.z=1,r1=r1,rd=rd,rc=rc,cens=1/5000,type=3) 
dtable(rr,~statusD)
rr <- cbind(rr,X[rr$id+1,])
###
out  <- phreg(Event(start,stop,statusD==1)~X1+X2+cluster(id),data=rr)
outs <- phreg(Event(start,stop,statusD==3)~X1+X2+cluster(id),data=rr)
## cox/cox
tsout <- twostageREC(outs,out,data=rr)
summary(tsout)
###
rr <- mets:::sim_GLcox(n,base1,dr,var.z=1,r1=r1,rd=rd,rc=rc,fz,cens=1/5000,type=3,share=0.5) 
rr <- cbind(rr,X[rr$id+1,])
###
out  <- phreg(Event(start,stop,statusD==1)~X1+X2+cluster(id),data=rr)
outs <- phreg(Event(start,stop,statusD==3)~X1+X2+cluster(id),data=rr)
#
tsout <- twostageREC(outs,out,data=rr,model="shared")
summary(tsout)
###
rr <- mets:::sim_GLcox(n,base1,dr,var.z=1,r1=r1,rd=rd,rc=rc,fz,cens=1/5000,type=2) 
rr <- cbind(rr,X[rr$id+1,])
outs  <- phreg(Event(start,stop,statusD==3)~X1+X2+cluster(id),data=rr)
outgl  <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),data=rr,twostage=TRUE,death.code=3)
##
## ghosh-lin/cox
glout <- twostageREC(outs,outgl,data=rr,theta=1)
summary(glout)
###
glout <- twostageREC(outs,outgl,data=rr,model="shared",theta=1,nu=0.9)
summary(glout)
glout$gradient

## -----------------------------------------------------------------------------
n <- 100
X <- matrix(rbinom(n*2,1,0.5),n,2)
colnames(X) <- paste("X",1:2,sep="")
###
r1 <- exp( X %*% c(0.3,-0.3))
rd <- exp( X %*% c(0.3,-0.3))
rc <- exp( X %*% c(0,0))
rr <- mets:::sim_GLcox(n,base1,dr,var.z=0,r1=r1,rd=rd,rc=rc,model="twostage",cens=3/5000) 
rr <- cbind(rr,X[rr$id+1,])

## -----------------------------------------------------------------------------
rr <- mets:::sim_GLcox(100,base1,dr,var.z=1,r1=r1,rd=rd,rc=rc,type=3,cens=3/5000) 
rr <- cbind(rr,X[rr$id+1,])
margsurv <- phreg(Surv(start,stop,statusD==3)~X1+X2+cluster(id),rr)
recurrent <- phreg(Surv(start,stop,statusD==1)~X1+X2+cluster(id),rr)
estimate(margsurv)
estimate(recurrent)
par(mfrow=c(1,2)); 
plot(margsurv); lines(dr,col=3); 
plot(recurrent); lines(base1,col=3)

## -----------------------------------------------------------------------------
simcoxcox <- sim_recurrent_ts(recurrent,margsurv,n=10,data=rr)

recurrentGL <- recreg(Event(start,stop,statusD)~X1+X2+cluster(id),rr,death.code=3)
simglcox <- sim_recurrent_ts(recurrentGL,margsurv,n=10,data=rr)

## -----------------------------------------------------------------------------
rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,cens=3/5000,dependence=4,var.z=1)
rr <- transform(rr,statusD=status)
rr <- dtransform(rr,statusD=3,death==1)
rr <-  count_history(rr)
dtable(rr,~statusD)

oo <- prob_exceed_recurrent(Event(entry,time,statusD)~cluster(id),rr,cause=1,death.code=3)
plot(oo,types=1:5)

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
with(oo,plot(time,meanN,col=2,type="l"))
with(oo,plot(time,varN,type="l"))

## -----------------------------------------------------------------------------
rr <- sim_recurrentII(200,base1,cumhaz2=base4,death.cumhaz=dr)
rr <-  count_history(rr)
dtable(rr,~death+status)

## -----------------------------------------------------------------------------
# Bivariate probability of exceeding 
## oo <- prob.exceedBiRecurrent(rr,1,2,exceed1=c(1,5),exceed2=c(1,2))
## with(oo, matplot(time,pe1e2,type="s"))
## nc <- ncol(oo$pe1e2)
## legend("topleft",legend=colnames(oo$pe1e2),lty=1:nc,col=1:nc)

## ----eval=FALSE---------------------------------------------------------------
#  data(CPH_HPN_CRBSI)
#  dr <- CPH_HPN_CRBSI$terminal
#  base1 <- CPH_HPN_CRBSI$crbsi
#  base4 <- CPH_HPN_CRBSI$mechanical
# 
#   par(mfrow=c(1,3))
#   var.z <- c(0.5,0.5,0.5)
#   # death related to  both causes in same way
#   cor.mat <- corM <- rbind(c(1.0, 0.0, 0.0), c(0.0, 1.0, 0.0), c(0.0, 0.0, 1.0))
#   rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
#   rr <- count_history(rr,types=1:2)
# ###  cor(attr(rr,"z"))
# ###  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
# ###  plot(coo,main ="Scenario I")

## ----eval=FALSE---------------------------------------------------------------
#   var.z <- c(0.5,0.5,0.5)
#   # death related to  both causes in same way
#   cor.mat <- corM <- rbind(c(1.0, 0.0, 0.5), c(0.0, 1.0, 0.5), c(0.5, 0.5, 1.0))
#   rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
#   rr <- count_history(rr,types=1:2)
# ###  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
# ###  par(mfrow=c(1,3))
# ###  plot(coo,main ="Scenario II")

## ----eval=FALSE---------------------------------------------------------------
#   var.z <- c(0.5,0.5,0.5)
#   # positive dependence for N1 and N2 all related in same way
#   cor.mat <- corM <- rbind(c(1.0, 0.5, 0.5), c(0.5, 1.0, 0.5), c(0.5, 0.5, 1.0))
#   rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
#   rr <- count_history(rr,types=1:2)
# ###  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
# ###  par(mfrow=c(1,3))
# ###  plot(coo,main="Scenario III")

## ----eval=FALSE---------------------------------------------------------------
#   var.z <- c(0.5,0.5,0.5)
#   # negative dependence for N1 and N2 all related in same way
#   cor.mat <- corM <- rbind(c(1.0, -0.4, 0.5), c(-0.4, 1.0, 0.5), c(0.5, 0.5, 1.0))
#   rr <- sim_recurrentII(200,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
#   rr <- count_history(rr,types=1:2)
# ###  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
# ###  par(mfrow=c(1,3))
# ###  plot(coo,main="Scenario IV")

## -----------------------------------------------------------------------------
sessionInfo()

