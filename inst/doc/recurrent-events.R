## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets)
set.seed(1000) # to control output in simulatins for p-values below.

## -----------------------------------------------------------------------------
data(base1cumhaz)
data(base4cumhaz)
data(drcumhaz)
ddr <- drcumhaz
base1 <- base1cumhaz
base4 <- base4cumhaz
rr <- simRecurrent(1000,base1,death.cumhaz=ddr)
rr$x <- rnorm(nrow(rr)) 
rr$strata <- floor((rr$id-0.01)/500)
dlist(rr,.~id| id %in% c(1,7,9))

## -----------------------------------------------------------------------------
#  to fit non-parametric models with just a baseline 
xr <- phreg(Surv(entry,time,status)~cluster(id),data=rr)
dr <- phreg(Surv(entry,time,death)~cluster(id),data=rr)
par(mfrow=c(1,3))
bplot(dr,se=TRUE)
title(main="death")
bplot(xr,se=TRUE)
# robust standard errors 
rxr <-   robust.phreg(xr,fixbeta=1)
bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=4)

# marginal mean of expected number of recurrent events 
out <- recurrentMarginal(xr,dr)
bplot(out,se=TRUE,ylab="marginal mean",col=2)

## -----------------------------------------------------------------------------
xr <- phreg(Surv(entry,time,status)~strata(strata)+cluster(id),data=rr)
dr <- phreg(Surv(entry,time,death)~strata(strata)+cluster(id),data=rr)
par(mfrow=c(1,3))
bplot(dr,se=TRUE)
title(main="death")
bplot(xr,se=TRUE)
rxr <-   robust.phreg(xr,fixbeta=1)
bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)

out <- recurrentMarginal(xr,dr)
bplot(out,se=TRUE,ylab="marginal mean",col=1:2)

## -----------------------------------------------------------------------------
# cox case
xr <- phreg(Surv(entry,time,status)~x+cluster(id),data=rr)
dr <- phreg(Surv(entry,time,death)~x+cluster(id),data=rr)
par(mfrow=c(1,3))
bplot(dr,se=TRUE)
title(main="death")
bplot(xr,se=TRUE)
rxr <- robust.phreg(xr)
bplot(rxr,se=TRUE,robust=TRUE,add=TRUE,col=1:2)

out <- recurrentMarginal(xr,dr)
bplot(out,se=TRUE,ylab="marginal mean",col=1:2)

# predictions witout se's 
outX <- recmarg(xr,dr,Xr=1,Xd=1)
bplot(outX,add=TRUE,col=3)

## -----------------------------------------------------------------------------
###cor.mat <- corM <- rbind(c(1.0, 0.6, 0.9), c(0.6, 1.0, 0.5), c(0.9, 0.5, 1.0))
###rr <- simRecurrent(1000,base1,cumhaz2=base4,death.cumhaz=ddr)
rr <-  count.history(rr)
dtable(rr,~death+status)

oo <- prob.exceedRecurrent(rr,1)
bplot(oo)

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
with(oo,plot(time,mu,col=2,type="l"))
#
with(oo,plot(time,varN,type="l"))

## -----------------------------------------------------------------------------
rr <- simRecurrent(1000,base1,cumhaz2=base4,death.cumhaz=ddr)
rr <-  count.history(rr)
dtable(rr,~death+status)

## -----------------------------------------------------------------------------
# Bivariate probability of exceeding 
oo <- prob.exceedBiRecurrent(rr,1,2,exceed1=c(1,5,10),exceed2=c(1,2,3))
with(oo, matplot(time,pe1e2,type="s"))
nc <- ncol(oo$pe1e2)
legend("topleft",legend=colnames(oo$pe1e2),lty=1:nc,col=1:nc)

## -----------------------------------------------------------------------------
rr$strata <- 1
dtable(rr,~death+status)

covrp <- covarianceRecurrent(rr,1,2,status="status",death="death",
                        start="entry",stop="time",id="id",names.count="Count")
par(mfrow=c(1,3)) 
plot(covrp)

# with strata, each strata in matrix column, provides basis for fast Bootstrap
covrpS <- covarianceRecurrentS(rr,1,2,status="status",death="death",
        start="entry",stop="time",strata="strata",id="id",names.count="Count")

## -----------------------------------------------------------------------------
times <- seq(500,5000,500)

coo1 <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
#
mug <- Cpred(cbind(coo1$time,coo1$EN1N2),times)[,2]
mui <- Cpred(cbind(coo1$time,coo1$EIN1N2),times)[,2]
mu2.1 <- Cpred(cbind(coo1$time,coo1$mu2.1),times)[,2]
mu2.i <- Cpred(cbind(coo1$time,coo1$mu2.i),times)[,2]
mu1.2 <- Cpred(cbind(coo1$time,coo1$mu1.2),times)[,2]
mu1.i <- Cpred(cbind(coo1$time,coo1$mu1.i),times)[,2]
cbind(times,mu2.1,mu2.i)
cbind(times,mu1.2,mu1.i)

## -----------------------------------------------------------------------------
bt1 <- BootcovariancerecurrenceS(rr,1,2,status="status",start="entry",stop="time",K=100,times=times)
#bt1 <- Bootcovariancerecurrence(rr,1,2,status="status",start="entry",stop="time",K=K,times=times)

output <- list(bt1=bt1,mug=mug,mui=mui,
               bse.mug=bt1$se.mug,bse.mui=bt1$se.mui,
               dmugi=mug-mui,
	bse.dmugi=apply(bt1$EN1N2-bt1$EIN1N2,1,sd),
	mu2.1 = mu2.1 , mu2.i = mu2.i , dmu2.i=mu2.1-mu2.i,
	mu1.2 = mu1.2 , mu1.i = mu1.i , dmu1.i=mu1.2-mu1.i,
	bse.mu2.1=apply(bt1$mu2.i,1,sd), bse.mu2.1=apply(bt1$mu2.1,1,sd),
	bse.dmu2.i=apply(bt1$mu2.1-bt1$mu2.i,1,sd),
	bse.mu1.2=apply(bt1$mu1.2,1,sd), bse.mu1.i=apply(bt1$mu1.i,1,sd),
	bse.dmu1.i=apply(bt1$mu1.2-bt1$mu1.i,1,sd)
	)

## -----------------------------------------------------------------------------
tt  <- output$dmugi/output$bse.dmugi
cbind(times,2*(1-pnorm(abs(tt))))

## -----------------------------------------------------------------------------
t21  <- output$dmu1.i/output$bse.dmu1.i 
t12  <- output$dmu2.i/output$bse.dmu2.i 
cbind(times,2*(1-pnorm(abs(t21))),2*(1-pnorm(abs(t12))))

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
matplot(bt1$time,bt1$EN1N2,type="l",lwd=0.3)
matplot(bt1$time,bt1$EIN1N2,type="l",lwd=0.3)

## -----------------------------------------------------------------------------
  data(base1cumhaz)
  data(base4cumhaz)
  data(drcumhaz)
  dr <- drcumhaz
  base1 <- base1cumhaz
  base4 <- base4cumhaz

  par(mfrow=c(1,3))
  var.z <- c(0.5,0.5,0.5)
  # death related to  both causes in same way 
  cor.mat <- corM <- rbind(c(1.0, 0.0, 0.0), c(0.0, 1.0, 0.0), c(0.0, 0.0, 1.0))
  rr <- simRecurrentII(3000,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
  rr <- count.history(rr,types=1:2)
  cor(attr(rr,"z"))
  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
  plot(coo,main ="Scenario I")

## -----------------------------------------------------------------------------
  var.z <- c(0.5,0.5,0.5)
  # death related to  both causes in same way 
  cor.mat <- corM <- rbind(c(1.0, 0.0, 0.5), c(0.0, 1.0, 0.5), c(0.5, 0.5, 1.0))
  rr <- simRecurrentII(3000,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
  rr <- count.history(rr,types=1:2)
  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
  par(mfrow=c(1,3))
  plot(coo,main ="Scenario II")

## -----------------------------------------------------------------------------
  var.z <- c(0.5,0.5,0.5)
  # positive dependence for N1 and N2 all related in same way
  cor.mat <- corM <- rbind(c(1.0, 0.5, 0.5), c(0.5, 1.0, 0.5), c(0.5, 0.5, 1.0))
  rr <- simRecurrentII(3000,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
  rr <- count.history(rr,types=1:2)
  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
  par(mfrow=c(1,3))
  plot(coo,main="Scenario III")

## -----------------------------------------------------------------------------
  var.z <- c(0.5,0.5,0.5)
  # negative dependence for N1 and N2 all related in same way
  cor.mat <- corM <- rbind(c(1.0, -0.4, 0.5), c(-0.4, 1.0, 0.5), c(0.5, 0.5, 1.0))
  rr <- simRecurrentII(3000,base1,base4,death.cumhaz=dr,var.z=var.z,cor.mat=cor.mat,dependence=2)
  rr <- count.history(rr,types=1:2)
  coo <- covarianceRecurrent(rr,1,2,status="status",start="entry",stop="time")
  par(mfrow=c(1,3))
  plot(coo,main="Scenario IV")

## -----------------------------------------------------------------------------
sessionInfo()

