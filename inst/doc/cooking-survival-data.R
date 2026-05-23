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

## -----------------------------------------------------------------------------
 library(mets)
 options(warn=-1)
 set.seed(10) # to control output in simulations

## -----------------------------------------------------------------------------
 nsim <- 1000
 chaz <-  c(0,1,1.5,2,2.1)
 breaks <- c(0,10,   20,  30,   40)
 cumhaz <- cbind(breaks,chaz)
 X <- rbinom(nsim,1,0.5)
 beta <- 0.2
 rrcox <- exp(X * beta)
 
 pctime <- rchaz(cumhaz,n=nsim)
 pctimecox <- rchaz(cumhaz,rrcox)

## -----------------------------------------------------------------------------
 library(mets)
 n <- nsim
 data(bmt)
 bmt$bmi <- rnorm(408)
 dcut(bmt) <- gage~age
 data <- bmt
 cox1 <- phreg(Surv(time,cause==1)~tcell+platelet+age,data=bmt)

 dd <- sim_phreg(cox1,n,data=bmt)
 dtable(dd,~cause)
 scox1 <- phreg(Surv(time,cause==1)~tcell+platelet+age,data=dd)
 cbind(coef(cox1),coef(scox1))
 par(mfrow=c(1,1))
 plot(scox1,col=2); plot(cox1,add=TRUE,col=1)

 ## modify regression coefficients 
 cox10 <- cox1
 cox10$coef <- c(0,0.4,0.3)
 dd <- sim_phreg(cox10,n,data=bmt)
 dtable(dd,~cause)
 scox1 <- phreg(Surv(time,cause==1)~tcell+platelet+age,data=dd)
 cbind(coef(cox10),coef(scox1))
 par(mfrow=c(1,1))
 plot(scox1,col=2); plot(cox10,add=TRUE,col=1)

## -----------------------------------------------------------------------------
 data(bmt); 
 cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~tcell+platelet,data=bmt)

 X1 <- bmt[,c("tcell","platelet")]
 n <- nsim
 xid <- sample(1:nrow(X1),n,replace=TRUE)
 Z1 <- X1[xid,]
 Z2 <- X1[xid,]
 rr1 <- exp(as.matrix(Z1) %*% cox1$coef)
 rr2 <- exp(as.matrix(Z2) %*% cox2$coef)

 d <-  rcrisk(cox1$cum,cox2$cum,rr1,rr2)
 dd <- cbind(d,Z1)

 scox1 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
 scox2 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
 par(mfrow=c(1,2))
 plot(cox1); plot(scox1,add=TRUE,col=2)
 plot(cox2); plot(scox2,add=TRUE,col=2)
 cbind(cox1$coef,scox1$coef,cox2$coef,scox2$coef)

## -----------------------------------------------------------------------------
 data(sTRACE)
 dtable(sTRACE,~chf+diabetes)
 coxs <-   phreg(Surv(time,status==9)~strata(diabetes,chf),data=sTRACE)
 strata <- sample(0:3,nsim,replace=TRUE)
 simb <- sim_phreg(coxs,nsim,data=NULL,strata=strata)

 cc <-   phreg(Surv(time,status)~strata(strata),data=simb)
 plot(coxs,col=1); plot(cc,add=TRUE,col=2)

 simb1 <- sim_phreg(coxs,nsim,data=sTRACE)
 cc1 <-   phreg(Surv(time,status)~strata(diabetes,chf),data=simb1)
 plot(cc1,add=TRUE,col=3)

## ----stratified---------------------------------------------------------------
 ## competing risks with phreg 
 cox0 <- phreg(Surv(time,cause==0)~tcell+platelet,data=bmt)
 cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~strata(tcell)+platelet,data=bmt)
 coxs <- list(cox0,cox1,cox2)
 dd <- sim_phregs(coxs,n,data=bmt)

 ## verify cause-specific hazards match fitted model; increase n for better agreement
 scox0 <- phreg(Surv(time,cause==1)~tcell+platelet,data=dd)
 scox1 <- phreg(Surv(time,cause==2)~tcell+platelet,data=dd)
 scox2 <- phreg(Surv(time,cause==3)~strata(tcell)+platelet,data=dd)
 cbind(cox0$coef,scox0$coef)
 cbind(cox1$coef,scox1$coef)
 cbind(cox2$coef,scox2$coef)
 par(mfrow=c(1,3))
 plot(cox0); plot(scox0,add=TRUE,col=2); 
 plot(cox1); plot(scox1,add=TRUE,col=2); 
 plot(cox2); plot(scox2,add=TRUE,col=2); 
 
 ########################################
 ## second example 
 ########################################

 cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
 coxs <- list(cox1,cox2)
 dd <- sim_phregs(coxs,n,data=bmt)
 scox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=dd)
 scox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=dd)
 cbind(cox1$coef,scox1$coef)
 cbind(cox2$coef,scox2$coef)
 par(mfrow=c(1,2))
 plot(cox1); plot(scox1,add=TRUE); 
 plot(cox2); plot(scox2,add=TRUE); 

## -----------------------------------------------------------------------------
 library(mets)
 n <- nsim
 data(bmt)
 bmt$bmi <- rnorm(408)
 dcut(bmt) <- gage~age
 data <- bmt
 cox1 <- phreg(Surv(time,cause==1)~strata(tcell,platelet),data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~strata(gage,tcell),data=bmt)
 cox3 <- phreg(Surv(time,cause==0)~strata(platelet)+bmi,data=bmt)
 coxs <- list(cox1,cox2,cox3)

 dd <- sim_phregs(coxs,n,data=bmt,extend=0.002)
 dtable(dd,~cause)
 scox1 <- phreg(Surv(time,cause==1)~strata(tcell,platelet),data=dd)
 scox2 <- phreg(Surv(time,cause==2)~strata(gage,tcell),data=dd)
 scox3 <- phreg(Surv(time,cause==3)~strata(platelet)+bmi,data=dd)
 cbind(coef(cox1),coef(scox1), coef(cox2),coef(scox2), coef(cox3),coef(scox3))
 par(mfrow=c(1,3))
 plot(scox1,col=2); plot(cox1,add=TRUE,col=1)
 plot(scox2,col=2); plot(cox2,add=TRUE,col=1)
 plot(scox3,col=2); plot(cox3,add=TRUE,col=1)

## ----sim-phregs, label="entry-stratified"-------------------------------------
data(bmt)
cox0 <- phreg(Surv(time, cause == 0) ~ tcell + platelet+age,          data = bmt)
cox1 <- phreg(Surv(time, cause == 1) ~ tcell + platelet+age,          data = bmt)
cox2 <- phreg(Surv(time, cause == 2) ~ strata(tcell) + platelet+age,  data = bmt)

nsim <- 800
entry <- rbinom(nsim, 1, 0.5) * runif(nsim) * 60

dd <- sim_phreg(cox0,nsim, data = bmt,entry=entry)
scox0 <- phreg(Surv(entry,time, cause == 1) ~ tcell + platelet+age,         data = dd)
cbind(cox0$coef, scox0$coef)
par(mfrow = c(2, 2))
plot(cox0); plot(scox0, add = TRUE, col = 2)

dd <- sim_phregs(list(cox0, cox1, cox2), nsim, data = bmt,entry=entry)

scox0 <- phreg(Surv(entry,time, cause == 1) ~ tcell + platelet+age,         data = dd)
scox1 <- phreg(Surv(entry,time, cause == 2) ~ tcell + platelet+age,         data = dd)
scox2 <- phreg(Surv(entry,time, cause == 3) ~ strata(tcell) + platelet+age, data = dd)

cbind(cox0$coef, scox0$coef)
cbind(cox1$coef, scox1$coef)
cbind(cox2$coef, scox2$coef)

plot(cox0); plot(scox0, add = TRUE, col = 2)
plot(cox1); plot(scox1, add = TRUE, col = 2)
plot(cox2); plot(scox2, add = TRUE, col = 2)

## ----weibull1-----------------------------------------------------------------
data(sTRACE, package = "mets")
dat <- sTRACE
cox1 <- phreg(Surv(time, status > 0) ~ strata(chf) + I(age - 67), data = sTRACE)
coxw <- phreg_weibull(Surv(time, status > 0) ~ chf + age,
    shape.formula = ~chf,
    data = sTRACE
    )
coxw

tt <- seq(0, max(sTRACE$time), length.out = 100)
newd <- data.frame(chf = c(1, 0), age=67)
pr <- predict(coxw, newdata = newd, times = tt, type="chaz")
plot(cox1, col = 1)
lines(tt, pr[, 1, 1], lty=2, lwd=2)
lines(tt, pr[, 1, 2], lty = 1, lwd = 2)

## ----weibull_sim--------------------------------------------------------------
n <- 5000
newd <- mets::dsample(size=n, sTRACE[,c("chf","age")]) # bootstrap covariates
lp <- predict(coxw, newdata=newd, type="lp") # linear-predictors
head(lp)

## simulate event times
tt <- rweibullcox(nrow(lp), rate = exp(lp[,1]), shape= exp(lp[,2]))

# censoring model
censw <- phreg_weibull(Surv(time, status==0) ~ 1, data=sTRACE)
censpar <- exp(coef(censw))
censtime <- pmin(8, rweibullcox(nrow(lp), censpar[1], censpar[2]))

# combined simulated data
newd <- transform(newd, time=pmin(tt, censtime), status=(tt<=censtime))
head(newd)

# estimate weibull model on new data
phreg_weibull(Surv(time,status) ~ chf + age, ~chf, data=newd)

## -----------------------------------------------------------------------------
# simulate(coxw, n = 5, cens.model = NULL, data=newd, var.names = c("time", "status"))
simulate(coxw, nsim = 5)

## -----------------------------------------------------------------------------
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- CPH_HPN_CRBSI$mechanical
 dr2 <- scalecumhaz(dr,1.5)
 cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))

 iddata <- sim_multistate(nsim,base1,base1,dr,dr2,cens=cens)
 dlist(iddata,.~id|id<3,n=0)
  
 ### estimating rates from simulated data  
 c0 <- phreg(Surv(start,stop,status==0)~+1,iddata)
 c3 <- phreg(Surv(start,stop,status==3)~+strata(from),iddata)
 c1 <- phreg(Surv(start,stop,status==1)~+1,subset(iddata,from==2))
 c2 <- phreg(Surv(start,stop,status==2)~+1,subset(iddata,from==1))
 ###
 par(mfrow=c(2,2))
 plot(c0)
 lines(cens,col=2) 
 plot(c3,main="rates 1-> 3 , 2->3")
 lines(dr,col=1,lwd=2)
 lines(dr2,col=2,lwd=2)
 ###
 plot(c1,main="rate 1->2")
 lines(base1,lwd=2)
 ###
 plot(c2,main="rate 2->1")
 lines(base1,lwd=2)


## -----------------------------------------------------------------------------
library(mets)
nsim <- 400
rho1 <- 0.4; rho2 <- 2
beta <- c(0.3,-0.3,-0.3,0.3)

dats <- simul_cifs(nsim,rho1,rho2,beta,rc=0.5,depcens=0,type="logistic")

par(mfrow=c(1,2))
# Fitting regression model with CIF logistic-link 
cif1 <- cifreg(Event(time,status)~Z1+Z2,dats)
summary(cif1)
plot(cif1)
lines(attr(dats,"Lam1"))

dats <- simul_cifs(nsim,rho1,rho2,beta,rc=0.5,depcens=0,type="cloglog")
ciff <- cifregFG(Event(time,status)~Z1+Z2,dats)
summary(ciff)
plot(ciff)
lines(attr(dats,"Lam1"))

## -----------------------------------------------------------------------------
 data(bmt)
 ################################################################
 #  simulating several causes with specific cumulatives 
 ################################################################
 ## two logistic link models 
 cif1 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=1)
 cif2 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=2)

 dd <- sim_cifs(list(cif1,cif2),nsim,data=bmt)

 ## still logistic link 
 scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=1)
 ## 2nd cause not on logistic form due to restriction
 scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dd,cause=2)
    
 cbind(cif1$coef,scif1$coef)
 cbind(cif2$coef,scif2$coef)
 par(mfrow=c(1,2))   
 plot(cif1); plot(scif1,add=TRUE,col=2)
 plot(cif2); plot(scif2,add=TRUE,col=2)

## -----------------------------------------------------------------------------
 data(bmt)
 ## two cloglog 
 cif1 <-  cifregFG(Event(time,cause)~tcell+platelet,data=bmt,cause=1)
 cif2 <-  cifregFG(Event(time,cause)~tcell+platelet,data=bmt,cause=2)

 nsim <- 800
 entry <- rbinom(nsim,1,0.5)*runif(nsim)*60
 dd <- sim_cif(cif1,nsim,data=bmt,entry=entry)

 scif1 <-  cif(Event(entry,time,cause)~strata(tcell,platelet),data=dd,cause=1)
 plot(scif1); 
 ###
 pcif1 <- predict(cif1,expand.grid(tcell=0:1,platelet=0:1))
 plot(pcif1,add=TRUE)

## -----------------------------------------------------------------------------
 dd <- sim_cifs(list(cif1,cif2),nsim,data=bmt,entry=entry)

 ## logistic link; nonparametric Aalen-Johansen with delayed entry
 scif1 <-  cif(Event(entry,time,cause)~strata(tcell,platelet),data=dd,cause=1)
 plot(scif1); 
 ###
 pcif1 <- predict(cif1,expand.grid(tcell=0:1,platelet=0:1))
 plot(pcif1,add=TRUE)

## -----------------------------------------------------------------------------
rho1 <- 0.3; rho2 <- 4
set.seed(100)
beta=c(0.3,-0.3,-0.3,0.3)
dep=0
rc <- 0.9

n <- nsim
entry <- rbinom(n,1,0.5)*runif(n)*6
data <- simul_cifs(n,rho1,rho2,beta,bin=1,rc=0.5,rate=c(3,7),entry=entry) 
###
scif1 <- cif(Event(entry,time,status)~strata(Z1,Z2),data)
plot(scif1,ylim=c(0,0.4))
###

## and without delayed entry for comparison
data <- simul_cifs(n,rho1,rho2,beta,bin=1,rc=0.5,rate=c(3,7)) 
scif1 <- cif(Event(entry,time,status)~strata(Z1,Z2),data)
plot(scif1,add=TRUE)

## true baseline cif, cloglog link 
baset <- attr(data,"Lam1")[,2]
timet <- attr(data,"Lam1")[,1]
F1base <- 1-exp(-baset)
lines(timet,F1base,lwd=3)

## -----------------------------------------------------------------------------
library(mets)
 data(hfactioncpx12)
 hf <- hfactioncpx12
 hf$x <- as.numeric(hf$treatment) 
 n <- 200

 ##  to fit Cox  models 
 xr <- phreg(Surv(entry,time,status==1)~treatment+cluster(id),data=hf)
 dr <- phreg(Surv(entry,time,status==2)~treatment+cluster(id),data=hf)
 estimate(xr)
 estimate(dr)

 simcoxcox <- sim_recurrent_ts(xr,dr,n=n,data=hf)

 xrs <- phreg(Surv(entry,time,status==1)~treatment+cluster(id),data=simcoxcox)
 drs <- phreg(Surv(entry,time,status==3)~treatment+cluster(id),data=simcoxcox)
 estimate(xrs)
 estimate(drs)

 par(mfrow=c(1,2))
 plot(xrs); 
 plot(xr,add=TRUE)
###
 plot(drs)
 plot(dr,add=TRUE)


## -----------------------------------------------------------------------------
 recGL <- recreg(Event(entry,time,status)~treatment+cluster(id),hf,death.code=2)
 estimate(recGL)
 estimate(dr)

 simglcox <- sim_recurrent_ts(recGL,dr,n=n,data=hf)

 GLs <- recreg(Event(entry,time,status)~treatment+cluster(id),data=simglcox,death.code=3)
 drs <- phreg(Surv(entry,time,status==3)~treatment+cluster(id),data=simglcox)
 estimate(GLs)
 estimate(drs)

 par(mfrow=c(1,2))
 plot(GLs); 
 plot(recGL,add=TRUE)
 plot(drs)
 plot(dr,add=TRUE)


## -----------------------------------------------------------------------------
data(hfactioncpx12)
hf <- hfactioncpx12
hf$x <- as.numeric(hf$treatment)
hf$age <- rnorm(741)[hf$id]
hf$Z1 <- rbinom(741,1,0.5)[hf$id]

xr <- phreg(Surv(entry,time,status==1)~strata(x)+age+cluster(id),data=hf)
dr <- phreg(Surv(entry,time,status==2)~x+strata(Z1)+age+cluster(id),data=hf)
n <- 400 
rr <- sim_recurrent_ts(xr,dr,n=n,data=hf)
rxr <- phreg(Surv(entry,time,status==1)~strata(x)+age+cluster(id),data=rr)
rdr <- phreg(Surv(entry,time,status==3)~x+strata(Z1)+age+cluster(id),data=rr)
estimate(xr)
estimate(rxr)
estimate(dr)
estimate(rdr)
plot(xr); plot(rxr,add=TRUE)
plot(dr,add=TRUE,col=2); plot(rdr,add=TRUE,col=2)
###

glr <- recreg(Event(entry,time,status)~strata(x)+age+cluster(id),data=hf,death.code=2)
dr <- phreg(Surv(entry,time,status==2)~x+strata(Z1)+age+cluster(id),data=hf)
n <- 400
rr <- sim_recurrent_ts(glr,dr,n=n,data=hf)
rxr <- recreg(Event(entry,time,status)~strata(x)+age+cluster(id),data=rr,death.code=3)
rdr <- phreg(Surv(entry,time,status==3)~x+strata(Z1)+age+cluster(id),data=rr)
estimate(xr)
estimate(rxr)
estimate(dr)
estimate(rdr)
plot(glr); plot(rxr,add=TRUE)
plot(dr,add=TRUE,col=2); plot(rdr,add=TRUE,col=2)

## -----------------------------------------------------------------------------
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- CPH_HPN_CRBSI$mechanical

 n <- 400
 rr <- sim_recurrent(n,base1,death.cumhaz=dr)
 ###
 mets:::showfitsim(causes=1,rr,dr,base1,base1,which=1:2)

 rr <- sim_recurrentII(n,base1,base4,death.cumhaz=dr)
 dtable(rr,~death+status)
 mets:::showfitsim(causes=2,rr,dr,base1,base4,which=1:2)

 cumhaz <- list(base1,base1,base4)
 drl <- list(dr,base4)
 rr <- sim_recurrent_list(n,cumhaz,death.cumhaz=drl)
 dtable(rr,~death+status)
 mets:::showfitsimList(rr,cumhaz,drl) 

## -----------------------------------------------------------------------------
sessionInfo()

