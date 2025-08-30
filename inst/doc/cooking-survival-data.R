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
 nsim <- 200
 chaz <-  c(0,1,1.5,2,2.1)
 breaks <- c(0,10,   20,  30,   40)
 cumhaz <- cbind(breaks,chaz)
 X <- rbinom(nsim,1,0.5)
 beta <- 0.2
 rrcox <- exp(X * beta)
 
 pctime <- rchaz(cumhaz,n=nsim)
 pctimecox <- rchaz(cumhaz,rrcox)

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
 simb <- sim.base(coxs$cumhaz,nsim,stratajump=coxs$strata.jumps,strata=strata)
 cc <-   phreg(Surv(time,status)~strata(strata),data=simb)
 plot(coxs,col=1); plot(cc,add=TRUE,col=2)

## -----------------------------------------------------------------------------
 cox <-  survival::coxph(Surv(time,status==9)~vf+chf+wmi,data=sTRACE)
 sim1 <- sim.cox(cox,nsim,data=sTRACE)
 cc <- survival::coxph(Surv(time,status)~vf+chf+wmi,data=sim1)
 cbind(cox$coef,cc$coef)
 cor(sim1[,c("vf","chf","wmi")])
 cor(sTRACE[,c("vf","chf","wmi")])
 
 cox <-  phreg(Surv(time, status==9)~vf+chf+wmi,data=sTRACE)
 sim3 <- sim.cox(cox,nsim,data=sTRACE)
 cc <-  phreg(Surv(time, status)~vf+chf+wmi,data=sim3)
 cbind(cox$coef,cc$coef)
 plot(cox,se=TRUE); plot(cc,add=TRUE,col=2)
 
 coxs <-  phreg(Surv(time,status==9)~strata(chf,vf)+wmi,data=sTRACE)
 sim3 <- sim.phreg(coxs,nsim,data=sTRACE)
 cc <-   phreg(Surv(time, status)~strata(chf,vf)+wmi,data=sim3)
 cbind(coxs$coef,cc$coef)
 plot(coxs,col=1); plot(cc,add=TRUE,col=2)

## -----------------------------------------------------------------------------
 data(bmt)
 # coxph          
 cox1 <- survival::coxph(Surv(time,cause==1)~tcell+platelet,data=bmt)
 cox2 <- survival::coxph(Surv(time,cause==2)~tcell+platelet,data=bmt)
 coxs <- list(cox1,cox2)
 dd <- sim.cause.cox(coxs,nsim,data=bmt)
 scox1 <- survival::coxph(Surv(time,status==1)~tcell+platelet,data=dd)
 scox2 <- survival::coxph(Surv(time,status==2)~tcell+platelet,data=dd)
 cbind(cox1$coef,scox1$coef)
 cbind(cox2$coef,scox2$coef)

## -----------------------------------------------------------------------------
 ## stratified with phreg 
 cox0 <- phreg(Surv(time,cause==0)~tcell+platelet,data=bmt)
 cox1 <- phreg(Surv(time,cause==1)~tcell+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~strata(tcell)+platelet,data=bmt)
 coxs <- list(cox0,cox1,cox2)
 dd <- sim.cause.cox(coxs,nsim,data=bmt)
 scox0 <- phreg(Surv(time,status==1)~tcell+platelet,data=dd)
 scox1 <- phreg(Surv(time,status==2)~tcell+platelet,data=dd)
 scox2 <- phreg(Surv(time,status==3)~strata(tcell)+platelet,data=dd)
 cbind(cox0$coef,scox0$coef)
 cbind(cox1$coef,scox1$coef)
 cbind(cox2$coef,scox2$coef)
 par(mfrow=c(1,3))
 plot(cox0); plot(scox0,add=TRUE,col=2); 
 plot(cox1); plot(scox1,add=TRUE,col=2); 
 plot(cox2); plot(scox2,add=TRUE,col=2); 
 
 cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~tcell+strata(platelet),data=bmt)
 coxs <- list(cox1,cox2)
 dd <- sim.cause.cox(coxs,nsim,data=bmt)
 scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
 scox2 <- phreg(Surv(time,status==2)~tcell+strata(platelet),data=dd)
 cbind(cox1$coef,scox1$coef)
 cbind(cox2$coef,scox2$coef)
 par(mfrow=c(1,2))
 plot(cox1); plot(scox1,add=TRUE); 
 plot(cox2); plot(scox2,add=TRUE); 

## -----------------------------------------------------------------------------
 library(mets)
 n <- 100
 data(bmt)
 bmt$bmi <- rnorm(408)
 dcut(bmt) <- gage~age
 data <- bmt
 cox1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
 cox2 <- phreg(Surv(time,cause==2)~strata(gage)+tcell+platelet,data=bmt)
 cox3 <- phreg(Surv(time,cause==0)~strata(platelet)+bmi,data=bmt)
 coxs <- list(cox1,cox2,cox3)

 dd <- sim.phregs(coxs,n,data=bmt,extend=0.002)
 scox1 <- phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
 scox2 <- phreg(Surv(time,status==2)~strata(gage)+tcell+platelet,data=dd)
 scox3 <- phreg(Surv(time,status==3)~strata(platelet)+bmi,data=dd)
 cbind(coef(cox1),coef(scox1), coef(cox2),coef(scox2), coef(cox3),coef(scox3))
 par(mfrow=c(1,3))
 plot(scox1,col=2); plot(cox1,add=TRUE,col=1)
 plot(scox2,col=2); plot(cox2,add=TRUE,col=1)
 plot(scox3,col=2); plot(cox3,add=TRUE,col=1)

 coxs1 <- phreg(Surv(time,cause==1)~strata(tcell),data=bmt)
 dd <- sim.phreg(coxs1,n,data=bmt)
 scoxs1 <-  phreg(Surv(time,status==1)~strata(tcell),data=dd)
 ###
 plot(coxs1)
 plot(scoxs1,add=TRUE)

 coxs1 <- phreg(Surv(time,cause==1)~strata(tcell)+platelet,data=bmt)
 dd <- sim.phreg(coxs1,n,data=bmt)
 scoxs1 <-  phreg(Surv(time,status==1)~strata(tcell)+platelet,data=dd)
 ###
 plot(coxs1)
 plot(scoxs1,add=TRUE)


## -----------------------------------------------------------------------------
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- CPH_HPN_CRBSI$mechanical
 dr2 <- scalecumhaz(dr,1.5)
 cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))

 iddata <- simMultistate(nsim,base1,base1,dr,dr2,cens=cens)
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
cif1 <- cbind(c(0,10,20,100),c(0,0.1,0.15,0.2))
cif2 <- cbind(c(0,10,20,100),c(0,0.4,0.45,0.5))

n <- 100; lrr1=c(0.2,0.1); lrr2=c(0.2,0.1); cens=NULL
### A binary, L binary
A <- rbinom(n,1,0.5)
L <- rbinom(n,1,0.5)
###
rr1 <- exp(cbind(A,L) %*% lrr1)
rr2 <- exp(cbind(A,L) %*% lrr2)
## model is fine
mmm<-max(rr1)*max(cif1[,2])+max(rr2)*max(cif2[,2])
mcif1 <- max(cif1[,2])
mcif2 <- max(cif2[,2])
if (mmm>1) warning(" models not satisfying sum <=1\n")
### here log-link model 
T1 <- simsubdist(cif1,rr1,type="cif")
T2 <- simsubdist(cif2,rr2,type="cif")
###
dies <- rbinom(n,1,rr1*mcif1+rr2*mcif2)
sel1 <- rbinom(n,1,mcif2/(mcif1+mcif2))+1
epsilon  <- dies*(sel1)
T1$epsilon <- epsilon
###
T1$A <- A; T1$L <- L
## times given 
T1$time <- T1$timecause
T1$time2 <- T2$timecause
T1$status <- epsilon
T1 <- dtransform(T1,time=100,epsilon==0)
T1 <- dtransform(T1,status=0,epsilon==0)
###
T1 <- dtransform(T1,time=time2,epsilon==2)
T1 <- dtransform(T1,status=2,epsilon==2)

dtable(T1,~status)

par(mfrow=c(1,2))
lrr1=c(0.2,0.1);lrr2=c(0.2,0.1)
pcif1 <- cif(Event(time,status)~strata(A,L),T1,cause=1)
pcif2 <- cif(Event(time,status)~strata(A,L),T1,cause=2)
###
newd <- data.frame(expand.grid(A=0:1,L=0:1))
rr1 <- c(exp(as.matrix(newd) %*% lrr1))
rr2 <- c(exp(as.matrix(newd) %*% lrr2))
###
cifm1 <- cbind(cif1[,1],cif1[,2] %o% rr1)
cifm2 <- cbind(cif2[,1],cif2[,2] %o% rr2)
###
par(mfrow=c(1,2))
plot(pcif1,ylim=c(0,0.3)); 
matlines(cifm1[,1],cifm1[,-1],col=1,lwd=2)
###
plot(pcif2,ylim=c(0,0.7))
matlines(cifm2[,1],cifm2[,-1],col=1,lwd=2)

## -----------------------------------------------------------------------------
 data(bmt)
 ################################################################
 #  simulating several causes with specific cumulatives 
 ################################################################
 cif1 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=1)
 cif2 <-  cifreg(Event(time,cause)~tcell+age,data=bmt,cause=2)

 ## dd <- sim.cifs(list(cif1,cif2),nsim,data=bmt)
 dds <- sim.cifsRestrict(list(cif1,cif2),nsim,data=bmt)

 scif1 <-  cifreg(Event(time,cause)~tcell+age,data=dds,cause=1)
 scif2 <-  cifreg(Event(time,cause)~tcell+age,data=dds,cause=2)
    
 cbind(cif1$coef,scif1$coef)
 cbind(cif2$coef,scif2$coef)
 par(mfrow=c(1,2))   
 plot(cif1); plot(scif1,add=TRUE,col=2)
 plot(cif2); plot(scif2,add=TRUE,col=2)

## -----------------------------------------------------------------------------
 set.seed(100)
 rho1 <- 0.2; rho2 <- 10
 n <- nsim
 beta=c(0.0,-0.1,-0.5,0.3)
 dats <- simul.cifs(n,rho1,rho2,beta,rc=0.2)
 dtable(dats,~status)
 dsort(dats) <- ~time
 fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
 summary(fg)

## -----------------------------------------------------------------------------
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- CPH_HPN_CRBSI$mechanical

 n <- 100
 rr <- simRecurrent(n,base1,death.cumhaz=dr)
 ###
 par(mfrow=c(1,3))
 showfitsim(causes=1,rr,dr,base1,base1,which=1:2)

 rr <- simRecurrentII(n,base1,base4,death.cumhaz=dr)
 dtable(rr,~death+status)
 par(mfrow=c(2,2))
 showfitsim(causes=2,rr,dr,base1,base4,which=1:2)

 cumhaz <- list(base1,base1,base4)
 drl <- list(dr,base4)
 rr <- simRecurrentList(n,cumhaz,death.cumhaz=drl)
 dtable(rr,~death+status)
 showfitsimList(rr,cumhaz,drl) 

## -----------------------------------------------------------------------------

 data(hfactioncpx12)
 hf <- hfactioncpx12
 hf$x <- as.numeric(hf$treatment) 
 n <- 100

 ##  to fit non-parametric models with just a baseline 
 xr <- phreg(Surv(entry,time,status==1)~cluster(id),data=hf)
 dr <- phreg(Surv(entry,time,status==2)~cluster(id),data=hf)

 simcoxcox <- sim.recurrent(xr,dr,n=n,data=hf)

 recGL <- recreg(Event(entry,time,status)~+cluster(id),hf,death.code=2)
 simglcox <- sim.recurrent(recGL,dr,n=n,data=hf)

## -----------------------------------------------------------------------------
sessionInfo()

