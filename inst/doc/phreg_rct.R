## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets) 
set.seed(100)

## Lu, Tsiatis simulation
data <- mets:::simLT(0.7,100)
dfactor(data) <- Z.f~Z
 
out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X)
summary(out)
###out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~X)
###out <- phreg_rct(Surv(time,status)~Z.f,data=data,augmentR0=~X,augmentC=~factor(Z):X,cens.model=~+1)

## -----------------------------------------------------------------------------
###library(speff2trial) 
library(mets)
data(ACTG175)
###
data <- ACTG175[ACTG175$arms==0 | ACTG175$arms==1, ]
data <- na.omit(data[,c("days","cens","arms","strat","cd40","cd80","age")])
data$days <- data$days+runif(nrow(data))*0.01
dfactor(data) <- arms.f~arms
notrun <- 1

if (notrun==0) { 
fit1 <- speffSurv(Surv(days,cens)~cd40+cd80+age,data=data,trt.id="arms",fixed=TRUE)
summary(fit1)
}
# 
# Treatment effect
#             Log HR       SE   LowerCI   UpperCI           p
# Prop Haz  -0.70375  0.12352  -0.94584  -0.46165  1.2162e-08
# Speff     -0.72430  0.12051  -0.96050  -0.48810  1.8533e-09

out <- phreg_rct(Surv(days,cens)~arms.f,data=data,augmentR0=~cd40+cd80+age,augmentC=~cd40+cd80+age)
summary(out)

## -----------------------------------------------------------------------------
dtable(data,~strat+arms)
dfactor(data) <- strat.f~strat
out <- phreg_rct(Surv(days,cens)~arms.f,data=data,augmentR0=~strat.f)
summary(out)

## -----------------------------------------------------------------------------
data(calgb8923)
calgt <- calgb8923

tm=At.f~factor(Count2)+age+sex+wbc
tm=At.f~factor(Count2)
tm=At.f~factor(Count2)*A0.f

head(calgt)
ll0 <- phreg_IPTW(Event(start,time,status==1)~strata(A0,A10)+cluster(id),calgt,treat.model=tm)
pll0 <- predict(ll0,expand.grid(A0=0:1,A10=0,id=1))
ll1 <- phreg_IPTW(Event(start,time,status==1)~strata(A0,A11)+cluster(id),calgt,treat.model=tm)
pll1 <- predict(ll1,expand.grid(A0=0:1,A11=1,id=1))
plot(pll0,se=1,lwd=2,col=1:2,lty=1,xlab="time (months)",xlim=c(0,30))
plot(pll1,add=TRUE,col=3:4,se=1,lwd=2,lty=1,xlim=c(0,30))
abline(h=0.25)
legend("topright",c("A1B1","A2B1","A1B2","A2B2"),col=c(1,2,3,4),lty=1)

summary(pll1,times=1:10)
summary(pll0,times=1:10)

## -----------------------------------------------------------------------------
library(mets)
data(calgb8923)
calgt <- calgb8923
calgt$treatvar <- 1

## making time-dependent indicators of going to B1/B2
calgt$A10t <- calgt$A11t <- 0
calgt <- dtransform(calgt,A10t=1,A1==0 & Count2==1)
calgt <- dtransform(calgt,A11t=1,A1==1 & Count2==1)
calgt0 <- subset(calgt,A0==0)

ss0 <- phreg_rct(Event(start,time,status)~A10t+A11t+cluster(id),data=subset(calgt,A0==0),
	 typesR=c("non","R1"),typesC=c("non","dynC"),
	 treat.var="treatvar",treat.model=At.f~factor(Count2),
	 augmentR1=~age+wbc+sex+TR,augmentC=~age+wbc+sex+TR+Count2)
summary(ss0)

ss1 <- phreg_rct(Event(start,time,status)~A10t+A11t+cluster(id),data=subset(calgt,A0==1),
	 typesR=c("non","R1"),typesC=c("non","dynC"),
	 treat.var="treatvar",treat.model=At.f~factor(Count2),
	 augmentR1=~age+wbc+sex+TR,augmentC=~age+wbc+sex+TR+Count2)
summary(ss1)

## -----------------------------------------------------------------------------
ssf <- phreg_rct(Event(start,time,status)~A0.f+A10t+A11t+cluster(id),
	 data=calgt,
	 typesR=c("non","R0","R1","R01"),typesC=c("non","C","dynC"),
	 treat.var="treatvar",treat.model=At.f~factor(Count2),
	 augmentR0=~age+wbc+sex,augmentR1=~age+wbc+sex+TR,augmentC=~age+wbc+sex+TR+Count2)

## -----------------------------------------------------------------------------
n <- 1000
beta <- 0.15; 
 data(CPH_HPN_CRBSI)
 dr <- CPH_HPN_CRBSI$terminal
 base1 <- CPH_HPN_CRBSI$crbsi 
 base4 <- scalecumhaz(CPH_HPN_CRBSI$mechanical,0.5)

cens <- rbind(c(0,0),c(2000,0.5),c(5110,3))
ce <- 3; betao1 <- 0

varz <- 1; dep=4; X <- z <- rgamma(n,1/varz)*varz
Z0 <- NULL
px <- 0.5
if (betao1!=0) px <- lava::expit(betao1*X)      
A0 <- rbinom(n,1,px)
r1 <- exp(A0*beta[1])
rd <- exp( A0 * 0.15)
rc <- exp( A0 * 0 )
###
rr <-    mets:::simLUCox(n,base1,death.cumhaz=dr,r1=r1,Z0=X,dependence=dep,var.z=varz,cens=ce/5000)
rr$A0 <- A0[rr$id]
rr$z1 <- attr(rr,"z")[rr$id]
rr$lz1 <- log(rr$z1)
rr$X <- rr$lz1 
rr$lX <- rr$z1
rr$statusD <- rr$status
rr <- dtransform(rr,statusD=2,death==1)
rr <- count.history(rr)
rr$Z <- rr$A0
data <- rr
data$Z.f <- as.factor(data$Z)
data$treattime <- 0
data <- dtransform(data,treattime=1,lbnr__id==1)
dlist(data,start+stop+statusD+A0+z1+treattime+Count1~id|id %in% c(4,5))

## -----------------------------------------------------------------------------
fit2 <- phreg_rct(Event(start,stop,statusD)~Z.f+cluster(id),data=data,
	 treat.var="treattime",typesR=c("non","R0"),typesC=c("non","C","dynC"),
	 augmentR0=~z1,augmentC=~z1+Count1)
summary(fit2)

## -----------------------------------------------------------------------------
n <- 500
beta=c(0.3,0.3);betatr=0.3;betac=0;betao=0;betao1=0;ce=3;fixed=1;sim=1;dep=4;varz=1;ztr=0; ce <- 3
## take possible frailty 
Z0 <- rgamma(n,1/varz)*varz
px0 <- 0.5; if (betao!=0) px0 <- expit(betao*Z0)
A0 <- rbinom(n,1,px0)
r1 <- exp(A0*beta[1])
#
px1 <- 0.5; if (betao1!=0) px1 <- expit(betao1*Z0)
A1 <- rbinom(n,1,px1)
r2 <- exp(A1*beta[2])
rtr <- exp(A0*betatr[1])
rr <-  mets:::simLUCox(n,base1,death.cumhaz=dr,cumhaz2=base1,rtr=rtr,betatr=0.3,A0=A0,Z0=Z0,
		r1=r1,r2=r2,dependence=dep,var.z=varz,cens=ce/5000,ztr=ztr)
rr$z1 <- attr(rr,"z")[rr$id]
rr$A1 <- A1[rr$id]
rr$A0 <- A0[rr$id]
rr$lz1 <- log(rr$z1)
rr <- count.history(rr,types=1:2)
rr$A1t <- 0
rr <- dtransform(rr,A1t=A1,Count2==1) 
rr$At.f <- rr$A0
rr$A0.f <- factor(rr$A0)
rr$A1.f <- factor(rr$A1)
rr <- dtransform(rr, At.f = A1, Count2 == 1)
rr$At.f <- factor(rr$At.f)
dfactor(rr)  <-  A0.f~A0
rr$treattime <- 0
rr <- dtransform(rr,treattime=1,lbnr__id==1)
rr$lagCount2 <- dlag(rr$Count2)
rr <- dtransform(rr,treattime=1,Count2==1 & (Count2!=lagCount2))
dlist(rr,start+stop+statusD+A0+A1+A1t+At.f+Count2+z1+treattime+Count1~id|id %in% c(5,10))

## -----------------------------------------------------------------------------
sse <- phreg_rct(Event(start,time,statusD)~A0.f+A1t+cluster(id),data=rr,
	 typesR=c("non","R0","R1","R01"),typesC=c("non","C","dynC"),treat.var="treattime",
	 treat.model=At.f~factor(Count2),
	 augmentR0=~z1,augmentR1=~z1,augmentC=~z1+Count1+A1t)
summary(sse)

## -----------------------------------------------------------------------------
fit2 <- phreg_rct(Event(start,stop,statusD)~Z.f+cluster(id),data=data,
	 typesR=c("non","R0"),typesC=c("non","C","dynC"),
         RCT=FALSE,treat.model=Z.f~z1,augmentR0=~z1,augmentC=~z1+Count1,
	 treat.var="treattime")
summary(fit2)

## -----------------------------------------------------------------------------
sse <- phreg_rct(Event(start,time,statusD)~A0.f+A1t+cluster(id),data=rr,
	 typesR=c("non","R0","R1","R01"),typesC=c("non","C","dynC"),
         treat.var="treattime",
	 RCT=FALSE,treat.model=At.f~z1*factor(Count2),
         augmentR0=~z1,augmentR1=~z1,augmentC=~z1+Count1+A1t)
summary(sse)

## -----------------------------------------------------------------------------
sessionInfo()

