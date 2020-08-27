## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 options(warn=-1)
 set.seed(1000) # to control output in simulatins for p-values below.

 rho1 <- 0.2; rho2 <- 10
 n <- 4000
 beta=c(0.0,-0.1,-0.5,0.3)
 ## beta1=c(0.0,-0.1); beta2=c(-0.5,0.3)
 dats <- simul.cifs(n,rho1,rho2,beta,rc=0.5,rate=7)
 dtable(dats,~status)
 dsort(dats) <- ~time

## -----------------------------------------------------------------------------
 par(mfrow=c(1,2))
 cifs1 <- cif(Event(time,status)~strata(Z1,Z2),dats,cause=1)
 plot(cifs1)

 cifs2 <- cif(Event(time,status)~strata(Z1,Z2),dats,cause=2)
 plot(cifs2)

## -----------------------------------------------------------------------------
 fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
 summary(fg)

 dd <- expand.grid(Z1=c(-1,1),Z2=0:1)
 pfg <- predict(fg,dd)
 plot(pfg,ylim=c(0,0.2))

## -----------------------------------------------------------------------------
### library(cmprsk)
### mm <- model.matrix(~Z1+Z2,dats)[,-1]
### cr <- with(dats,crr(time,status,mm))
### 
### cbind(cr$coef,diag(cr$var)^.5,fg$coef,fg$se.coef,cr$coef-fg$coef,diag(cr$var)^.5-fg$se.coef)
#           [,1]       [,2]        [,3]       [,4]         [,5]          [,6]
# Z1  0.06882201 0.08184642  0.06882199 0.08184642 1.398233e-08 -9.805955e-11
# Z2 -0.36807861 0.16606441 -0.36807951 0.16606444 9.045060e-07 -2.529264e-08

### fg0 <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,beta=cr$coef,no.opt=TRUE)
### cbind(diag(cr$var)^.5,fg0$se.coef,diag(cr$var)^.5-fg0$se.coef)
#            [,1]       [,2]          [,3]
# [1,] 0.08184642 0.08184642 -8.975096e-11
# [2,] 0.16606441 0.16606441 -1.091264e-10

## -----------------------------------------------------------------------------
datsnc <- dtransform(dats,status=2,status==0)
dtable(datsnc,~status)
datsnc$id <- 1:n
datsnc$entry <- 0
max <- max(dats$time)+1
## for cause 2 add risk interaval 
datsnc2 <- subset(datsnc,status==2)
datsnc2 <- transform(datsnc2,entry=time)
datsnc2 <- transform(datsnc2,time=max)
datsncf <- rbind(datsnc,datsnc2)
#
cifnc <- cifreg(Event(time,status)~Z1+Z2,data=datsnc,cause=1,propodds=NULL)
cc <- coxph(Surv(entry,time,status==1)~Z1+Z2+cluster(id),datsncf,robust=TRUE)
cbind(cifnc$coef,cifnc$se.coef, cc$coef, diag(cc$var)^.5)
cbind(cifnc$coef-cc$coef, diag(cc$var)^.5-cifnc$se.coef)

## -----------------------------------------------------------------------------
### library(cmprsk)
### mm <- model.matrix(~Z1+Z2,datsnc)[,-1]
### cr <- with(datsnc,crr(time,status,mm))
### cbind(cifnc$coef-cr$coef, diag(cr$var)^.5-cifnc$se.coef)
#             [,1]          [,2]
# Z1  1.194188e-10 -6.742754e-11
# Z2 -8.828051e-07 -2.488335e-08

## -----------------------------------------------------------------------------
 fgcm <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,
		cens.model=~strata(Z1,Z2))
 summary(fgcm)
 summary(fg)

## -----------------------------------------------------------------------------
  fgaugS <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fg$E)
  summary(fgaugS)

  fgaugS2 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS$E)
  summary(fgaugS2)

  fgaugS3 <- FG_AugmentCifstrata(Event(time,status)~Z1+Z2+strata(Z1,Z2),data=dats,cause=1,E=fgaugS2$E)
  summary(fgaugS3)

## -----------------------------------------------------------------------------
 rho1 <- 0.2; rho2 <- 10
 n <- 4000
 beta=c(0.0,-0.1,-0.5,0.3)
 dats <- simul.cifs(n,rho1,rho2,beta,rc=0.5,rate=7,type="logistic")
 dtable(dats,~status)
 dsort(dats) <- ~time

## -----------------------------------------------------------------------------
 or <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1)
 summary(or)

## -----------------------------------------------------------------------------
sessionInfo()

