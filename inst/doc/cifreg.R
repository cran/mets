## ----include = FALSE----------------------------------------------------------
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
 n <- 400
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
gofFG(Event(time,status)~Z1+Z2,data=dats,cause=1)

## -----------------------------------------------------------------------------
### predictions with CI based on iid decomposition of baseline and beta
fg <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,cox.prep=TRUE)
Biid <- IIDbaseline.cifreg(fg,time=5)
pfgse <- FGprediid(Biid,dd)
pfgse
plot(pfg,ylim=c(0,0.2))
for (i in 1:4) lines(c(5,5)+i/10,pfgse[i,3:4],col=i,lwd=2)

## -----------------------------------------------------------------------------
run <- 0
if (run==1) {
 library(cmprsk)
 mm <- model.matrix(~Z1+Z2,dats)[,-1]
 cr <- with(dats,crr(time,status,mm))
 cbind(cr$coef,diag(cr$var)^.5,fg$coef,fg$se.coef,cr$coef-fg$coef,diag(cr$var)^.5-fg$se.coef)
#          [,1]      [,2]       [,3]      [,4]          [,5]          [,6]
# Z1  0.6968603 0.3876029  0.6968603 0.3876029 -1.534155e-09 -7.811395e-10
# Z2 -0.8592892 0.6245258 -0.8592892 0.6245258  8.537615e-14  4.430734e-11
}

if (run==1) {
 fg0 <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL,beta=cr$coef,no.opt=TRUE)
 cbind(diag(cr$var)^.5,fg0$se.coef,diag(cr$var)^.5-fg0$se.coef)
#           [,1]      [,2]          [,3]
# [1,] 0.3876029 0.3876029 -8.881784e-16
# [2,] 0.6245258 0.6245258  3.330669e-16
}

## -----------------------------------------------------------------------------
if (run==1) {
 dats$id <- 1:nrow(dats)
 dats$event <- factor(dats$status,0:2, labels=c("censor", "death", "other"))
 fgdats <- finegray(Surv(time,event)~.,data=dats)
 coxfg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ Z1+Z2 + cluster(id), weight=fgwt, data=fgdats)

 fg0 <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1,propodds=NULL)
 cbind( coxfg$coef,fg0$coef, coxfg$coef-fg0$coef)
#          [,1]       [,2]          [,3]
# Z1  0.6968603  0.6968603 -1.534153e-09
# Z2 -0.8592892 -0.8592892  8.726353e-14
 cbind(diag(coxfg$var)^.5,fg0$se.coef,diag(coxfg$var)^.5-fg0$se.coef)
#           [,1]      [,2]          [,3]
# [1,] 0.3889129 0.3876029  0.0013099907
# [2,] 0.6241225 0.6245258 -0.0004033148
 cbind(diag(coxfg$var)^.5,fg0$se1.coef,diag(coxfg$var)^.5-fg0$se1.coef)
#           [,1]      [,2]          [,3]
# [1,] 0.3889129 0.3889129 -7.839487e-10
# [2,] 0.6241225 0.6241225  4.448153e-11
}


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
if (run==1) {
 library(cmprsk)
 mm <- model.matrix(~Z1+Z2,datsnc)[,-1]
 cr <- with(datsnc,crr(time,status,mm))
 cbind(cifnc$coef-cr$coef, diag(cr$var)^.5-cifnc$se.coef)
#            [,1]          [,2]
# Z1 2.053291e-09 -1.048037e-09
# Z2 1.926237e-13  7.571876e-11
}

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
 n <- 400
 beta=c(0.0,-0.1,-0.5,0.3)
 dats <- simul.cifs(n,rho1,rho2,beta,rc=0.5,rate=7,type="logistic")
 dtable(dats,~status)
 dsort(dats) <- ~time

## -----------------------------------------------------------------------------
 or <- cifreg(Event(time,status)~Z1+Z2,data=dats,cause=1)
 summary(or)

## -----------------------------------------------------------------------------
sessionInfo()

