## ----include = FALSE, label=setup---------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  ##dev="png",
  dpi=50,
  fig.width=7.15, fig.height=5.5,
  out.width="600px",
  fig.retina=1,
  comment = "#>"
  )
fullVignette <- Sys.getenv("_R_FULL_VIGNETTE_") %in% c("1","TRUE") || length(list.files("data"))==0
saveobj <- function(obj, null_elements) {
  x <- get(obj, envir=parent.frame())
  if (length(null_elements)>0) {
    print(object.size(x))
    x[null_elements] <- NULL
    attributes(x)[null_elements] <- NULL
    print(object.size(x))
  }
  saveRDS(x, paste0("data/", obj, ".rds"), compress="xz")
}

library(mets)
 cols <- c("darkred","darkblue","black")
 ltys <- c(1,3,2)
 fig_w <- 5
 fig_h <- 5
 savefig <- TRUE

## ----results="hide", echo=FALSE, eval=!fullVignette---------------------------
## To save time building the vignettes on CRAN, we cache time consuming computations
fitco1 <- readRDS("data/fitco1.rds")
fitco2 <- readRDS("data/fitco2.rds")
fitco3 <- readRDS("data/fitco3.rds")
fitco4 <- readRDS("data/fitco4.rds")
fitace <- readRDS("data/fitace.rds")
fitde <- readRDS("data/fitde.rds")
cse <- readRDS("data/cse.rds")
slr <- readRDS("data/slr.rds")
outacem <- readRDS("data/outacem.rds")
b0 <- readRDS("data/b0.rds")
b1 <- readRDS("data/b1.rds")
b2 <- readRDS("data/b2.rds")
a1 <- readRDS("data/a1.rds")
h2 <- readRDS("data/h2.rds")
concMZ <- readRDS("data/concMZ.rds")
s_mz_country <- readRDS("data/s_mz_country.rds")
s_dz_country <- readRDS("data/s_dz_country.rds")

## ----label=data-prt-----------------------------------------------------------
library(mets)
 set.seed(122)
 data(prt)
 
 dtable(prt,~status+cancer)
 dtable(prt,~zyg+country,level=1)

## ----label=survival-marginal--------------------------------------------------
 # Marginal Cox model here stratified on country without covariates 
 margph <- phreg(Surv(time,cancer)~strata(country)+cluster(id),data=prt)
 plot(margph)

## ----label=survival-pairwise1, eval=fullVignette------------------------------
# # Clayton-Oakes, MLE , overall variance
# fitco1<-twostageMLE(margph,data=prt,theta=2.7)

## -----------------------------------------------------------------------------
 summary(fitco1)

## ----label=survival-pairwise2, eval=fullVignette------------------------------
# fitco2 <- survival.twostage(margph,data=prt,theta=2.7,clusters=prt$id,var.link=0)

## -----------------------------------------------------------------------------
 summary(fitco2)

## ----label=survival-pairwise3, eval=fullVignette------------------------------
#  mm <- model.matrix(~-1+factor(zyg),prt)
#  fitco3<-twostageMLE(margph,data=prt,theta=1,theta.des=mm)

## -----------------------------------------------------------------------------
 summary(fitco3)

## ----label=survival-pairwise4, eval=fullVignette------------------------------
# fitco4 <- survival.twostage(margph,data=prt,theta=1,clusters=prt$id,var.link=0,theta.des=mm)

## -----------------------------------------------------------------------------
 summary(fitco4)
 round(estimate(coef=fitco4$coef,vcov=fitco4$var.theta)$coefmat[,c(1,3:4)],2)

 ## mz kendalls tau
 kendall.ClaytonOakes.twin.ace(fitco4$theta[2],0,K=1000)$mz.kendall
 ## dz kendalls tau
 kendall.ClaytonOakes.twin.ace(fitco4$theta[1],0,K=1000)$mz.kendall

## ----label=survival-polygenic1, eval=fullVignette-----------------------------
#  ### setting up design for random effects and parameters of random effects
#  desace <- twin.polygen.design(prt,type="ace")
# 
#  ### ace model
#  fitace <- survival.twostage(margph,data=prt,theta=1,
#        clusters=prt$id,var.link=0,model="clayton.oakes",
#        numDeriv=1,random.design=desace$des.rv,theta.des=desace$pardes)

## -----------------------------------------------------------------------------
 summary(fitace)

## ----label=survival-polygenic2, eval=fullVignette-----------------------------
#  ### ace model with positive random effects variances
#  # fitacee <- survival.twostage(margph,data=prt,theta=1,
#  #      clusters=prt$id,var.link=1,model="clayton.oakes",
#  #      numDeriv=1,random.design=desace$des.rv,theta.des=desace$pardes)
#  #summary(fitacee)
# 
#  ### ae model
#  #desae <- twin.polygen.design(prt,type="ae")
#  #fitae <- survival.twostage(margph,data=prt,theta=1,
#  #      clusters=prt$id,var.link=0,model="clayton.oakes",
#  #      numDeriv=1,random.design=desae$des.rv,theta.des=desae$pardes)
#  #summary(fitae)
# 
#  ### de model
#  desde <- twin.polygen.design(prt,type="de")
#  fitde <- survival.twostage(margph,data=prt,theta=1,clusters=prt$id,var.link=0,model="clayton.oakes",
# numDeriv=1,random.design=desde$des.rv,theta.des=desde$pardes)
# 

## -----------------------------------------------------------------------------
summary(fitde)

## -----------------------------------------------------------------------------
prt <-  force.same.cens(prt,cause="status")

 dtable(prt,~status+cancer)
 dtable(prt,~status+country)
 dtable(prt,~zyg+country)

## ----label=concordance--------------------------------------------------------
 ## cumulative incidence with cluster standard errors.
 cif1 <- cif(Event(time,status)~strata(country)+cluster(id),prt,cause=2)
 plot(cif1,se=1)

 cifa <- cif(Event(time,status)~+1,prt,cause=2)

 ### concordance estimator, ignoring country differences. 
 p11 <- bicomprisk(Event(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2))
p11mz <- p11$model$"MZ"
p11dz <- p11$model$"DZ" 

## -----------------------------------------------------------------------------
 par(mfrow=c(1,2))
 ## Concordance
 plot(p11mz,ylim=c(0,0.1));
 plot(p11dz,ylim=c(0,0.1));

## ----label=concordance2-------------------------------------------------------
 library(prodlim)
 outm <- prodlim(Hist(time,status)~+1,data=prt)

 cifzyg <- cif(Event(time,status)~+strata(zyg)+cluster(id),data=prt,cause=2)
 cifprt <- cif(Event(time,status)~country+cluster(id),data=prt,cause=2)
     
 times <- 70:100
 cifmz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="MZ")) ## cause is 2 (second cause) 
 cifdz <- predict(outm,cause=2,time=times,newdata=data.frame(zyg="DZ"))
    
 ### concordance for MZ and DZ twins<
 cc <- bicomprisk(Event(time,status)~strata(zyg)+id(id),data=prt,cause=c(2,2),prodlim=TRUE)
 ccdz <- cc$model$"DZ"
 ccmz <- cc$model$"MZ"
     
 cdz <- casewise(ccdz,outm,cause.marg=2) 
 cmz <- casewise(ccmz,outm,cause.marg=2)

 dd <- bicompriskData(Event(time,status)~country+strata(zyg)+id(id),data=prt,cause=c(2,2))
 conczyg <- cif(Event(time,status)~strata(zyg)+cluster(id),data=dd,cause=1)

 par(mfrow=c(1,2))
 plot(conczyg,se=TRUE,col=cols[2:1], lty=ltys[2:1], legend=FALSE,xlab="Age",ylab="Concordance")
 legend("topleft",c("concordance-MZ","concordance-DZ"),col=cols[1:2],lty=ltys[1:2])

 plot(cmz,ci=NULL,ylim=c(0,.8),xlim=c(70,97),legend=FALSE,col=cols[c(1,3,3)],lty=ltys[c(1,3,3)],
      ylab="Casewise",xlab="Age")
  plot(cdz,ci=NULL,ylim=c(0,.8),xlim=c(70,97),legend=FALSE,ylab="Casewise",xlab="Age",
      col=c(cols[2],NA,NA), lty=ltys[c(2,3,3)], add=TRUE)
 with(data.frame(cmz$casewise),plotConfRegionSE(time,casewise.conc,se.casewise,col=cols[1]))
 with(data.frame(cdz$casewise),plotConfRegionSE(time,casewise.conc,se.casewise,col=cols[2]))
 legend("topleft",c("casewise-MZ","casewise-DZ","marginal"),col=cols, lty=ltys, bg="white")

 summary(cdz)
 summary(cmz)

 cpred(cmz$casewise,c(70,80))
 cpred(cdz$casewise,c(70,80))

## ----label=concordance3-------------------------------------------------------

 dd <- bicompriskData(Event(time,status)~country+strata(zyg)+id(id),data=prt,cause=c(2,2))
 conczyg <- cif(Event(time,status)~strata(zyg)+cluster(id),data=dd,cause=1)

 par(mfrow=c(1,2))
 plot(conczyg,se=TRUE,legend=FALSE,xlab="Age",ylab="Concordance")
 legend("topleft",c("concordance-DZ","concordance-MZ"),col=c(1,2),lty=1)
 plot(cmz,ci=NULL,ylim=c(0,0.6),xlim=c(70,100),legend=FALSE,col=c(2,3,3),ylab="Casewise",xlab="Age",lty=c(1,3))
 plot(cdz,ci=NULL,ylim=c(0,0.6),xlim=c(70,100),legend=FALSE,ylab="Casewise",xlab="Age",
      col=c(1,3,3), add=TRUE, lty=c(2,3))
 legend("topleft",c("casewise-MZ","casewise-DZ","marginal"),col=c(2,1,3),lty=1)
 with(data.frame(cmz$casewise),plotConfRegionSE(time,casewise.conc,se.casewise,col=2))
 with(data.frame(cdz$casewise),plotConfRegionSE(time,casewise.conc,se.casewise,col=1))


## ----label=concordance4, eval=fullVignette------------------------------------
#  ### new version of Casewise for specific time-point based on binreg
#  dd <- bicompriskData(Event(time,status)~country+strata(zyg)+id(id),data=prt,cause=c(2,2))
#  newdata <- data.frame(zyg=c("DZ","MZ"),id=1)
# 
#  ## concordance
#  bcif1 <- binreg(Event(time,status)~-1+factor(zyg)+cluster(id),dd,time=80,cause=1,cens.model=~strata(zyg))
#  pconc <- predict(bcif1,newdata)
# 
#  ## marginal estimates
#  mbcif1 <- binreg(Event(time,status)~cluster(id),prt,time=80,cause=2)
#  mc <- predict(mbcif1,newdata)
# 
#  ### casewise with improved se's from log-scale
#  cse <- binregCasewise(bcif1,mbcif1)

## -----------------------------------------------------------------------------
 cse 

## ----label=semiparconc, eval=fullVignette-------------------------------------
#  ### semi-parametric modelling of concordance
#  dd <- bicompriskData(Event(time,status)~country+strata(zyg)+id(id),data=prt,cause=c(2,2))
#  regconc <- cifreg(Event(time,status)~country*zyg,data=dd,prop=NULL)
#  regconc
#  ### interaction test
#  wald.test(regconc,coef.null=5:7)
# 
#  regconc <- cifreg(Event(time,status)~country+zyg,data=dd,prop=NULL)
#  regconc
# 
#  ## logistic link
#  logitregconc <- cifreg(Event(time,status)~country+zyg,data=dd)
#  slr <- summary(logitregconc)

## -----------------------------------------------------------------------------
slr
### library(Publish)
### publish(round(slr$exp.coef[,-c(2,5)],2),latex=TRUE,digits=2)

## ----label=additive_gamma, eval=fullVignette----------------------------------
#  timereg <- 0
# if (timereg==1) {
#   times <- seq(50,90,length.out=5)
#   cif1 <- timereg::comp.risk(Event(time,status)~-1+factor(country)+cluster(id),prt,
# 		   cause=2,times=times,max.clust=NULL)
# 
#   mm <- model.matrix(~-1+factor(zyg),prt)
#   out1<-random.cif(cif1,data=prt,cause1=2,cause2=2,theta=1,
# 		  theta.des=mm,same.cens=TRUE,step=0.5)
#   summary(out1)
#   round(estimate(coef=out1$theta,vcov=out1$var.theta)$coefmat[,c(1,3:4)],2)
# 
#   desace <- twin.polygen.design(prt,type="ace")
# 
#   outacem <- Grandom.cif(cif1,data=prt,cause1=2,cause2=2,
#   	 same.cens=TRUE,theta=c(0.45,0.15),var.link=0,
#          step=0.5,theta.des=desace$pardes,random.design=desace$des.rv)
#   ##outacem$score
# }

## -----------------------------------------------------------------------------
timereg <- 0
if (timereg==1) {
  summary(outacem)

 ###  variances
 estimate(coef=outacem$theta,vcov=outacem$var.theta,f=function(p) p/sum(p)^2)

 ## AE polygenic model
 # desae <- twin.polygen.design(prt,type="ae")
 # outaem <- Grandom.cif(cif1,data=prt,cause1=2,cause2=2,
 #    same.cens=TRUE,theta=c(0.45,0.15),var.link=0,
 #        step=0.5,theta.des=desae$pardes,random.design=desae$des.rv)
 # outaem$score
 # summary(outaem)
 # estimate(coef=outaem$theta,vcov=outaem$var.theta,f=function(p)     p/sum(p)^2)

 ## AE polygenic model
 # desde <- twin.polygen.design(prt,type="de")
 # outaem <- Grandom.cif(cif1,data=prt,cause1=2,cause2=2,
 #   same.cens=TRUE,theta=c(0.35),var.link=0,
 #   step=0.5,theta.des=desde$pardes,random.design=desde$des.rv)
 # outaem$score
 # summary(outaem)
 # estimate(coef=outaem$theta,vcov=outaem$var.theta,f=function(p) p/sum(p)^2)

  times <- 90
  cif1 <- timereg::comp.risk(Event(time,status)~-1+factor(country)+cluster(id),prt,
		   cause=2,times=times,max.clust=NULL)

  mm <- model.matrix(~-1+factor(zyg),prt)
  out1<-random.cif(cif1,data=prt,cause1=2,cause2=2,theta=1,
		  theta.des=mm,same.cens=TRUE,step=0.5)
  summary(out1)
  round(estimate(coef=out1$theta,vcov=out1$var.theta)$coefmat[,c(1,3:4)],2)

 desde <- twin.polygen.design(prt,type="de")
 outaem <- Grandom.cif(cif1,data=prt,cause1=2,cause2=2,
	same.cens=TRUE,theta=c(0.35),var.link=0,
        step=0.5,theta.des=desde$pardes,random.design=desde$des.rv)
 outaem$score
 summary(outaem)
 estimate(coef=outaem$theta,vcov=outaem$var.theta,f=function(p) p/sum(p)^2)
}

## ----label=probit1------------------------------------------------------------
rm(prt)
data(prt)
prt0 <-  force.same.cens(prt, cause="status", cens.code=0, time="time", id="id")
prt0$country <- relevel(prt0$country, ref="Sweden")
prt_wide <- fast.reshape(prt0, id="id", num="num", varying=c("time","status","cancer"))
prt_time <- subset(prt_wide,  cancer1 & cancer2, select=c(time1, time2, zyg))
tau <- 95
tt <- seq(70, tau, length.out=5) ## Time points to evaluate model in

## ----b0, eval=fullVignette----------------------------------------------------
# b0 <- bptwin.time(cancer ~ 1, data=prt0, id="id", zyg="zyg", DZ="DZ", type="cor",
#               cens.formula=Surv(time,status==0)~zyg, breaks=tau)

## -----------------------------------------------------------------------------
summary(b0)

## ----label=liability_ace1, eval=fullVignette----------------------------------
# b1 <- bptwin.time(cancer ~ 1, data=prt0, id="id", zyg="zyg", DZ="DZ", type="ace",
#               cens.formula=Surv(time,status==0)~zyg, breaks=tau)

## -----------------------------------------------------------------------------
summary(b1)

## -----------------------------------------------------------------------------
AIC(b0, b1)

## ----label=liability_ace_country, eval=fullVignette---------------------------
# b2 <- bptwin.time(cancer ~ country, data=prt0, id="id", zyg="zyg", DZ="DZ", type="ace",
#               cens.formula=Surv(time,status==0)~zyg+country, breaks=95)

## -----------------------------------------------------------------------------
summary(b2)

## ----label=bptime1, eval=fullVignette-----------------------------------------
# bt0 <- bptwin.time(cancer ~ 1, data=prt0, id="id", zyg="zyg", DZ="DZ", type="ace",
#               cens.formula=Surv(time,status==0)~zyg,
#               summary.function=function(x) x, breaks=tt)
# h2 <- Reduce(rbind, lapply(bt0$coef, function(x) x$heritability))[,c(1,3,4),drop=FALSE]
# concMZ <- Reduce(rbind, lapply(bt0$coef, function(x) x$probMZ["Concordance",,drop=TRUE]))
# 

## -----------------------------------------------------------------------------
par(mfrow=c(1,2))
plot(tt, h2[,1], type="s", lty=1, col=cols[3], xlab="Age", ylab="Heritability", ylim=c(0,1))
lava::confband(tt, h2[,2], h2[,3],polygon=TRUE, step=TRUE, col=lava::Col(cols[3], 0.1), border=NA)
plot(tt, concMZ[,1], type="s", lty=1, col=cols[1], xlab="Age", ylab="Concordance", ylim=c(0,.1))
lava::confband(tt, concMZ[,2], concMZ[,3],polygon=TRUE, step=TRUE, col=lava::Col(cols[1], 0.1), border=NA)

## ----label=biprobittime1------------------------------------------------------
system.time(a.mz <- biprobit.time(cancer~1, id="id", data=subset(prt0, zyg=="MZ"),
                               cens.formula = Surv(time,status==0)~1, pairs.only=TRUE,
                                breaks=tt))
system.time(a.dz <- biprobit.time(cancer~1, id="id", data=subset(prt0, zyg=="DZ"),
                               cens.formula = Event(time,status==0)~1, pairs.only=TRUE,
                               breaks=tt))

#system.time(a.zyg <- biprobit.time(cancer~1, rho=~1+zyg, id="id", data=prt, 
#                               cens.formula = Event(time,status==0)~1,
#                               eqmarg=FALSE, fix.cens.weight
#                               breaks=seq(75,100,by=10)))

a.mz
a.dz

plot(conczyg,se=TRUE,legend=FALSE,xlab="Age",ylab="Concordance", ylim=c(0,0.07))
plot(a.mz, ylim=c(0,.07), col=cols[1], lty=ltys[1], legend=FALSE, add=TRUE)
plot(a.dz, col=cols[2], lty=ltys[2], add=TRUE)

## ----label=biprobittime2, eval=fullVignette-----------------------------------
# a.mz_country <- biprobit.time(cancer~country, id="id", data=subset(prt0, zyg=="MZ"),
#                                cens.formula = Surv(time,status==0)~country, pairs.only=TRUE,
#                                 breaks=tt)
# system.time(a.dz_country <- biprobit.time(cancer~country, id="id", data=subset(prt0, zyg=="DZ"),
#                                cens.formula = Event(time,status==0)~country, pairs.only=TRUE,
#                                breaks=tt))
# 
# s_mz_country <- summary(a.mz_country)
# s_dz_country <- summary(a.dz_country)

## -----------------------------------------------------------------------------
s_mz_country
s_dz_country

## ----label=liability_ace_time1, eval=fullVignette-----------------------------
# ## ACE model (time-varying) with and without adjustment for country
# a1 <- bptwin.time(cancer~1, id="id", data=prt0, type="ace",
#                               zyg="zyg", DZ="DZ",
#                               cens.formula=Surv(time,status==0)~zyg,
#                               breaks=tt)
# 
# #a2 <- bptwin.time(cancer~country, id="id", data=prt0, #type="ace",
# #                              zyg="zyg", DZ="DZ",
# #                              #cens.formula=Surv(time,status==0)~country+zyg,
# #                              breaks=tt)

## -----------------------------------------------------------------------------
plot(a.mz, which=c(6), xlab="Age", ylab="Correlation", ylim=c(0,1), col=cols[1], lty=ltys[1], legend=NULL, alpha=.1)
plot(a.dz, which=c(6), col=cols[2], lty=ltys[2], legend=NULL, add=TRUE, alpha=.1)
legend("topleft", c("MZ tetrachoric correlation", "DZ tetrachoric correlation"),
       col=cols, lty=ltys, lwd=2)

plot(a.mz, which=c(4), xlab="Age", ylab="Relative Recurrence Risk",
     ylim=c(1,20), col=cols[1], lty=ltys[1], legend=NULL, lwd=2, alpha=.1)
plot(a.dz, which=c(4), col=cols[2], lty=ltys[2], legend=NULL, add=TRUE, lwd=2, alpha=.1)
legend("topright", c("MZ relative recurrence risk", "DZ relative recurrence risk"),
       col=cols, lty=ltys, lwd=2)

plot(a1, which=c(5,6), xlab="Age", ylab="Correlation", ylim=c(0,1), col=cols[1:2], lty=ltys[1:2], lwd=2, alpha=0.1,
     legend=c("MZ tetrachoric correlation", "DZ tetrachoric correlation"))

plot(a1, which=c(1), xlab="Age", ylim=c(0,1), col="black", lty=1, ylab="Heritability", legend=NULL, alpha=.1)

## -----------------------------------------------------------------------------
sessionInfo()

## ----saveobj, results="hide", echo=FALSE, eval=fullVignette-------------------
# ## To save time building the vignettes on CRAN, we cache time consuming computations
# 
# rms <- c('id','theta.iid','theta.des','marginal.trunc',
#   'loglikeiid','marginal.surv','theta.iid.naive',
#   'antclust','secluster','cluster.call','trunclikeiid',
#   'logl.iid','score.iid','clusters')
# tmp <- lapply(as.list(paste0("fitco", 1:4)),
#               function(x) saveobj(x, rms))
# 
# rms <- c('random.design','marginal.surv','marginal.trunc','theta.iid','score.iid','loglikeiid','trunclikeiid','antclust','cluster.call','secluster','clusters')
# saveobj("fitace", rms)
# saveobj("fitde", rms)
# 
# saveobj("cse", NULL)
# saveobj("slr", NULL)
# 
# rms <- c("theta.iid", "Clusters", "p11")
# saveobj("outacem", rms)
# 
# rms <- c("model.frame", "score", "id", "logLik")
# saveobj("b0", rms)
# saveobj("b1", rms)
# saveobj("b2", rms)
# saveobj("a1", rms)
# saveobj("h2", NULL)
# saveobj("concMZ", NULL)
# saveobj("s_mz_country", NULL)
# saveobj("s_dz_country", NULL)

