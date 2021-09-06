## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets)

data(ttpd) 
dtable(ttpd,~entry+time2)
out <- interval.logitsurv.discrete(Interval(entry,time2)~X1+X2+X3+X4,ttpd)
summary(out)

## -----------------------------------------------------------------------------
set.seed(1000) # to control output in simulatins for p-values below.
n <- 10000
Z <- matrix(rbinom(n*4,1,0.5),n,4)
outsim <- simlogitSurvd(out$coef,Z)
outsim <- transform(outsim,left=time,right=time+1)
outsim <- dtransform(outsim,right=Inf,status==0)
outss <- interval.logitsurv.discrete(Interval(left,right)~+X1+X2+X3+X4,outsim)
summary(outss)

pred <- predictlogitSurvd(out,se=TRUE)
plotSurvd(pred,se=TRUE)

## -----------------------------------------------------------------------------
test <- 0 
if (test==1) {

require(icenReg)
data(IR_diabetes)
IRdia <- IR_diabetes
## removing fully observed data in continuous version, here making it a discrete observation 
IRdia <- dtransform(IRdia,left=left-1,left==right)
dtable(IRdia,~left+right,level=1)

ints <- with(IRdia,dInterval(left,right,cuts=c(0,5,10,20,30,40,Inf),show=TRUE) )
}

## -----------------------------------------------------------------------------
if (test==1) {
ints$Ileft <- ints$left
ints$Iright <- ints$right
IRdia <- cbind(IRdia,data.frame(Ileft=ints$Ileft,Iright=ints$Iright))
dtable(IRdia,~Ileft+Iright)
# 
#       Iright   1   2   3   4   5 Inf
# Ileft                               
# 0             10   1  34  25   4   0
# 1              0  55  19  17   1   1
# 2              0   0 393  16   4   0
# 3              0   0   0 127   1   0
# 4              0   0   0   0  21   0
# 5              0   0   0   0   0   2

outss <- interval.logitsurv.discrete(Interval(Ileft,Iright)~+gender,IRdia)
#            Estimate Std.Err    2.5%    97.5%   P-value
# time1        -3.934  0.3316 -4.5842 -3.28418 1.846e-32
# time2        -2.042  0.1693 -2.3742 -1.71038 1.710e-33
# time3         1.443  0.1481  1.1530  1.73340 1.911e-22
# time4         3.545  0.2629  3.0295  4.06008 1.976e-41
# time5         6.067  0.7757  4.5470  7.58784 5.217e-15
# gendermale   -0.385  0.1691 -0.7165 -0.05351 2.283e-02
summary(outss)
outss$ploglik
# [1] -646.1946

fit <- ic_sp(cbind(Ileft, Iright) ~ gender, data = IRdia, model = "po")
# 
# Model:  Proportional Odds
# Dependency structure assumed: Independence
# Baseline:  semi-parametric 
# Call: ic_sp(formula = cbind(Ileft, Iright) ~ gender, data = IRdia, 
#     model = "po")
# 
#            Estimate Exp(Est)
# gendermale    0.385     1.47
# 
# final llk =  -646.1946 
# Iterations =  6 
# Bootstrap Samples =  0 
# WARNING: only  0  bootstrap samples used for standard errors. 
# Suggest using more bootstrap samples for inference
summary(fit)

## sometimes NR-algorithm needs modifications of stepsize to run 
## outss <- interval.logitsurv.discrete(Interval(Ileft,Iright)~+gender,IRdia,control=list(trace=TRUE,stepsize=1.0))
}


## -----------------------------------------------------------------------------
sessionInfo()

