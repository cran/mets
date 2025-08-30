## ----include = FALSE----------------------------------------------------------
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

dfactor(ttpd) <- entry.f~entry
out <- cumoddsreg(entry.f~X1+X2+X3+X4,ttpd)
summary(out)

## -----------------------------------------------------------------------------
set.seed(1000) # to control output in simulatins for p-values below.
n <- 200
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

data(ttpd) 
dtable(ttpd,~entry+time2)
ttpd <- dfactor(ttpd,fentry~entry)
out <- cumoddsreg(fentry~X1+X2+X3+X4,ttpd)
summary(out)

out$ploglik

if (test==1) {
### library(ordinal)
### out1 <- clm(fentry~X1+X2+X3+X4,data=ttpd)
### summary(out1)

# formula: fentry ~ X1 + X2 + X3 + X4
# data:    ttpd
# 
#  link  threshold nobs logLik   AIC     niter max.grad cond.H 
#  logit flexible  1000 -1676.46 3372.91 6(2)  1.17e-12 5.3e+02
# 
# Coefficients:
#    Estimate Std. Error z value Pr(>|z|)    
# X1  -0.9913     0.1171  -8.465  < 2e-16 ***
# X2  -0.6962     0.1156  -6.021 1.74e-09 ***
# X3  -0.3466     0.1150  -3.013  0.00259 ** 
# X4  -0.3223     0.1147  -2.810  0.00495 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Threshold coefficients:
#     Estimate Std. Error z value
# 0|1  -2.0064     0.1461 -13.733
# 1|2  -1.3940     0.1396  -9.984
# 2|3  -0.7324     0.1347  -5.435
# 3|4  -0.6266     0.1343  -4.667
# 4|5  -0.1814     0.1333  -1.361
# 5|6   0.2123     0.1342   1.582
}


## -----------------------------------------------------------------------------
sessionInfo()

