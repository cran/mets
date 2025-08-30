## ----include=FALSE,echo=FALSE,message=FALSE,warning=FALSE---------------------
options(warn=-1, family="Times")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=50,
  fig.width=7.15, fig.height=5.5,
  out.width="600px",
  fig.retina=1
  ##dev="png",
  ##dpi=72,
  ## out.width = "70%")
)
library("mets")

## ----install, eval=FALSE, echo=FALSE------------------------------------------
# # install.packages("remotes")
# remotes::install_github("kkholst/mets", dependencies="Suggests")

## ----twinbmi------------------------------------------------------------------
library(mets)
data("twinbmi")
head(twinbmi)

## ----twinwide-----------------------------------------------------------------
twinwide <- fast.reshape(twinbmi, id="tvparnr",varying=c("bmi"))
head(twinwide)

## ----scatter1, warning=FALSE,message=FALSE,fig.cap="Scatter plot of logarithmic BMI measurements in MZ twins"----
mz <- log(subset(twinwide, zyg=="MZ")[,c("bmi1","bmi2")])
plot_twin(mz)

## ----scatter2, warning=FALSE,message=FALSE,fig.cap="Scatter plot of logarithmic BMI measurements in DZ twins"----
dz <- log(subset(twinwide, zyg=="DZ")[,c("bmi1","bmi2")])
plot_twin(dz)

## -----------------------------------------------------------------------------
cor.test(mz[,1],mz[,2], method="spearman")

## -----------------------------------------------------------------------------
cor.test(dz[,1],dz[,2], method="spearman")

## ----gee----------------------------------------------------------------------
l0 <- lm(bmi ~ gender + I(age-40), data=twinbmi)
estimate(l0, id=twinbmi$tvparnr)

## -----------------------------------------------------------------------------
library("splines")
l1 <- lm(bmi ~ gender*ns(age,3), data=twinbmi)
marg1 <- estimate(l1, id=twinbmi$tvparnr)

## ----marg1, warning=FALSE,message=FALSE,fig.cap="Marginal association between BMI and Age for males and females."----
dm <- lava::Expand(twinbmi,
	    bmi=0,
	    gender=c("male"),
	    age=seq(33,61,length.out=50))
df <- lava::Expand(twinbmi,
	    bmi=0,
	    gender=c("female"),
	    age=seq(33,61,length.out=50))

plot(marg1, function(p) model.matrix(l1,data=dm)%*%p,
     data=dm["age"], ylab="BMI", xlab="Age",
     ylim=c(22,26.5))
plot(marg1, function(p) model.matrix(l1,data=df)%*%p,
     data=df["age"], col="red", add=TRUE)
legend("bottomright", c("Male","Female"),
       col=c("black","red"), lty=1, bty="n")

## -----------------------------------------------------------------------------
dd <- na.omit(twinbmi)

## ----lmsat, eval=FALSE--------------------------------------------------------
# l0 <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="sat")

## ----lmflex, eval=FALSE-------------------------------------------------------
# lf <- twinlm(bmi ~ age+gender, data=dd,DZ="DZ", zyg="zyg", id="tvparnr", type="flex")

## ----lmeqmarg-----------------------------------------------------------------
lu <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="eqmarg")
estimate(lu)

## -----------------------------------------------------------------------------
estimate(lu,lava::contr(5:6,6))

## ----ace----------------------------------------------------------------------
ace0 <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="ace")
summary(ace0)

