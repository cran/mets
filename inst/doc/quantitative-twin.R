## ----include=FALSE,echo=FALSE,message=FALSE,warning=FALSE---------------------
options(warn=-1, family="Times")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  #dev="png",
  #dpi=72,
  out.width = "70%")
library("mets")

## ----install, eval=FALSE------------------------------------------------------
#  # install.packages("remotes")
#  remotes::install_github("kkholst/mets", dependencies="Suggests")

## ----twinbmi------------------------------------------------------------------
library(mets)
data("twinbmi")
head(twinbmi)

## ----twinwide-----------------------------------------------------------------
twinwide <- fast.reshape(twinbmi, id="tvparnr",varying=c("bmi"))
head(twinwide)

## ----scatterdens, echo=FALSE,message=FALSE,warning=FALSE----------------------
library("cowplot")

scatterdens <- function(x) {
    require(ggplot2)
    sp <- ggplot(x,
                aes_string(colnames(x)[1], colnames(x)[2])) +
        theme_minimal() +
        geom_point(alpha=0.3) + geom_density_2d()
    xdens <- ggplot(x, aes_string(colnames(x)[1],fill=1)) +
        theme_minimal() +
        geom_density(alpha=.5)+
        theme(axis.text.x = element_blank(),
	      legend.position = "none") + labs(x=NULL)
    ydens <- ggplot(x, aes_string(colnames(x)[2],fill=1)) +
        theme_minimal() +
        geom_density(alpha=.5) +
        theme(axis.text.y = element_blank(),
	      axis.text.x = element_text(angle=90, vjust=0),
	      legend.position = "none") +
        labs(x=NULL) +
        coord_flip()
    g <- plot_grid(xdens,NULL,sp,ydens,
                  ncol=2,nrow=2,
                  rel_widths=c(4,1.4),rel_heights=c(1.4,4))
    return(g)
}

## ----scatter1, warning=FALSE,message=FALSE,fig.cap="Scatter plot of logarithmic BMI measurements in MZ twins"----
mz <- log(subset(twinwide, zyg=="MZ")[,c("bmi1","bmi2")])
scatterdens(mz)

## ----scatter2, warning=FALSE,message=FALSE,fig.cap="Scatter plot of logarithmic BMI measurements in DZ twins"----
dz <- log(subset(twinwide, zyg=="DZ")[,c("bmi1","bmi2")])
scatterdens(dz)

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
dm <- Expand(twinbmi,
	    bmi=0,
	    gender=c("male"),
	    age=seq(33,61,length.out=50))
df <- Expand(twinbmi,
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
l0 <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="sat")

# different marginals (but within pair)
lf <- twinlm(bmi ~ age+gender, data=dd,DZ="DZ", zyg="zyg", id="tvparnr", type="flex")

# same marginals but free correlation with MZ, DZ 
lu <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="u")
estimate(lu,contr(5:6,6))
estimate(lu)

lf <- twinlm(bmi ~ zyg, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="flex")
coef(lf)


###sink("lu-est-summary.txt")
lu <- twinlm(bmi ~ zyg, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="u")
summary(lu)
estimate(lu)
crossprod(iid(lu))^.5
###sink()

vcov(lu)

estimate(lu)
dim(iid(lu))

estimate(lu,contr(4:5,5))

estimate(coef=coef(lu),vcov=vcov(lu),contr(4:5,5))

wald.test(coef=coef(lu),vcov=vcov(lu),contrast=c(0,0,0,1,-1))

## -----------------------------------------------------------------------------

l <- twinlm(bmi ~ ns(age,1)+gender, data=twinbmi,
           DZ="DZ", zyg="zyg", id="tvparnr", type="cor", missing=TRUE)
summary(l)

## -----------------------------------------------------------------------------
estimate(l,contr(5:6,6))

## -----------------------------------------------------------------------------
l <- twinlm(bmi ~ ns(age,1)+gender, data=twinbmi,
           DZ="DZ", zyg="zyg", id="tvparnr", type="cor", missing=TRUE)
summary(l)

## -----------------------------------------------------------------------------
ace0 <- twinlm(bmi ~ age+gender, data=dd, DZ="DZ", zyg="zyg", id="tvparnr", type="ace")

