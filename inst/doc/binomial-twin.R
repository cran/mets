## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
library(mets)
data(twinstut)
twinstut$binstut <- 1*(twinstut$stutter=="yes")
twinsall <- twinstut
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
head(twinstut)
twinstut <- subset(twinstut,tvparnr < 3000)

## -----------------------------------------------------------------------------
margbin <- glm(binstut~factor(sex)+age,data=twinstut,family=binomial())
summary(margbin)

## ----twostage1----------------------------------------------------------------
bina <- binomial_twostage(margbin,data=twinstut,var.link=1,
                       clusters=twinstut$tvparnr,detail=0)
summary(bina)

## ----twostage2----------------------------------------------------------------
# design for OR dependence 
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
bin <- binomial_twostage(margbin,data=twinstut,var.link=1,
                          clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bin)

## ----twostage3----------------------------------------------------------------
twinstut$cage <- scale(twinstut$age)
theta.des <- model.matrix( ~-1+factor(zyg)+factor(sex),data=twinstut)
bina <- binomial_twostage(margbin,data=twinstut,var.link=1,
                          clusters=twinstut$tvparnr,theta.des=theta.des)
summary(bina)

## -----------------------------------------------------------------------------
library(mets)
data(twinstut)
twinstut <- subset(twinstut,zyg%in%c("mz","dz"))
twinstut$binstut <- 1*(twinstut$stutter=="yes")
head(twinstut)
twinstut <- subset(twinstut,tvparnr < 2000)

## ----biprobit1----------------------------------------------------------------
b1 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="un")
summary(b1)

## ----bptwin1------------------------------------------------------------------
b1 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="ace")
summary(b1)

## ----bptwin2------------------------------------------------------------------
b0 <- bptwin(binstut~sex,data=twinstut,id="tvparnr",zyg="zyg",DZ="dz",type="ae")
summary(b0)

## ----addgamma1----------------------------------------------------------------
theta.des <- model.matrix( ~-1+factor(zyg),data=twinstut)
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin <- binomial_twostage(margbin,data=twinstut,model="gamma",
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=1,
     theta.des=theta.des)
summary(bintwin)

# test for same dependence in MZ and DZ 
theta.des <- model.matrix( ~factor(zyg),data=twinstut)
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin <- binomial_twostage(margbin,data=twinstut,model="gamma",
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=1,
     theta.des=theta.des)
summary(bintwin)

## ----polygenic1---------------------------------------------------------------
out <- twin_polygen_design(twinstut,id="tvparnr",zygname="zyg",zyg="dz",type="ace")
head(cbind(out$des.rv,twinstut$tvparnr),10)
out$pardes

## ----polygenic2---------------------------------------------------------------
margbin <- glm(binstut~sex,data=twinstut,family=binomial())
bintwin1 <- binomial_twostage(margbin,data=twinstut,
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin1)

## -----------------------------------------------------------------------------
concordanceTwinACE(bintwin1,type="ace")

## ----polygenic_ae-------------------------------------------------------------
out <- twin.polygen.design(twinstut,id="tvparnr",zygname="zyg",zyg="dz",type="ae")

bintwin <- binomial_twostage(margbin,data=twinstut,
     clusters=twinstut$tvparnr,detail=0,theta=c(0.1)/1,var.link=0,
     random.design=out$des.rv,theta.des=out$pardes)
summary(bintwin)

## -----------------------------------------------------------------------------
concordanceTwinACE(bintwin,type="ae")

## -----------------------------------------------------------------------------
sessionInfo()

