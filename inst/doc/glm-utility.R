## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
set.seed(100)

library(mets)
data(bmt); 
bmt$id <- sample(1:100,408,replace=TRUE)

glm1 <- glm(tcell~platelet+age,bmt,family=binomial)
summaryGLM(glm1)

## GEE robust standard errors
summaryGLM(glm1,id=bmt$id)

## -----------------------------------------------------------------------------
age <- seq(-2,2,by=0.1)
nd <- data.frame(platelet=0,age=seq(-2,2,by=0.1))
pnd <- predictGLM(glm1,nd)
head(pnd$pred)
plot(age,pnd$pred[,1],type="l",ylab="predictions",xlab="age",ylim=c(0,0.3))
plotConfRegion(age,pnd$pred[,2:3],col=2)
###matlines(age,pnd$pred[,-1],col=2,lty=1)

## -----------------------------------------------------------------------------
sessionInfo()

