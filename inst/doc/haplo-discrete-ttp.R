## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
types <- c("DCGCGCTCACG","DTCCGCTGACG","ITCAGTTGACG","ITCCGCTGAGG")

## some haplotypes frequencies for simulations 
data(haplo)
hapfreqs <- haplo$hapfreqs 
print(hapfreqs)

## -----------------------------------------------------------------------------
www <-which(hapfreqs$haplotype %in% types)
hapfreqs$freq[www]

baseline=hapfreqs$haplotype[9]
baseline

## -----------------------------------------------------------------------------
haploX  <- haplo$haploX
dlist(haploX,.~id|id %in% c(1,4,7))

## -----------------------------------------------------------------------------
ghaplos <- haplo$ghaplos
head(ghaplos)

## -----------------------------------------------------------------------------
designftypes <- function(x,sm=0) {
hap1=x[1]
hap2=x[2]
if (sm==0) y <- 1*( (hap1==types) | (hap2==types))
if (sm==1) y <- 1*(hap1==types) + 1*(hap2==types)
return(y)
}


## -----------------------------------------------------------------------------
haploX$time <- haploX$times
Xdes <- model.matrix(~factor(time),haploX)
colnames(Xdes) <- paste("X",1:ncol(Xdes),sep="")
X <- dkeep(haploX,~id+y+time)
X <- cbind(X,Xdes)
Haplos <- dkeep(ghaplos,~id+"haplo*"+p)
desnames=paste("X",1:6,sep="")   # six X's related to 6 cycles 
head(X)


## -----------------------------------------------------------------------------
out <- haplo.surv.discrete(X=X,y="y",time.name="time",
      Haplos=Haplos,desnames=desnames,designfunc=designftypes) 
names(out$coef) <- c(desnames,types)
out$coef
summary(out)

## -----------------------------------------------------------------------------
tcoef=c(-1.93110204,-0.47531630,-0.04118204,-1.57872602,-0.22176426,-0.13836416,
0.88830288,0.60756224,0.39802821,0.32706859)

cbind(out$coef,tcoef)

## -----------------------------------------------------------------------------
head(out$X,10)

## -----------------------------------------------------------------------------
sessionInfo()

