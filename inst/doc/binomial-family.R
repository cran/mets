## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 library(timereg)
 set.seed(100)
 data <- simbinClaytonOakes.family.ace(500,2,1,beta=NULL,alpha=NULL)
 data$number <- c(1,2,3,4)
 data$child <- 1*(data$number==3)
 head(data)

## -----------------------------------------------------------------------------
 aa <- margbin <- glm(ybin~x,data=data,family=binomial())
 summary(aa)

## -----------------------------------------------------------------------------
# make ace random effects design
out <- ace.family.design(data,member="type",id="cluster")
out$pardes
head(out$des.rv,4)

## -----------------------------------------------------------------------------
# fitting ace model for family structure
ts <- binomial.twostage(margbin,data=data,clusters=data$cluster,
theta=c(2,1),random.design=out$des.rv,theta.des=out$pardes)
summary(ts)
# true variance parameters
c(2,1)
# total variance 
3

## -----------------------------------------------------------------------------
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
pairs <- mm$pairs
dim(pairs)
head(pairs,12)

## -----------------------------------------------------------------------------
tsp <- binomial.twostage(margbin,data=data,clusters=data$cluster,theta=c(2,1),detail=0,
        random.design=out$des.rv,theta.des=out$pardes,pairs=pairs)
summary(tsp)

## -----------------------------------------------------------------------------
set.seed(100)
ssid <- sort(sample(1:nrow(pairs),nrow(pairs)/2))
tsd <- binomial.twostage(aa,data=data,clusters=data$cluster,
               theta=c(2,1),step=1.0,
               random.design=out$des.rv,iid=1,Nit=10,
  	           theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)

## -----------------------------------------------------------------------------
head(pairs[ssid,])
ids <- sort(unique(c(pairs[ssid,])))

pairsids <- c(pairs[ssid,])
pair.new <- matrix(fast.approx(ids,c(pairs[ssid,])),ncol=2)
head(pair.new)

dataid <- dsort(data[ids,],"cluster")
outid <- ace.family.design(dataid,member="type",id="cluster")
outid$pardes
head(outid$des.rv)

## -----------------------------------------------------------------------------
aa <- glm(ybin~x,data=dataid,family=binomial())
tsdid <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
         theta=c(2,1),random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)

## -----------------------------------------------------------------------------
pair.types <-  matrix(dataid[c(t(pair.new)),"type"],byrow=T,ncol=2)
head(pair.new,7)
head(pair.types,7)
###
theta.des  <- rbind( c(rbind(c(1,0),  c(1,0),  c(0,1),  c(0,0))),
		     c(rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))))
random.des <- rbind( 
        c(1,0,1,0),c(0,1,1,0),
        c(1,1,0,1),c(1,0,1,1))
mf <- 1*(pair.types[,1]=="mother" & pair.types[,2]=="father")
##          pair, rv related to pairs,  theta.des related to pair 
pairs.new <- cbind(pair.new,(mf==1)*1+(mf==0)*3,(mf==1)*2+(mf==0)*4,(mf==1)*1+(mf==0)*2,(mf==1)*3+(mf==0)*4)

## -----------------------------------------------------------------------------
# 3 rvs here 
random.des
theta.des

head(pairs.new)

## -----------------------------------------------------------------------------
tsdid2 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster, theta=c(2,1),
           random.design=random.des,theta.des=theta.des,pairs=pairs.new,dim.theta=2)
summary(tsdid2)

## -----------------------------------------------------------------------------
kinship  <- rep(0.5,nrow(pair.types))
kinship[pair.types[,1]=="mother" & pair.types[,2]=="father"] <- 0
head(kinship,n=10)

out <- make.pairwise.design(pair.new,kinship,type="ace") 

## -----------------------------------------------------------------------------
tsdid3 <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
             theta=c(2,1)/9,random.design=out$random.design,
             theta.des=out$theta.des,pairs=out$new.pairs,dim.theta=2)
summary(tsdid3)

## -----------------------------------------------------------------------------
library(mets)
set.seed(1000)
data <- simbinClaytonOakes.family.ace(500,2,1,beta=NULL,alpha=NULL)
head(data)
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=20)
pairs <- mm$pairs
dim(pairs)
head(pairs,12)

## -----------------------------------------------------------------------------
 dtypes <- interaction( data[pairs[,1],"type"], data[pairs[,2],"type"])
 dtypes <- droplevels(dtypes)
 table(dtypes)
 dm <- model.matrix(~-1+factor(dtypes))

## -----------------------------------------------------------------------------
aa <- glm(ybin~x,data=data,family=binomial())

tsp <- binomial.twostage(aa,data=data, clusters=data$cluster,
		 theta.des=dm,pairs=cbind(pairs,1:nrow(dm)))
summary(tsp)

## -----------------------------------------------------------------------------
tdp <-cbind( dataid[pair.new[,1],],dataid[pair.new[,2],])
names(tdp) <- c(paste(names(dataid),"1",sep=""),
		paste(names(dataid),"2",sep=""))
tdp <-transform(tdp,tt=interaction(type1,type2))
dlevel(tdp)
drelevel(tdp,newlevels=list(mother.father=4:9)) <-  obs.types~tt
dtable(tdp,~tt+obs.types)
tdp <- model.matrix(~-1+factor(obs.types),tdp)

## -----------------------------------------------------------------------------
###porpair <- binomial.twostage(aa,data=dataid,clusters=dataid$cluster,
###           theta.des=tdp,pairs=pair.new,model="or",var.link=1)
###summary(porpair)

## -----------------------------------------------------------------------------
sessionInfo()

