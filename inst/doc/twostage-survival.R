## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  dpi=50,
  fig.width=7.15, fig.height=5.5,
  out.width="600px",
  fig.retina=1,
  comment = "#>"
)
library(mets)

## -----------------------------------------------------------------------------
 library(mets)
 data(diabetes)
 set.seed(100)
 
 # Marginal Cox model  with treat as covariate
 margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 # Clayton-Oakes, MLE 
 fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
 summary(fitco1)
 
 # Clayton-Oakes
 fitco2 <- survival.twostage(margph,data=diabetes,theta=0.0,
                  clusters=diabetes$id,var.link=1,model="clayton.oakes")
 summary(fitco2)
 fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,
                  clusters=diabetes$id,var.link=0,model="clayton.oakes")
 summary(fitco3)

## -----------------------------------------------------------------------------
  # without covariates but marginal model stratified 
  marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
 fitco<-twostageMLE(marg,data=diabetes,theta=1.0)
 summary(fitco)

  fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
		   model="clayton.oakes",var.link=0)
  summary(fitcoa)

## -----------------------------------------------------------------------------
 d <- simClaytonOakes(200,2,0.5,0,3)
  margph <- phreg(Surv(time,status)~x+cluster(cluster),data=d)
 # Clayton-Oakes, MLE 
 fitco1<-twostageMLE(margph,data=d)
 summary(fitco1)

## -----------------------------------------------------------------------------
 udp <- piecewise.twostage(c(0,0.5,2),data=d,id="cluster",timevar="time",status="status",model="clayton.oakes",silent=0)
 summary(udp)

## -----------------------------------------------------------------------------
 data <- simClaytonOakes.twin.ace(200,2,1,0,3)

 out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 pardes <- out$pardes
 pardes 

## -----------------------------------------------------------------------------
 des.rv <- out$des.rv
 # MZ
 head(des.rv,2)
 # DZ 
 tail(des.rv,2)

## -----------------------------------------------------------------------------
### data <- simClaytonOakes.twin.ace(2000,2,1,0,3)
### out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
 ts <- twostage(aa,data=data,clusters=data$cluster,
      theta=c(2,1),var.link=0,random.design=out$des.rv,theta.des=out$pardes)
 summary(ts)

## -----------------------------------------------------------------------------
run <- 0
if (run==1) {
 data <- simClaytonOakes.twin.ace(1000000,2,1,0,3)

 out <- twin.polygen.design(data,id="cluster",zyg="DZ",zygname="zyg",type="ace")
 pardes <- out$pardes
 aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data)
 system.time(
 ts <- twostage(aa,data=data,clusters=data$cluster,
      theta=c(2,1),var.link=0,random.design=out$des.rv,theta.des=out$pardes)
 )
 summary(ts)
}

## -----------------------------------------------------------------------------
kendall.ClaytonOakes.twin.ace(ts$theta[1],ts$theta[2],K=10000) 

## -----------------------------------------------------------------------------
library(mets)
set.seed(1000)
data <- simClaytonOakes.family.ace(200,2,1,0,3)
head(data)
data$number <- c(1,2,3,4)
data$child <- 1*(data$number==3)

## -----------------------------------------------------------------------------
out <- ace.family.design(data,member="type",id="cluster")
out$pardes
head(out$des.rv,4)

## -----------------------------------------------------------------------------
pa <- phreg(Surv(time,status)~+1+cluster(cluster),data=data)

# make ace random effects design
ts <- twostage(pa,data=data,clusters=data$cluster,var.par=1,var.link=0,theta=c(2,1),
        random.design=out$des.rv,theta.des=out$pardes)
summary(ts)

## -----------------------------------------------------------------------------
# now specify fitting via specific pairs 
# first all pairs 
mm <- familycluster.index(data$cluster)
head(mm$familypairindex,n=12)
pairs <- matrix(mm$familypairindex,ncol=2,byrow=TRUE)
head(pairs,n=6)

## -----------------------------------------------------------------------------
ts <- twostage(pa,data=data,clusters=data$cluster, theta=c(2,1),var.link=0,step=1.0,
        random.design=out$des.rv, theta.des=out$pardes,pairs=pairs)
summary(ts)

## -----------------------------------------------------------------------------
ssid <- sort(sample(1:nrow(pairs),200))
tsd <- twostage(pa,data=data,clusters=data$cluster,
    theta=c(2,1)/10,var.link=0,random.design=out$des.rv,
   theta.des=out$pardes,pairs=pairs[ssid,])
summary(tsd)

## -----------------------------------------------------------------------------
ids <- sort(unique(c(pairs[ssid,])))

pairsids <- c(pairs[ssid,])
pair.new <- matrix(fast.approx(ids,c(pairs[ssid,])),ncol=2)
head(pair.new)

# this requires that pair.new refers to id's in dataid (survival, status and so forth)
# random.design and theta.des are constructed to be the array 3 dims via individual specfication from ace.family.design
dataid <- dsort(data[ids,],"cluster")
outid <- ace.family.design(dataid,member="type",id="cluster")
outid$pardes
head(outid$des.rv)

## -----------------------------------------------------------------------------
tsdid <- twostage(pa,data=dataid,clusters=dataid$cluster,theta=c(2,1)/10,var.link=0,baseline.iid=0,
          random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdid)

paid <- phreg(Surv(time,status)~+1+cluster(cluster),data=dataid)
tsdidb <- twostage(paid,data=dataid,clusters=dataid$cluster,theta=c(2,1)/10,
  var.link=0,random.design=outid$des.rv,theta.des=outid$pardes,pairs=pair.new)
summary(tsdidb)
coef(tsdid)

## -----------------------------------------------------------------------------
pair.types <-  matrix(dataid[c(t(pair.new)),"type"],byrow=T,ncol=2)
head(pair.new)
head(pair.types)

theta.des  <- rbind( c(rbind(c(1,0),c(1,0),c(0,1),c(0,0))),
		c(rbind(c(0.5,0),c(0.5,0),c(0.5,0),c(0,1))))
random.des <- rbind( 
        c(1,0,1,0),c(0,1,1,0),
        c(1,1,0,1),c(1,0,1,1))
mf <- 1*(pair.types[,1]=="mother" & pair.types[,2]=="father")
##          pair, rv related to pairs,  theta.des related to pair 
pairs.new <- cbind(pair.new,(mf==1)*1+(mf==0)*3,(mf==1)*2+(mf==0)*4,(mf==1)*1+(mf==0)*2,(mf==1)*3+(mf==0)*4)

## -----------------------------------------------------------------------------
head(pairs.new[1:3,])
head(dataid)

## -----------------------------------------------------------------------------
random.des[1,]
random.des[2,]
matrix(theta.des[1,],4,2)

## -----------------------------------------------------------------------------
head(dataid)
matrix(theta.des[2,],4,2)
random.des[3,]
random.des[4,]

## -----------------------------------------------------------------------------
tsdid2 <- twostage(pa,data=dataid,clusters=dataid$cluster,
       theta=c(2,1)/10,var.link=0,step=1.0,random.design=random.des,
       baseline.iid=0, theta.des=theta.des,pairs=pairs.new,dim.theta=2)
summary(tsdid2)
tsd$theta
tsdid2$theta
tsdid$theta

## -----------------------------------------------------------------------------
kinship  <- rep(0.5,nrow(pair.types))
kinship[pair.types[,1]=="mother" & pair.types[,2]=="father"] <- 0
head(kinship,n=10)

out <- make.pairwise.design(pair.new,kinship,type="ace") 

## -----------------------------------------------------------------------------
tsdid3 <- twostage(pa,data=dataid,clusters=dataid$cluster,
   theta=c(2,1)/10,var.link=0,step=1.0,random.design=out$random.design,
   baseline.iid=0,theta.des=out$theta.des,pairs=out$new.pairs,dim.theta=2)
summary(tsdid3)
tsdid2$theta
tsdid$theta

## -----------------------------------------------------------------------------
outae <- make.pairwise.design(pair.new,kinship,type="ae") 
tsdid4 <- twostage(pa,data=dataid,clusters=dataid$cluster,
   theta=c(2,1)/10,var.link=0,random.design=outae$random.design,
   baseline.iid=0,theta.des=outae$theta.des,pairs=outae$new.pairs,dim.theta=1)
summary(tsdid4)

## -----------------------------------------------------------------------------
library(mets)
set.seed(1000)
data <- simClaytonOakes.family.ace(200,2,1,0,3)
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
 head(dm)

## -----------------------------------------------------------------------------
pa <- phreg(Surv(time,status)~cluster(cluster),data)

tsp <- twostage(pa,data=data,theta.des=dm,pairs=cbind(pairs,1:nrow(dm)),se.clusters=data$clust)
summary(tsp)

## -----------------------------------------------------------------------------
 library(mets)
 data(diabetes)
 
 # Marginal Cox model  with treat as covariate
 margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 # Clayton-Oakes, MLE 
 fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
 summary(fitco1)
 
 # Plackett model
 mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
 fitp <- survival.twostage(mph,data=diabetes,theta=3.0,
                clusters=diabetes$id,var.link=1,model="plackett")
 summary(fitp)
 
 # without covariates but with stratafied 
 marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
 fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
                 clusters=diabetes$id)
 summary(fitpa)
 
 fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
                  model="clayton.oakes")
 summary(fitcoa)

## -----------------------------------------------------------------------------
 mm <- model.matrix(~-1+factor(adult),diabetes)
 fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
                clusters=diabetes$id,var.link=1,model="plackett",
		theta.des=mm)
 summary(fitp)

## -----------------------------------------------------------------------------
 # Piecewise constant cross hazards ratio modelling

 d <- subset(simClaytonOakes(1000,2,0.5,0,stoptime=2,left=0),!truncated)
 udp <- piecewise.twostage(c(0,0.5,2),data=d,id="cluster",timevar="time",
                           status="status",model="plackett",silent=0)
 summary(udp)

## -----------------------------------------------------------------------------
sessionInfo()

