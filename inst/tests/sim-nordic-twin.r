
library(timereg)

F1addfg<-function(t,lam0=0.5,beta=c(-0.5,-0.005,-0.004),x=rep(0,3)) # FG
{ ## {{{
baset <- 0.13*pnorm((t-.70)/.13)
xm <- matrix(x,ncol=3) 
return( 1 - exp(-baset*exp(xm %*% matrix(beta,3,1)))) 
} ## }}}

corsim.prostate <- function(n,theta=1,thetaslope=0,crate=2,test=0,pcens=0,mt=1, same.cens=TRUE) 
{ ## {{{
###n <- 10; theta <- 1; thetaslope <- 0; mt <- 1
xl <- sample(1:4,n,replace=TRUE)
###xl <- rep(xl,each=2)
x<-cbind(xl==2,xl==3,xl==4)*1
tt<-seq(0,mt,length=mt*100)
###
###n=100;theta=1;lam0=0.5;beta=0.3;crate=2
thetat <- exp(log(theta))
F11x<-F1addfg(mt,x=x)
F12x<-F1addfg(mt,x=x)
###
thetaslut <- exp(log(theta)+thetaslope*(mt-mt/2))
p11 <- thetaslut*F11x*F12x/((1-F11x)+thetaslut*F11x)
p12 <- F11x-p11
p21 <- F12x-p11
p22 <- 1-F12x-F11x+p11
###apply(cbind(p11,p12,p21),1,sum)
if (test==1) { ## {{{
for (i in 1:2) {
print(x[i,]); 
F11xt<-F1addfg(tt,x=x[i,])
F12xt<-F1addfg(tt,x=x[i,])
p11t <- thetat* F11xt*F12xt/((1-F11xt)+thetat*F11xt)
cortt <- ((p11t)/(F12xt-p11t))/(F11xt/(1-F11xt))
###plot(tt,log(cortt))
if (i==1) { 
plot(tt,p11t,type="l",ylim=c(0,0.1),xlim=c(0,mt))
###lines(tt,F11x[i]-p11t,col=2)
###lines(tt,F12x[i]-p11t,col=2)
} else lines(tt,p11t,col=2);
###if (sum(diff(p11t<0))>0) stop("dec\n"); 
###p11 <- max(p11t)
###p12 <- F11x[i]-p11
###p21 <- F12x[i]-p11
###p22 <- 1- F12x[i]-F11x[i]+p11
###pnn <- 1- F12x[i]-F11x[i]+p11
}
} ## }}}
###apply(cbind(p11,p12,p21,p22),1,sum)
###
types <- rep(0,n)
causes <- matrix(0,n,2)
stime<-matrix(mt+1,n,2); 
for (i in 1:n)
{
ptype <- runif(1)
if (ptype<=p11[i]) {
	types[i] <- 1
	myhazx<-F1addfg(tt,x=x[i,])/F12x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist\n"); 
	stime[i,2]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	f1<- F1addfg(tt,x=x[i,])
       	myhazx<- (F12x[i]/p11[i]) * (thetat*f1/((1-f1)+thetat*f1))
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist\n"); 
	stime[i,1]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(1,1)
}
if ((ptype>p11[i]) & (ptype<=p12[i]+p11[i])) {
	types[i] <- 2
	f1 <- F1addfg(tt,x=x[i,])
	myhazx<- ( f1 - thetat*F12x[i]*f1/((1-f1)+thetat*f1))/p12[i]; 
	myhazx <- f1/F11x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist 2 \n"); 
	stime[i,1]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(1,2)
	stime[i,2] <- runif(1)*mt
}
if ((ptype>p11[i]+p12[i]) && (ptype<=p21[i]+p12[i]+p11[i])) {
	types[i] <- 3
        f2 <- F1addfg(tt,x=x[i,])
	myhazx <-  (f2 - (thetat*F11x[i]*f2/((1-F11x[i])+thetat*F11x[i])))/p21[i]; 
	myhazx <- f2/F12x[i]
###	if (abs(max(myhazx)-1)> 0.001) stop("not dist3 \n"); 
	stime[i,2]<-Cpred(cbind(myhazx,tt),runif(1))[1,2]+runif(1,0,0.001)
	causes[i,] <- c(2,1)
	stime[i,1] <- runif(1)*mt
}
if (ptype>p11[i]+p12[i]+p21[i] ) {
	types[i] <- 4
	causes[i,] <- c(2,2)
	stime[i,1:2] <- runif(2)*mt
}

}
###stime
###causes
stime <- c(t(stime))
cause <- c(t(causes))

###same.cens=TRUE
if (same.cens==TRUE) {
	ctime <- rep(rbinom(n/2,1,pcens),each=2) 
        ctime[ctime==1] <- rep(runif(sum(ctime==1)/2),each=2)*(mt+2)
}
else {
	ctime<- rbinom(n,1,pcens)
        ctime[ctime==1] <- runif(sum(ctime==1))*(mt+2)
}

ctime[ctime==0] <- mt;

cens <- (ctime< stime)
time <- ifelse(cens,ctime,stime)
cause <- ifelse(cens,0,cause)
id <- rep(1:n,rep(2,n))

country <- c()
country[xl==1] <- "SWE"
country[xl==2] <- "DK"
country[xl==3] <- "FIN"
country[xl==4] <- "NOR"

data<-data.frame(time=time,cause=cause,xl=rep(xl,each=2),
		 country=rep(country,each=2),id=id,cens=cens,stime=stime,type=rep(types,each=2),
		 f1=rep(F11x,each=2),p11=rep(p11,each=2),p12=rep(p12,each=2),p21=rep(p21,each=2),
		 p22=rep(p22,each=2))
return(data)
} ## }}}

simnordic <- function(n) 
{

outdz <- corsim.prostate(n,theta=1.7,crate=2,pcens=0.8,mt=1,same.cens=TRUE,test=0) 
outmz <- corsim.prostate(n,theta=3,crate=2,pcens=0.8,mt=1,same.cens=TRUE,test=0) 
outdz$zyg <- 0
outmz$zyg <- 1
outmz$id <- outmz$id+nrow(outdz)
###
out <- rbind(outdz,outmz)
out$time <- out$time*100
table(out$type,out$country)
table(out$type,out$cause)
out$country <- relevel(factor(out$country),ref="SWE")
table(out$country)
outk <- out[,c("country","cause","id","time","zyg","type")]
table(outk$cause)
table(outk$type,outk$country)
table(outk$cause,outk$country)
###

###times <- seq(40,100,length=30)
###adds <- comp.risk(Surv(time,cause==0)~const(country)+cluster(id),data=outk,outk$cause,
###	  causeS=1,conservative=1,times=times,max.clust=NULL,n.sim=0,cens.model="aalen",model="fg")
###summary(adds)
###plot(adds)
######adds$cum
######
###padds <- predict(adds,X=rep(1,2),Z=rbind(c(1,0,0),rep(0,3)),se=0,uniform=0)
###plot(padds,ylim=c(0,0.15),multiple=1,col=1:2)
###padds$P1
return(outk)
}

outk <- simnordic(5000)

write.table(outk,"nordic-prostate.txt")

