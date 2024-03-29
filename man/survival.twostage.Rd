% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/twostage.R
\name{survival.twostage}
\alias{survival.twostage}
\alias{twostage.aalen}
\alias{twostage.cox.aalen}
\alias{twostage.coxph}
\alias{twostage.phreg}
\alias{randomDes}
\alias{readmargsurv}
\title{Twostage survival model for multivariate survival data}
\usage{
survival.twostage(
  margsurv,
  data = parent.frame(),
  method = "nr",
  detail = 0,
  clusters = NULL,
  silent = 1,
  weights = NULL,
  theta = NULL,
  theta.des = NULL,
  var.link = 1,
  baseline.iid = 1,
  model = "clayton.oakes",
  marginal.trunc = NULL,
  marginal.survival = NULL,
  strata = NULL,
  se.clusters = NULL,
  numDeriv = 1,
  random.design = NULL,
  pairs = NULL,
  dim.theta = NULL,
  numDeriv.method = "simple",
  additive.gamma.sum = NULL,
  var.par = 1,
  no.opt = FALSE,
  ...
)
}
\arguments{
\item{margsurv}{Marginal model}

\item{data}{data frame}

\item{method}{Scoring method "nr", for lava NR optimizer}

\item{detail}{Detail}

\item{clusters}{Cluster variable}

\item{silent}{Debug information}

\item{weights}{Weights}

\item{theta}{Starting values for variance components}

\item{theta.des}{design for dependence parameters, when pairs are given the indeces of the
theta-design for this pair, is given in pairs as column 5}

\item{var.link}{Link function for variance:  exp-link.}

\item{baseline.iid}{to adjust for baseline estimation, using phreg function on same data.}

\item{model}{model}

\item{marginal.trunc}{marginal left truncation probabilities}

\item{marginal.survival}{optional vector of marginal survival probabilities}

\item{strata}{strata for fitting, see example}

\item{se.clusters}{for clusters for se calculation with iid}

\item{numDeriv}{to get numDeriv version of second derivative, otherwise uses sum of squared scores for each pair}

\item{random.design}{random effect design for additive gamma model, when pairs are given the
indeces of the pairs random.design rows are given as columns 3:4}

\item{pairs}{matrix with rows of indeces (two-columns) for the pairs considered in the pairwise
composite score, useful for case-control sampling when marginal is known.}

\item{dim.theta}{dimension of the theta parameter for pairs situation.}

\item{numDeriv.method}{uses simple to speed up things and second derivative not so important.}

\item{additive.gamma.sum}{for two.stage=0, this is specification of the lamtot in the models via
a matrix that is multiplied onto the parameters theta (dimensions=(number random effects x number
of theta parameters), when null then sums all parameters.}

\item{var.par}{is 1 for the default parametrization with the variances of the random effects,
var.par=0 specifies that the \eqn{\lambda_j}'s are used as parameters.}

\item{no.opt}{for not optimizng}

\item{...}{Additional arguments to maximizer NR of lava.
and ascertained sampling}
}
\description{
Fits Clayton-Oakes or bivariate Plackett models for bivariate survival data
using marginals that are on Cox form. The dependence can be modelled via
\enumerate{
\item  Regression design on dependence parameter.
\item  Random effects, additive gamma model.
}

If clusters contain more than two subjects, we use a composite likelihood
based on the pairwise bivariate models, for full MLE see twostageMLE.

The two-stage model is constructed such that
given the gamma distributed random effects it is assumed that the survival functions
are indpendent, and that the marginal survival functions are on Cox form (or additive form)
\deqn{
P(T > t| x) = S(t|x)= exp( -exp(x^T \beta) A_0(t) )
}

One possibility is to model the variance within clusters via a regression design, and
then one can specify a regression structure for the independent gamma distributed
random effect for each cluster, such that the variance is given by
\deqn{
 \theta = h( z_j^T \alpha)
}
where \eqn{z} is specified by theta.des, and a possible link function var.link=1 will
will use the exponential link \eqn{h(x)=exp(x)}, and var.link=0 the identity link \eqn{h(x)=x}.
The reported standard errors are based on the estimated information from the
likelihood assuming that the marginals are known (unlike the twostageMLE and for the
additive gamma model below).

Can also fit a structured additive gamma random effects model, such
as the ACE, ADE model for survival data.  In this case the
random.design specificies the random effects for each subject within a cluster. This is
a matrix of 1's and 0's with dimension n x d.  With d random effects.
For a cluster with two subjects, we let the random.design rows be
 \eqn{v_1} and \eqn{v_2}.
Such that the random effects for subject
1 is \deqn{v_1^T (Z_1,...,Z_d)}, for d random effects. Each random effect
has an associated parameter \eqn{(\lambda_1,...,\lambda_d)}.
By construction subjects 1's random effect are Gamma distributed with
mean \eqn{\lambda_j/v_1^T \lambda}
and variance \eqn{\lambda_j/(v_1^T \lambda)^2}. Note that the random effect
\eqn{v_1^T (Z_1,...,Z_d)} has mean 1 and variance \eqn{1/(v_1^T \lambda)}.
It is here asssumed that  \eqn{lamtot=v_1^T \lambda} is fixed within clusters
as it would be for the ACE model below.

Based on these parameters the relative contribution (the heritability, h) is
equivalent to  the expected values of the random effects: \eqn{\lambda_j/v_1^T \lambda}

The DEFAULT parametrization (var.par=1) uses the variances of the random effecs
\deqn{
\theta_j  = \lambda_j/(v_1^T \lambda)^2
}
For alternative parametrizations one can specify how the parameters relate to \eqn{\lambda_j}
with the argument var.par=0.

For both types of models the basic model assumptions are that
given the random effects of the clusters the survival distributions within a cluster
are independent and ' on the form
\deqn{
P(T > t| x,z) = exp( -Z \cdot Laplace^{-1}(lamtot,lamtot,S(t|x)) )
}
with the inverse laplace of the gamma distribution with mean 1 and variance 1/lamtot.

The parameters \eqn{(\lambda_1,...,\lambda_d)} are related to the parameters of the model
by a regression construction \eqn{pard} (d x k), that links the \eqn{d}
\eqn{\lambda} parameters
with the (k) underlying \eqn{\theta} parameters
\deqn{
\lambda = theta.des  \theta
}
here using theta.des to specify these low-dimension association. Default is a diagonal matrix.
This can be used to make structural assumptions about the variances of the random-effects
as is needed for the ACE model for example.

The case.control option that can be used with the pair specification of the pairwise parts
of the estimating equations. Here it is assumed that the second subject of each pair is the proband.
}
\examples{
data(diabetes)

# Marginal Cox model  with treat as covariate
margph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
### Clayton-Oakes, MLE
fitco1<-twostageMLE(margph,data=diabetes,theta=1.0)
summary(fitco1)

### Plackett model
mph <- phreg(Surv(time,status)~treat+cluster(id),data=diabetes)
fitp <- survival.twostage(mph,data=diabetes,theta=3.0,Nit=40,
               clusters=diabetes$id,var.link=1,model="plackett")
summary(fitp)

### Clayton-Oakes
fitco2 <- survival.twostage(mph,data=diabetes,theta=0.0,detail=0,
                 clusters=diabetes$id,var.link=1,model="clayton.oakes")
summary(fitco2)
fitco3 <- survival.twostage(margph,data=diabetes,theta=1.0,detail=0,
                 clusters=diabetes$id,var.link=0,model="clayton.oakes")
summary(fitco3)

### without covariates but with stratafied
marg <- phreg(Surv(time,status)~+strata(treat)+cluster(id),data=diabetes)
fitpa <- survival.twostage(marg,data=diabetes,theta=1.0,
                clusters=diabetes$id,model="clayton.oakes")
summary(fitpa)

fitcoa <- survival.twostage(marg,data=diabetes,theta=1.0,clusters=diabetes$id,
                 model="clayton.oakes")
summary(fitcoa)

### Piecewise constant cross hazards ratio modelling
########################################################

d <- subset(simClaytonOakes(2000,2,0.5,0,stoptime=2,left=0),!truncated)
udp <- piecewise.twostage(c(0,0.5,2),data=d,method="optimize",
                          id="cluster",timevar="time",
                          status="status",model="clayton.oakes",silent=0)
summary(udp)

\donttest{ ## Reduce Ex.Timings
### Same model using the strata option, a bit slower
########################################################
## makes the survival pieces for different areas in the plane
##ud1=surv.boxarea(c(0,0),c(0.5,0.5),data=d,id="cluster",timevar="time",status="status")
##ud2=surv.boxarea(c(0,0.5),c(0.5,2),data=d,id="cluster",timevar="time",status="status")
##ud3=surv.boxarea(c(0.5,0),c(2,0.5),data=d,id="cluster",timevar="time",status="status")
##ud4=surv.boxarea(c(0.5,0.5),c(2,2),data=d,id="cluster",timevar="time",status="status")

## everything done in one call
ud <- piecewise.data(c(0,0.5,2),data=d,timevar="time",status="status",id="cluster")
ud$strata <- factor(ud$strata);
ud$intstrata <- factor(ud$intstrata)

## makes strata specific id variable to identify pairs within strata
## se's computed based on the id variable across strata "cluster"
ud$idstrata <- ud$id+(as.numeric(ud$strata)-1)*2000

marg2 <- timereg::aalen(Surv(boxtime,status)~-1+factor(num):factor(intstrata),
               data=ud,n.sim=0,robust=0)
tdes <- model.matrix(~-1+factor(strata),data=ud)
fitp2 <- survival.twostage(marg2,data=ud,se.clusters=ud$cluster,clusters=ud$idstrata,
                model="clayton.oakes",theta.des=tdes,step=0.5)
summary(fitp2)

### now fitting the model with symmetry, i.e. strata 2 and 3 same effect
ud$stratas <- ud$strata;
ud$stratas[ud$strata=="0.5-2,0-0.5"] <- "0-0.5,0.5-2"
tdes2 <- model.matrix(~-1+factor(stratas),data=ud)
fitp3 <- survival.twostage(marg2,data=ud,clusters=ud$idstrata,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5)
summary(fitp3)

### same model using strata option, a bit slower
fitp4 <- survival.twostage(marg2,data=ud,clusters=ud$cluster,se.cluster=ud$cluster,
                model="clayton.oakes",theta.des=tdes2,step=0.5,strata=ud$strata)
summary(fitp4)
}

\donttest{ ## Reduce Ex.Timings
### structured random effects model additive gamma ACE
### simulate structured two-stage additive gamma ACE model
data <- simClaytonOakes.twin.ace(4000,2,1,0,3)
out <- twin.polygen.design(data,id="cluster")
pardes <- out$pardes
pardes
des.rv <- out$des.rv
head(des.rv)
aa <- phreg(Surv(time,status)~x+cluster(cluster),data=data,robust=0)
ts <- survival.twostage(aa,data=data,clusters=data$cluster,detail=0,
	       theta=c(2,1),var.link=0,step=0.5,
	       random.design=des.rv,theta.des=pardes)
summary(ts)
}

}
\references{
Twostage estimation of additive gamma frailty models for survival data.
Scheike (2019), work in progress

Shih and Louis (1995) Inference on the association parameter in copula models for bivariate
survival data, Biometrics, (1995).

Glidden (2000), A Two-Stage estimator of the dependence
parameter for the Clayton Oakes model, LIDA, (2000).

Measuring early or late dependence for bivariate twin data
Scheike, Holst, Hjelmborg (2015), LIDA

Estimating heritability for cause specific mortality based on twins studies
Scheike, Holst, Hjelmborg (2014), LIDA

Additive Gamma frailty models for competing risks data, Biometrics (2015)
Eriksson and Scheike (2015),
}
\author{
Thomas Scheike
}
\keyword{survival}
