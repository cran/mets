##' @export
`lifetable` <- function(x,...) UseMethod("lifetable")
    

##' Create simple life table 
##' 
##' @title Life table
##' @param x time formula (Surv) or matrix/data.frame with columns time,status or entry,exit,status
##' @param strata Strata
##' @param data data.frame
##' @param breaks Time intervals
##' @param confint If TRUE 95% confidence limits are calculated
##' @param ... Additional arguments to lower level functions
##' @author Klaus K. Holst
##' @aliases lifetable lifetable.matrix lifetable.formula
##' @usage
##'  \method{lifetable}{matrix}(x, strata = list(), breaks = c(),
##'    confint = FALSE, ...)
##'
##'  \method{lifetable}{formula}(x, data=parent.frame(), breaks = c(),
##'    confint = FALSE, ...)
##' @examples
##' library(timereg)
##' data(TRACE)
##' \donttest{
##'     lifetable(Surv(time,status==9)~sex+I(cut(wmi,c(-Inf,1,1.5,Inf))),
##'               data=TRACE,breaks=c(0.2,0.5),confint=TRUE)
##' }
##' 
##' d <- with(TRACE,lifetable(Surv(time,status==9)~sex+vf,breaks=c(0.2,0.5)))
##' summary(glm(events ~ offset(log(atrisk))+interval*vf + sex*vf,
##'             data=d,poisson))
##' @S3method lifetable matrix
lifetable.matrix <- function(x,strata=list(),breaks=c(),confint=FALSE,...) {
    if (ncol(x)==3) {
        status <- x[,3]
        entry <- x[,1]
        time <- x[,2]
    } else {
        status <- x[,2]
        time <- x[,1]
        entry <- rep(0,length(time))
    }
    LifeTable(time,status,entry,strata,breaks,confint,...)
}

##' @S3method lifetable formula
lifetable.formula <- function(x,data=parent.frame(),breaks=c(),confint=FALSE,...) {
    cl <- match.call()
    mf <- model.frame(x,data)
    Y <- model.extract(mf, "response")
    Terms <- terms(x, data = data)
    if (!is.Surv(Y)) stop("Expected a 'Surv'-object")
    if (ncol(Y)==2) {
        exit <- Y[,1]
        entry <- NULL ## rep(0,nrow(Y))
        status <- Y[,2]
    } else {
        entry <- Y[,1]
        exit <- Y[,2]
        status <- Y[,3]
    }
    strata <- list()
    X <- model.matrix(Terms, data)
    if (!is.null(intpos  <- attributes(Terms)$intercept))
        X <- X[,-intpos,drop=FALSE]
    if (ncol(X)>0) {
        strata <- as.list(model.frame(Terms,data)[,-1,drop=FALSE])
    }
    LifeTable(exit,status,entry,strata,breaks,confint,...)       
}


LifeTable <- function(time,status,entry=NULL,strata=list(),breaks=c(),confint=FALSE) {
    if (is.null(entry)) entry <- rep(0,NROW(time))
    if ((is.matrix(time) || is.data.frame(time)) && ncol(time)>1) {
        if (ncol(time)==3) {
            status <- time[,3]
            entry <- time[,1]
            time <- time[,2]
        } else {
            status <- time[,2]
            time <- time[,1]
            entry <- rep(0,length(time))
        }
    }
    if (length(strata)>0) {
        a <- by(cbind(entry,time,status), strata,
                FUN=LifeTable, breaks=breaks, confint=confint)
        nn <- do.call("expand.grid",attributes(a)$dimnames)
        nam <- nn[rep(seq(NROW(nn)),each=NROW(a[[1]])),,drop=FALSE]
        res <- cbind(Reduce("rbind",a),nam)
        return(res)
    }    
    if (length(breaks)>0) breaks <- sort(unique(breaks))
    en <- matrix(unlist(lapply(c(-Inf, breaks),function(x) pmax(x,entry))),
                  ncol=length(breaks)+1)
    ex <- matrix(unlist(lapply(c(breaks, Inf),function(x) pmin(x,time))),
                  ncol=length(breaks)+1)    
    dur <- ex-en
    endur <- rbind(c(breaks,Inf))%x%cbind(rep(1,nrow(en)))-en
    dur[dur<0] <- NA
    enter <- colSums(!is.na(dur))
    atrisk <- colSums(dur,na.rm=TRUE)
    eventcens <- dur<endur
    eventcens <- rbind(apply(dur<endur,2,function(x) x*(status+1)))
    lost <- colSums(eventcens==1,na.rm=TRUE)
    events <- colSums(eventcens==2,na.rm=TRUE)
    res <- data.frame(enter=enter,
                      atrisk=atrisk,
                      lost=lost,
                      events=events,
                      interval=as.factor(c(breaks,Inf)),
                      rate=events/atrisk)
    if (confint) {
        f <- events ~ offset(log(atrisk))
        if (length(breaks)>0) f <- update(f,.~.+factor(interval)-1)
        g <- glm(f,data=res,poisson)
        ci <- rbind(exp(stats::confint(g)))
        res[,"2.5%"] <- ci[,1]
        res[,"97.5%"] <- ci[,2]
    }
    res
}