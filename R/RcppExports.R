# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

.ApplyBy2 <- function(idata, icluster, F, Env, Argument = "x", Columnwise = 0L, Reduce = 0L, epsilon = 1.0e-16) {
    .Call(`_mets_ApplyBy2`, idata, icluster, F, Env, Argument, Columnwise, Reduce, epsilon)
}

.ApplyBy <- function(idata, icluster, f) {
    .Call(`_mets_ApplyBy`, idata, icluster, f)
}

# Register entry points for exported C++ functions
methods::setLoadAction(function(ns) {
    .Call(`_mets_RcppExport_registerCCallable`)
})
