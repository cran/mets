#ifndef MVN_H
#define MVN_H

#include <RcppArmadillo.h>
#include <Rmath.h>
#include <cmath>
/* #include <iostream> */
/* #include <sstream> */
/* #include <cstring> */


using namespace std;
using namespace Rcpp;
using namespace arma;


const double twopi = 2*M_PI; //datum::pi;
const double itwopi = 0.5/M_PI; //datum::pi;
const double log2pi = log(twopi);

struct vecmat
{
  vec V;
  mat M;
};


RcppExport SEXP pmvn(SEXP lower, SEXP upper, SEXP mu, SEXP sigma, SEXP cor);

RcppExport SEXP Dpmvn(SEXP lower, SEXP upper, SEXP mu, SEXP sigma, SEXP std);

/* RcppExport SEXP loglikMVN(SEXP yl, SEXP yu,  */
/* 			  SEXP status, */
/* 			  SEXP mu, SEXP dmu, */
/* 			  SEXP s, SEXP ds, */
/* 			  SEXP z, SEXP su, SEXP dsu, */
/* 			  SEXP threshold, SEXP dthreshold, */
/* 			  SEXP score); */


RcppExport SEXP bvncdf(SEXP a, SEXP b, SEXP r);

double Sbvn(double &l1, double &l2,double &r);

inline double Fbvn(double u1, double u2, double r) { 
  u1 *= -1; u2 *= -1; 
  return(Sbvn(u1,u2,r)); 
}

vecmat Dbvn(double y1, double y2, double R);

double dbvnorm(double y1, double y2, double R);

double mvtdst(int* n,
	      int* nu,
	      double* lower,
	      double* upper,
	      int* infin,
	      double* correl,
	      double* delta,
	      int* maxpts,  
	      double* abseps,
	      double* releps,
	      double* error, 
	      double* value, 
	      int* inform);  

extern "C" double bvnd_(const double *dh, const double *dk, const double *r);

#endif /* MVN_H */
