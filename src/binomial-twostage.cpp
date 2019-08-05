#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <complex.h>
#include "twostage.h"

using namespace arma;
using namespace Rcpp;

double claytonoakesP(double theta,int status1,int status2,double cif1,double cif2,vec &dp,vec &ps,vec &dp00) 
{ // {{{
double valr=1,x,y,z;
double p00,p10,p01,p11; 

//double cifs=cif1+cif2; //double S=1+(cifs*(theta-1)); 
x=theta; y=cif1; z=cif2; 

valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);

p11=valr; 
p10=cif1-p11; 
p01=cif2-p11; 
p00=1-cif1-cif2+p11; 

ps(0)=p00; ps(1)=p10; ps(2)=p01; ps(3)=p11; 
dp00(0)=dp(0); 
//printf(" %lf %lf %lf %lf %lf %lf %lf \n",theta,y,z,p11,p10,p01,p00); 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 

//printf(" %lf \n",valr); 

//if (status1==0 && status2==0) { // {{{
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//} // }}}
//
//if (status1==1 && status2==0) { // {{{
//valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==0 && status2==1) { // {{{
//valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==1 && status2==1) { // {{{
//valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
//dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
//} // }}}

return(valr); 
} // }}}

double placklikeP(double theta,int status1,int status2,double cif1,double cif2,vec &dp,vec &ps,vec &dp00) 
{ // {{{
//double S,S2,a;
//S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
double x,y,z,valr=1,p11,p10,p01,p00; 
//double cifs=cif1+cif2; 
//a=(1+(theta-1)*(cifs)); 
x=theta; y=cif1; z=cif2; 

dp(0)=0; 

if (theta!=1) {
p11=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
} else p11=cif1*cif2;

p10=y-p11; 
p01=z-p11; 
p00=1-y-z+p11; 
ps(0)=p00; ps(1)=p10; ps(2)=p01; ps(3)=p11; 
dp00(0)=dp(0); 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 

return(valr); 
} // }}}

//RcppExport SEXP claytonoakesPR(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2) 
//{ // {{{
// colvec theta = Rcpp::as<colvec>(itheta);
// colvec cif1 = Rcpp::as<colvec>(icif1);
// colvec cif2 = Rcpp::as<colvec>(icif2);
// colvec status1 = Rcpp::as<colvec>(iistatus1);
// colvec status2 = Rcpp::as<colvec>(iistatus2);
//
// colvec L=theta; 
// colvec dL=theta; 
// int n=cif1.size(); 
//
//double valr=1,x,y,z;
//double p00,p10,p01,p11; 
//
////double cifs=cif1+cif2; //double S=1+(cifs*(theta-1)); 
//x=theta; y=cif1; z=cif2; 
//
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//
//p11=valr; 
//p10=x-p11; 
//p01=y-p11; 
//p00=1-x-y+p11; 
//
//if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
//if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
//if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
//if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//
//if (status1==0 && status2==0) { // {{{
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//} // }}}
//
//if (status1==1 && status2==0) { // {{{
//valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==0 && status2==1) { // {{{
//valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==1 && status2==1) { // {{{
//valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
//dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
//} // }}}
//
//return(valr); 
//} // }}}
//
//RcppExport SEXP placklikePR(SEXP itheta,SEXP istatus1,SEXP istatus2,SEXP icif1,SEXP icif2)
//{ // {{{
////double S,S2,a;
////S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
//double x,y,z,valr=1,p11,p10,p01,p00; 
////double cifs=cif1+cif2; 
////a=(1+(theta-1)*(cifs)); 
//x=theta; y=cif1; z=cif2; 
//
//dp(0)=0; 
//
//if (theta!=1) {
//p11=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
//dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
//} else p11=cif1*cif2;
//
//p11=p11;
//p10=y-p11; 
//p01=z-p11; 
//p00=1-y-z+p11; 
//
//if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
//if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
//if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
//if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//
//return(valr); 
//} // }}}


// double CclaytonoakesP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
//{ // {{{
//double valr=1,x,y,z;
//double p00,p10,p01,p11; 
//
////double cifs=cif1+cif2; //double S=1+(cifs*(theta-1)); 
//x=theta; y=cif1; z=cif2; 
//
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//
//p11=valr; 
//p10=x-p11; 
//p01=y-p11; 
//p00=1-x-y+p11; 
//
//double epsilon=1E-20; 
//cx_double Ctheta,Cvalr,Cy,Cz; 
//Ctheta=cx_double(theta,epsilon); 
//Cy=cx_double(y,0); 
//Cz=cx_double(z,0); 
//
//printf(" mig \n"); 
////Cvalr=  pow((1/pow(Cy,1/Ctheta) + 1/pow(Cz,1/Ctheta)) - 1,-Ctheta);
//double dd=imag(Cvalr)/epsilon; 
//
//printf("complex  %lf ",dd); 
//
//if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
//if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
//if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
//if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//
//if (status1==0 && status2==0) { // {{{
//valr=  pow((1/pow(y,1/x) + 1/pow(z,1/x)) - 1,-x);
//dp(0)= (-((x*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x))) - log(-1 + pow(y,-1/x) + pow(z,-1/x)))/pow(-1 + pow(y,-1/x) + pow(z,-1/x),x);
//} // }}}
//
//if (status1==1 && status2==0) { // {{{
//valr=pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(y))/pow(x,2) + pow(y,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==0 && status2==1) { // {{{
//valr=pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x);
//dp(0)=(pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*log(z))/pow(x,2) + pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-1 - x)*(((-1 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x)));
//} // }}}
//
//if (status1==1 && status2==1) { // {{{
//valr= -(((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x);
//dp(0)=((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/pow(x,2) + (pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x))/x - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(y))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*log(z))/pow(x,3) - ((-1 - x)*pow(y,-1 - 1/x)*pow(z,-1 - 1/x)*pow(-1 + pow(y,-1/x) + pow(z,-1/x),-2 - x)*(((-2 - x)*(log(y)/(pow(x,2)*pow(y,1/x)) + log(z)/(pow(x,2)*pow(z,1/x))))/(-1 + pow(y,-1/x) + pow(z,-1/x)) - log(-1 + pow(y,-1/x) + pow(z,-1/x))))/x;
//} // }}}
//
//printf(" %lf \n",dp(0)); 
//return(valr); 
//} // }}}

cx_double Cpij(cx_double x, cx_double y, cx_double z,int status1, int status2)
{ // {{{ 
cx_double p11,one,two,four; 
one=cx_double(1,0); two=cx_double(2,0); four=cx_double(4,0); 
//p11=(1+(y+z)*(x-1)-sqrt(exp(ln(1+(y+z)*(x-1))*2)-4*x*(x-1)*y*z))/(2*(x-1));
p11=(one+(y+z)*(x-one)-sqrt(pow((one+(y+z)*(x-one)),two)-four*x*(x-one)*y*z));
//p11=one+(y+z)*(x-one); 
p11=p11/(two*(x-one));

// calculates probs depending on status 
if (status1==1 && status2==0) p11=y-p11;
if (status1==0 && status2==1) p11=z-p11; 
if (status1==0 && status2==0) p11=one-y-z+p11; 

return(p11); 
} // }}}

double CplacklikeP(double theta,int status1,int status2,double cif1,double cif2,vec &dp) 
{ // {{{
//double S,S2,a;
//S=1+cifs*(theta-1); S2=4*cif1*cif2*theta*(theta-1);
double x,y,z,valr=1,p11,p10,p01,p00; 
//double cifs=cif1+cif2; 
//a=(1+(theta-1)*(cifs)); 
x=theta; y=cif1; z=cif2; 

dp(0)=0; 

if (theta!=1) {
p11=(1+(y+z)*(x-1)-sqrt(pow(1+(y+z)*(x-1),2)-4*x*(x-1)*y*z))/(2*(x-1));
dp(0)= (y + z - (-4*(-1 + x)*y*z - 4*x*y*z + 2*(y + z)*(1 + (-1 + x)*(y + z)))/(2.*sqrt(-4*(-1 + x)*x*y*z + pow(1 + (  -1 + x)*(y + z),2))))/(2.*(-1 + x)) - (1 + (-1 + x)*(y + z) - sqrt(-4*(-1 + x)*x*y*z + pow(1 + (-1 + x)*(y + z  ),2)))/(2.*pow(-1 + x,2));
} else p11=cif1*cif2;

////det komplekse trick, derivative wrt y og zi, dvs D_1 P(y,z,theta) og D_2 P
//cx_double CCp11,Ctheta,Cy,Cz; 
//Ctheta=cx_double(theta,0); 
//Cy=cx_double(y,1E-20); 
//Cz=(cx_double) z; 
//CCp11=Cpij(Ctheta,Cy,Cz,status1,status2); 
//dp(1)=imag(CCp11)/1E-20; 
//Cz=cx_double(z,1E-20); 
//Cy=(cx_double) y; 
//CCp11=Cpij(Ctheta,Cy,Cz,status1,status2); 
//dp(2)=imag(CCp11)/1E-20; 
//printf(" %lf  ",imag(CCp11)/1E-20); 

p10=y-p11; 
p01=z-p11; 
p00=1-y-z+p11; 

if (status1==1 && status2==1) { valr=p11; dp(0)= dp(0); }
if (status1==1 && status2==0) { valr=p10; dp(0)=-dp(0); }
if (status1==0 && status2==1) { valr=p01; dp(0)=-dp(0); }
if (status1==0 && status2==0) { valr=p00; dp(0)= dp(0); } 
//printf(" %lf \n",dp(0)); 

return(valr); 
} // }}}
               
//double min(double a, double b) { if (a<b) return(a); else return(b); }
//double max(double a, double b) { if (a>b) return(a); else return(b); }

RcppExport SEXP twostageloglikebin( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex, SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid , SEXP irvdes,
	       	SEXP iags, SEXP ibetaiid, SEXP ipairascertained, SEXP itwostage
) // {{{
{
	try {
// {{{ setting matrices and vectors, and exporting to armadillo matrices
 mat thetades = Rcpp::as<mat>(ithetades); 
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 mat ags= Rcpp::as<colvec>(iags);
 int pt=theta.n_rows; 
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 int antclust = clusterindex.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 colvec pmargsurv = Rcpp::as<colvec>(ipmargsurv);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);
 mat rvdes= Rcpp::as<mat>(irvdes); 
 int depmodel= Rcpp::as<int>(idepmodel); 
 int twostage= Rcpp::as<int>(itwostage); 
 int pairascertained = Rcpp::as<int>(ipairascertained); 

// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);
 IntegerVector strata(istrata);

 int varlink= Rcpp::as<int>(ivarlink);
 int silent = Rcpp::as<int>(isilent);
 int iid= Rcpp::as<int>(iiid); 
 int antiid = Rcpp::as<int>(iantiid);
 double loglikecont=0; 

 mat Xtheta = Rcpp::as<mat>(iXtheta);

  int udtest=0; 
  if (udtest==1) { // {{{
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//  Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//        est.print("est"); 
//	est2.print("est2"); 
//        z.print("z"); 
//	zsem.print("zsemi"); 
//	z2.print("z2"); 
        thetades.print("theta.des"); 
        clusterindex.print("clusterindex"); 
//        rvdes.print("rvdes"); 
	theta.print("theta"); 
	Xtheta.print("Xtheta"); 
//	  y.print("y-times"); 
	  clustsize.print("clustsize"); 
	  pmargsurv.print("margsurv"); 
	  cause.print("cause"); 
	  cluster.print("cluster"); 
//	  Zgamma.print("zgam"); 
//	  Z2gamma2.print("zgam2"); 
//	  KMtimes.print("KMtimes"); 
//	  KMc.print("KMc"); 
	  weights.print("weights"); 
//	  entryage.print("entryage"); 
//	  cif1entry.print("cif1entry"); 
//	  cif2entry.print("cif2entry"); 
	  trunkp.print("trunkp"); 
  } else if (udtest==2) 
  {
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//     Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//      Rprintf("est %lf \n",mean(mean(est))); 
//      Rprintf("est2 %lf \n",mean(mean(est2))); 
//      Rprintf("z %lf \n",mean(mean(z))); 
//      Rprintf("zsem %lf \n",mean(mean(zsem))); 
//      Rprintf("z2 %lf \n",mean(mean(z2))); 
      mat mt=mean(thetades); 
      mt.print("meancol thetades"); 
//      Rprintf("theatdes %lf \n",mean(mean(thetades))); 
      Rprintf("ci %lf \n",mean(mean(clusterindex))); 
//      Rprintf("rvdes %lf \n",mean(mean(rvdes))); 
      Rprintf("theta %lf \n",mean(theta)); 
      Rprintf("Xtheta %lf \n",mean(mean(Xtheta))); 
//      Rprintf("y %lf \n",mean(y)); 
      Rprintf("ci %lf \n",mean(clustsize)); 
//      Rprintf("times %lf \n",mean(times)); 
      Rprintf("cause %lf \n",mean(cause)); 
      Rprintf("cluster %lf \n",mean(cluster)); 
//      Rprintf("Zgamma %lf \n",mean(Zgamma)); 
//      Rprintf("Z2gamma2 %lf \n",mean(Z2gamma2)); 
//      Rprintf("KMtimes %lf \n",mean(KMtimes)); 
//      Rprintf("KMc %lf \n",mean(KMc)); 
      Rprintf("weights %lf \n",mean(weights)); 
//      Rprintf("entry %lf \n",mean(entryage)); 
//      Rprintf("cif1entry %lf \n",mean(cif1entry)); 
//      Rprintf("cif2entry %lf \n",mean(cif2entry)); 
      Rprintf("trunkp %lf \n",mean(trunkp)); 
  } // }}}

  int ci,ck,i,j,c,s=0,k,v,c1; 
  double ll=1,llt=1,Li,Lk,diff=0;
  //double sdj=0;
  //  double Lit=1,Lkt=1,llt=1;
  double deppar=1,ssf=0,thetak=0; 
//  double plack(); 
  vec dplack(pt); dplack.fill(pt);
  vec dp00(pt); 
  vec ps(8); 
  vec ckij(4),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

    mat betaiid= Rcpp::as<mat>(ibetaiid);
  int dimbeta=betaiid.n_cols; 
//  printf("%d %d \n",pt,dimbeta); 
  mat DbetaDtheta(pt,dimbeta); 
  vec vDbetaDtheta(2*pt); 
  DbetaDtheta.fill(0); 
  vec Du1(pt),Du2(pt); 

  colvec p11tvec(antclust); 
//  p11tvec=0; 
//  Rprintf(" %d \n",pt); 
  int scoredim; 
  if (twostage==1) scoredim=pt; else scoredim=pt+dimbeta; 
//  printf("========= %d %d \n",pt,scoredim); 

  colvec Utheta(scoredim); 
  colvec vthetascore(scoredim); 
  colvec pthetavec(scoredim); 
  vec vtheta2(scoredim); 
  mat DUtheta(scoredim,scoredim); 
  DUtheta.fill(0); Utheta.fill(0); 

  mat thetiid(antiid,scoredim); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); loglikeiid.fill(0); }

  int nr=rvdes.n_cols; 
//  if  (depmodel==3) nr=arrayDD[2]; 
  vec rv2(nr),rv1(nr);

  vec etheta=theta; 
//  if (!Utheta.is_finite()) {  Rprintf(" NA's i def U\n"); Utheta.print("U"); }
//  if (!DUtheta.is_finite()) { Rprintf(" NA's i def DU\n"); DUtheta.print("DU"); }
  // }}}
  
//  double claytonoakesbinRVC(); 
  colvec likepairs(antclust); 

for (j=0;j<antclust;j++) if (clustsize(j)>=2) { 

    R_CheckUserInterrupt(); diff=0; //sdj=0; 

  for (c=0;c<clustsize(j)-1;c++) for (v=c+1;v<clustsize(j);v++) // {{{ //  if ((c<v)) 
  { 
     i=clusterindex(j,c); k=clusterindex(j,v); 
     if (strata(i)==strata(k)) { // {{{

     ci=cause(i); ck=cause(k); Li=pmargsurv(i); Lk=pmargsurv(k); 
         
     int flexfunc=0; 
      if (flexfunc==0) {
	  if (depmodel!=3) {
		  thetak=Xtheta(i,0);  
		  pthetavec= trans(thetades.row(i)); 
		  vthetascore=1*pthetavec; 
	  }
      } else { // {{{ 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      } // }}} 

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }
//  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{
	   ll=claytonoakesP(deppar,ci,ck,Li,Lk,dplack,ps,dp00);
	   loglikecont=log(ll); 
	   ssf+=weights(i)*log(ll); 
	   diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	if (j<-10) {
	Rprintf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	dplack.print("dtheta"); 
	}

	if (varlink==1) diff=-pow(deppar,1)*diff;  
	if (varlink==0) diff=-1*pow(deppar,2)*diff; 

          if (pairascertained==1) {
                ssf-=weights(i)*log(1-ps(0)); 
	        loglikecont=log(ll)-log(1-ps(0));
		diff+=dp00(0)/(1-ps(0));
	   if (varlink==1) diff+=-pow(deppar,1)*dp00(0)/(1-ps(0));  
	   if (varlink==0) diff+=-1*pow(deppar,2)*dp00(0)/(1-ps(0)); 
	   }

	   } // }}}
	   else if (depmodel==3) { // clayton-oakes addtive gamma  // {{{

	   rv1=trans(rvdes.row(i)); rv2=trans(rvdes.row(k)); 
	   ll=claytonoakesbinRVC(etheta,thetades,ags,ci,ck,Li,Lk,rv1,rv2,dplack,vDbetaDtheta,ps,dp00);

	   if (j<-10) {
		   rv1.print("rv1"); rv2.print("rv1"); ps.print("ps"); 
		   thetades.print("theta.des"); etheta.print("theta"); 
	           Rprintf("%d  %d %d %lf %lf %lf \n",j,ci,ck,Li,Lk,ll); 
	   }

           ssf+=weights(i)*log(ll); 
	   loglikecont=log(ll);
//	   if (varlink==1) dplackt=dplackt % etheta;  
//	   vthetascore=dplack/ll; 
//	   ps.print("ps"); 
//	   vec vvv=trans(betaiid.row(i)); vvv.print("vvv"); 
//	   vec vvv1=trans(betaiid.row(k)); vvv1.print("vvv1"); 
	   vthetascore.subvec(0,pt-1)=dplack/ll; 
	   if (twostage==0) 
	      vthetascore.subvec(pt,pt+dimbeta-1)=(ps(4)*trans(betaiid.row(i))+ps(5)*trans(betaiid.row(k)))/ll; 
	   
	if (pairascertained==1) {
           ssf-=weights(i)*log(1-ps(0)); 
	   loglikecont=log(ll)-log(1-ps(0));
           vthetascore.subvec(0,pt-1)+=dp00/(1-ps(0)); 
	   if (twostage==0) {
	      llt=(1-ps(0)); 
              vthetascore.subvec(pt,pt+dimbeta-1)+=(ps(6)*trans(betaiid.row(i))+ps(7)*trans(betaiid.row(k)))/llt; 
	   }
	}
//	   ll=claytonoakesP(deppar,ci,ck,Li,Lk,dplack);
//	   ssf+=weights(i)*log(ll); 
//	   diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
//	   if (varlink==1) diff=pow(deppar,1)*diff;  
//	   if (varlink==0) diff=1*pow(deppar,2)*diff; 
	   //sdj=-pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
           ll=placklikeP(deppar,ci,ck,Li,Lk,dplack,ps,dp00);
	   loglikecont=log(ll); 
	   ssf+=weights(i)*log(ll); 
	   // sdj=pow(dplack(0)/ll,2); 
	   if (varlink==1) diff=deppar*dplack(0)/ll; 
	   if (varlink==0) diff=dplack(0)/ll;

	   if (pairascertained==1) {
                ssf-=weights(i)*log(1-ps(0)); 
	        loglikecont-=log(1-ps(0));
	   if (varlink==1) diff+=deppar*dp00(0)/(1-ps(0)); 
	   if (varlink==0) diff+=dp00(0)/(1-ps(0));
	   }
	} // }}}


     if (depmodel!=3) {
	DUtheta+=weights(i)*pow(diff,2)*vthetascore*trans(vthetascore);
        vthetascore=weights(i)*diff*vthetascore; 
        Utheta=Utheta+vthetascore; 
     } else  { // additive gamma structure 
//		printf(" mig 1\n"); 
//		vthetascore.print("vvv"); 
	     DUtheta+=weights(i)*vthetascore*trans(vthetascore);
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 
	     if (twostage==1) {
		     Du1=vDbetaDtheta.subvec(0,pt-1); 
		     Du2=vDbetaDtheta.subvec(pt,2*pt-1); 
//	     if (j<-10) {
//	     vDbetaDtheta.print("dtds"); 
//	     Du1.print("Du1"); 
//	     Du2.print("Du1"); 
//	     }
	            DbetaDtheta+= Du1*betaiid.row(i)+ Du2*betaiid.row(k); 
	     }
//		vthetascore.print("vvv 2"); 
	}

     if (iid==1) { 
	 for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	 loglikeiid(secluster(i))+=loglikecont; 

     }
     } // }}} strata(i)==strata(k) indenfor strata

  } /* for (c=0....... */   // }}}

if (iid==1)  likepairs(j)=loglikecont; 
} /* j in antclust */ 

//printf("Sum of squares %lf \n",ssf);theta.print("theta");Utheta.print("Utheta");DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
res["DbetaDtheta"]=DbetaDtheta; 
if (iid==1) { 
             res["theta.iid"]=thetiid; 
             res["DbetaDtheta.iid"]=thetiid; 
	     res["loglikeiid"]=loglikeiid; 
             res["likepairs"]=likepairs; 
            }

return(res); 
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall



} // }}}

RcppExport SEXP twostageloglikebinpairs( 
		SEXP icause, SEXP ipmargsurv, 
		SEXP itheta, SEXP iXtheta, SEXP iDXtheta, SEXP idimDX, SEXP ithetades,
		SEXP icluster,SEXP iclustsize,SEXP iclusterindex, SEXP ivarlink, 
                SEXP iiid, SEXP  iweights, SEXP isilent, 
		SEXP idepmodel, // SEXP ientryage,
		SEXP itrunkp , SEXP istrata, SEXP isecluster, SEXP  iantiid , SEXP irvdes,
		SEXP idimthetades, SEXP idimrvdes, SEXP inrvs, SEXP iags, 
		SEXP ibetaiid, SEXP ipairascertained, SEXP itwostage,SEXP icasecontrol 
) // {{{
{
  try {
// {{{ setting matrices and vectors, and exporting to armadillo matrices
// mat thetades = Rcpp::as<mat>(ithetades); 
 IntegerVector nrvs(inrvs);
 mat clusterindex = Rcpp::as<mat>(iclusterindex);
 colvec theta = Rcpp::as<colvec>(itheta);
 mat ags= Rcpp::as<colvec>(iags);
 int pt=theta.n_rows; 
 colvec clustsize = Rcpp::as<colvec>(iclustsize);
 // this is number of pairs (rather than clusters)
 int antclust = clusterindex.n_rows; 
 colvec cause = Rcpp::as<colvec>(icause);
 colvec pmargsurv = Rcpp::as<colvec>(ipmargsurv);
 colvec cluster = Rcpp::as<colvec>(icluster);
 colvec weights = Rcpp::as<colvec>(iweights);
// colvec entryage = Rcpp::as<colvec>(ientryage);
 colvec trunkp = Rcpp::as<colvec>(itrunkp);
 colvec secluster = Rcpp::as<colvec>(isecluster);
// mat rvdes= Rcpp::as<mat>(irvdes); 
 int depmodel= Rcpp::as<int>(idepmodel); 
 IntegerVector strata(istrata);

//    uvec                 cause = Rcpp::as<uvec>(icause); 
 uvec pairascertained = Rcpp::as<uvec>(ipairascertained); 
 uvec casecontrol     = Rcpp::as<uvec>(icasecontrol); 

// array for derivative of flexible design
 NumericVector DXthetavec(iDXtheta);
 IntegerVector arrayDims(idimDX);
 arma::cube DXtheta(DXthetavec.begin(), arrayDims[0], arrayDims[1], arrayDims[2], false);

// mat thetades = Rcpp::as<mat>(ithetades); 
// array for parameter restrictions (one for each pair) pairs * (ant random effects)* (ant par)
 NumericVector thetadesvec(ithetades);
 IntegerVector arrayDims1(idimthetades); 
 IntegerVector arrayDD(3); 
//printf(" mig thetades 222222\n"); 
 if (depmodel==3)  {
	 arrayDD[0]=arrayDims1[0]; arrayDD[1]=arrayDims1[1]; arrayDD[2]=arrayDims1[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
arma::cube thetadesi(thetadesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);
//printf(" mig cube thetades 222222\n"); 

mat thetades=mat(thetadesvec.begin(),arrayDims1[0],arrayDims1[1]*arrayDD[2],false); 

// array for pairwise random effects (two vectors for each pair)  pairs * 2* (ant random effects)
 NumericVector rvdesvec(irvdes);
 IntegerVector arrayDims2(idimrvdes); 
 if (depmodel==3)  {
	 arrayDD[0]=arrayDims2[0]; arrayDD[1]=arrayDims2[1]; arrayDD[2]=arrayDims2[2];  }
 else { arrayDD[0]=1; arrayDD[1]=1; arrayDD[2]=1;  }
// printf(" %d %d %d \n",arrayDD[0], arrayDD[1], arrayDD[2]); 
arma::cube rvdesC(rvdesvec.begin(), arrayDD[0], arrayDD[1], arrayDD[2], false);

mat rvdes=mat(rvdesvec.begin(),arrayDims2[0],arrayDims2[1]*arrayDD[2],false); 
// mat rvdes= Rcpp::as<mat>(irvdes); 

 int varlink= Rcpp::as<int>(ivarlink);
 int silent = Rcpp::as<int>(isilent);
 int iid= Rcpp::as<int>(iiid); 
 int antiid = Rcpp::as<int>(iantiid);
 double loglikecont=0; 
 mat Xtheta = Rcpp::as<mat>(iXtheta);
 int twostage= Rcpp::as<int>(itwostage); 

 int udtest=0; 
  if (udtest==1) { // {{{
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//  Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//        est.print("est"); 
//	est2.print("est2"); 
//        z.print("z"); 
//	zsem.print("zsemi"); 
//	z2.print("z2"); 
        thetades.print("theta.des"); 
        clusterindex.print("clusterindex"); 
//        rvdes.print("rvdes"); 
	theta.print("theta"); 
	Xtheta.print("Xtheta"); 
//	  y.print("y-times"); 
	  clustsize.print("clustsize"); 
	  pmargsurv.print("margsurv"); 
	  cause.print("cause"); 
	  cluster.print("cluster"); 
//	  Zgamma.print("zgam"); 
//	  Z2gamma2.print("zgam2"); 
//	  KMtimes.print("KMtimes"); 
//	  KMc.print("KMc"); 
	  weights.print("weights"); 
//	  entryage.print("entryage"); 
//	  cif1entry.print("cif1entry"); 
//	  cif2entry.print("cif2entry"); 
	  trunkp.print("trunkp"); 
  } else if (udtest==2) 
  {
//  Rprintf(" %d %d %d %d %d %d %d \n",samecens,inverse,semi,semi2,flexfunc,stabcens,silent); 
//     Rprintf(" %d %d %d %d %d %d %d \n",cifmodel,CA1,CA2,sym,depmodel,estimator,iid); 
//      Rprintf("est %lf \n",mean(mean(est))); 
//      Rprintf("est2 %lf \n",mean(mean(est2))); 
//      Rprintf("z %lf \n",mean(mean(z))); 
//      Rprintf("zsem %lf \n",mean(mean(zsem))); 
//      Rprintf("z2 %lf \n",mean(mean(z2))); 
      mat mt=mean(thetades); 
      mt.print("meancol thetades"); 
//      Rprintf("theatdes %lf \n",mean(mean(thetades))); 
      Rprintf("ci %lf \n",mean(mean(clusterindex))); 
//      Rprintf("rvdes %lf \n",mean(mean(rvdes))); 
      Rprintf("theta %lf \n",mean(theta)); 
      Rprintf("Xtheta %lf \n",mean(mean(Xtheta))); 
//      Rprintf("y %lf \n",mean(y)); 
      Rprintf("ci %lf \n",mean(clustsize)); 
//      Rprintf("times %lf \n",mean(times)); 
      Rprintf("cause %lf \n",mean(cause)); 
      Rprintf("cluster %lf \n",mean(cluster)); 
//      Rprintf("Zgamma %lf \n",mean(Zgamma)); 
//      Rprintf("Z2gamma2 %lf \n",mean(Z2gamma2)); 
//      Rprintf("KMtimes %lf \n",mean(KMtimes)); 
//      Rprintf("KMc %lf \n",mean(KMc)); 
      Rprintf("weights %lf \n",mean(weights)); 
//      Rprintf("entry %lf \n",mean(entryage)); 
//      Rprintf("cif1entry %lf \n",mean(cif1entry)); 
//      Rprintf("cif2entry %lf \n",mean(cif2entry)); 
      Rprintf("trunkp %lf \n",mean(trunkp)); 
  } // }}}

  int sss=1,ci,ck,i,j,s=0,k,c1; 
  double llt=1,ll=1,Li,Lk,diff=0;
  double deppar=1,ssf=0,thetak=0; 
  vec dplack(pt); dplack.fill(0);
  vec dp00(pt); dp00.fill(0);
  vec ps(8); ps.fill(0);
  vec ckij(4),dckij(4),ckijvv(4),dckijvv(4),ckijtv(4),dckijtv(4),ckijvt(4),dckijvt(4);
  i=silent+1; 

  colvec p11tvec(antclust); 
//  p11tvec=0; 
//  Rprintf(" %d \n",pt); 
  int scoredim; 
  mat betaiid= Rcpp::as<mat>(ibetaiid);
  int dimbeta=betaiid.n_cols; 
  if (twostage==1) scoredim=pt; else scoredim=pt+dimbeta; 
  colvec Utheta(scoredim); 
  colvec vthetascore(scoredim); 
  colvec pthetavec(scoredim); 
  vec vtheta2(scoredim); 
  mat DUtheta(scoredim,scoredim); 
  DUtheta.fill(0); Utheta.fill(0); 

  mat thetiid(antiid,scoredim); 
  colvec loglikeiid(antiid); 
  if (iid==1) { thetiid.fill(0); loglikeiid.fill(0); }

  int nr=1; 
  if  (depmodel==3) nr=arrayDD[2]; 
  vec rv2(nr),rv1(nr);

  mat DbetaDtheta(pt,dimbeta); 
  vec vDbetaDtheta(2*pt); 
  DbetaDtheta.fill(0); 
  vec Du1(pt),Du2(pt); 

  vec etheta=theta; 
//  if (!Utheta.is_finite()) {  Rprintf(" NA's i def U\n"); Utheta.print("U"); }
//  if (!DUtheta.is_finite()) { Rprintf(" NA's i def DU\n"); DUtheta.print("DU"); }
  // }}}
  
  colvec likepairs(antclust); 
//  double claytonoakesbinRVC(); 

//    Rprintf("--------%d  \n",antclust); 

for (j=0;j<antclust;j++) { 

    R_CheckUserInterrupt(); diff=0; //sdj=0; 

//    printf("--------%d %d %d %lf %lf \n",j,ci,ck,Li,Lk); 
// index of subject's in pair "j"
   i=clusterindex(j,0); k=clusterindex(j,1); 
//    printf("--------%d %d %d %lf %lf \n",j,i,k,Li,Lk); 

    if (strata(i)==strata(k)) { // {{{

     ci=cause(i); ck=cause(k); Li=pmargsurv(i); Lk=pmargsurv(k); 
//     printf("%d %d %d %d %lf %lf \n",j,(int) secluster(i),ci,ck,Li,Lk); 
         
     int flexfunc=0; 
      if (flexfunc==0) {
	  if (depmodel!=3) {
//		  thetak=Xtheta(i,0);  
//		  pthetavec= trans(thetades.row(i)); 
//		  vthetascore=1*pthetavec; 
                  thetak=Xtheta(j,0);  
		  pthetavec= trans(thetades.row(j)); 
		  vthetascore=1*pthetavec; 
	  }
      } else { // {{{ 
	  thetak=Xtheta(i,s); 
	  pthetavec = DXtheta(span(s),span(i),span::all); 
      } // }}} 

//  pthetavec.print("pthetavec");
//  printf(" %lf \n",thetak);

  if (depmodel==1){ if (varlink==1) deppar=1/exp(thetak); else deppar=1/thetak;}
  if (depmodel==2){ if (varlink==1) deppar=exp(thetak); else deppar=thetak; }
//  if (depmodel==3){ if (varlink==1) etheta=exp(theta); else etheta=theta; }

	if (depmodel==1) { // clayton-oakes  // {{{
	   ll=claytonoakesP(deppar,ci,ck,Li,Lk,dplack,ps,dp00);
	   loglikecont=log(ll); 
	   ssf+=weights(i)*log(ll); 
	   diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   
	   if (pairascertained[j]==1) {
                ssf-=weights(i)*log(1-ps(0)); 
	        loglikecont=log(ll)-log(1-ps(0));
		diff+=dp00(0)/(1-ps(0));
	   if (varlink==1) diff=-pow(deppar,1)*diff;  
	   if (varlink==0) diff=-1*pow(deppar,2)*diff; 
	   }

	   } // }}}
	   else if (depmodel==3) { // clayton-oakes addtive gamma  // {{{

	// 3-dimensional array pairs*(2xrandom effects)
        int lnrv= nrvs(j)-1; // number of random effects for this cluster 	
//	Rprintf(" %d \n",lnrv); 
	mat rv=rvdesC.slice(j); 
        vec rv1= trans(rv.submat(span(0),span(0,lnrv)));
        vec rv2= trans(rv.submat(span(1),span(0,lnrv)));

//        rv1= rvdesC.subcube( span(0),span(0,lnrv),span(j));
//        rv2= rvdesC.subcube( span(1),span(0,lnrv),span(j));

	// takes parameter relations for each pair
	// 3-dimensional array pairs*(random effects* pars )
//	mat thetadesv=thetadesi.slice(j); 
	mat thetadesvv=thetadesi.slice(j); 
	mat thetadesv=thetadesvv.rows(0,lnrv); 

	if (j< -10)  {/*{{{*/
	   Rprintf(" %d %d \n",lnrv,pt); 
           rv1.print("rv1");    rv2.print("rv2"); 
	   thetadesvv.print("thetades "); 
	   etheta.print("e-theta"); 
//	   mat test=mat(thetades.begin(),3,1); 
//	   test.print("test"); 
	}/*}}}*/

	ll=claytonoakesbinRVC(etheta,thetadesv,ags,ci,ck,Li,Lk,rv1,rv2,dplack,vDbetaDtheta,ps,dp00);
        if (j<-10) { 
		Rprintf("==============================\n"); 
		rv1.print("rv1"); rv2.print("rv1"); thetadesv.print("theta.des"); etheta.print("theta"); 
		Rprintf("%d  %d %d %lf %lf %lf \n",j,ci,ck,Li,Lk,ll); ps.print("ps"); 
	}

        ssf+=weights(i)*log(ll); 
	loglikecont=log(ll);
	vthetascore.subvec(0,pt-1)=dplack/ll; 
	if (twostage==0) 
	   vthetascore.subvec(pt,pt+dimbeta-1)=(ps(4)*trans(betaiid.row(i))+ps(5)*trans(betaiid.row(k)))/ll; 

	if (pairascertained(j)==1) {
           ssf-=weights(i)*log(1-ps(0)); 
	   loglikecont=log(ll)-log(1-ps(0));
	   vthetascore.subvec(0,pt-1)+=dp00/(1-ps(0)); 
	   if (twostage==0) {
	      llt=(1-ps(0)); 
              vthetascore.subvec(pt,pt+dimbeta-1)+=(ps(6)*trans(betaiid.row(i))+ps(7)*trans(betaiid.row(k)))/llt; 
	   }
	}

	if (casecontrol(j)==1) { //         marginal of second person Lk
	   llt=(ck==1)*(Lk)+(ck==0)*(1-Lk); 
           ssf-=weights(i)*log(llt); 
	   loglikecont=log(ll)-log(llt);
	   if (twostage==0) {
	      sss=-1; 
	      if (ck==0) sss=1; 
              vthetascore.subvec(pt,pt+dimbeta-1)=vthetascore.subvec(pt,pt+dimbeta-1)+sss*trans(betaiid.row(k))/llt; 
	   }
	}

//	   ll=claytonoakesP(deppar,ci,ck,Li,Lk,dplack);
//	   ssf+=weights(i)*log(ll); 
//	   diff=dplack(0)/ll; 
//	printf(" %d %d %d %d %d  \n",j,c,v,i,k); 
//	printf(" %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf \n",j,c,v,i,k,ci,ck,thetak,Li,Lk,weights(i),ll,log(ll)); 
//	   if (varlink==1) diff=pow(deppar,1)*diff;  
//	   if (varlink==0) diff=1*pow(deppar,2)*diff; 
	   //sdj=-pow(diff,2); 
	   // }}}
	} else if (depmodel==2) { // plackett model  // {{{
           ll=placklikeP(deppar,ci,ck,Li,Lk,dplack,ps,dp00);
	   loglikecont=log(ll); 
	   ssf+=weights(i)*log(ll); 
	   // sdj=pow(dplack(0)/ll,2); 
	   
	   if (varlink==1) diff=deppar*dplack(0)/ll; 
	   if (varlink==0) diff=dplack(0)/ll;

	   if (pairascertained(j)==1) {
                ssf-=weights(i)*log(1-ps(0)); 
	        loglikecont-=log(1-ps(0));
	   if (varlink==1) diff+=deppar*dp00(0)/(1-ps(0)); 
	   if (varlink==0) diff+=dp00(0)/(1-ps(0));
	   }
	   
	} // }}}


     if (depmodel!=3) {
	     DUtheta+=weights(i)*pow(diff,2)*vthetascore*trans(vthetascore);
	     vthetascore=weights(i)*diff*vthetascore; 
	     Utheta=Utheta+vthetascore; 
     } else  { // additive gamma structure 
//		printf(" mig 1\n"); 
//		vthetascore.print("vvv"); 
	     DUtheta+=weights(i)*vthetascore*trans(vthetascore);
//		printf(" mig 1\n"); 
	     vthetascore=weights(i)*vthetascore; 
	     Utheta=Utheta+vthetascore; 

	     if (twostage==1) {
		     Du1=vDbetaDtheta.subvec(0,pt-1); 
		     Du2=vDbetaDtheta.subvec(pt,2*pt-1); 
		     DbetaDtheta+= Du1*betaiid.row(i)+ Du2*betaiid.row(k); 
	     }
	}

     if (iid==1) { 
	 for (c1=0;c1<pt;c1++) thetiid((int) secluster(i),c1)+=vthetascore(c1); 
	 loglikeiid((int) secluster(i))+=loglikecont; 
     }
     } // }}} strata(i)==strata(k) indenfor strata

//   printf(" hej 2 \n"); 

if (iid==1)  likepairs(j)=ll; 
} /* j in antclust */ 

//printf("Sum of squares %lf \n",ssf);theta.print("theta");Utheta.print("Utheta");DUtheta.print("DUtheta"); 
List res; 
res["loglike"]=ssf; 
res["score"]=Utheta; 
res["Dscore"]=DUtheta; 
res["DbetaDtheta"]=DbetaDtheta; 
if (iid==1) { res["theta.iid"]=thetiid; 
	      res["loglikeiid"]=loglikeiid; 
              res["likepairs"]=likepairs; 
            }

return(res); 
} catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {  
    ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall

} // }}}

