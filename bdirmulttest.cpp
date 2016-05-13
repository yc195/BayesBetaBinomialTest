#include <math.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;


//multvariate log beta function
// [[Rcpp::export]]
double mlbeta(NumericVector x) {
  int K=x.size();
  double nx = 0;
  double nx2 = 0;
  for(int k = 0; k < K; k++) {
    nx = nx + lgamma(x(k));
    nx2 = nx2 + x(k);
  }
  return (nx-lgamma(nx2));
}

//find the determinant of a matrix 
// [[Rcpp::export]]
double determ( NumericMatrix x) {
  arma::mat a;
  a=Rcpp::as<arma::mat>(x);
  double b=det(a);
  return b;
}

//define a bijective function h that maps R to (-1,1)
// [[Rcpp::export]]
double h(double x) 
{
  return x*pow(1+x*x, -0.5);
}

//find the derivative of function h
// [[Rcpp::export]]
double hprime(double x) 
{
  return pow(1+x*x, -1.5);
}


//find loglikelihood under H1 with fixed eta and baseline logit transform and with priors
// [[Rcpp::export]]
double logitnoeta1(NumericMatrix n1, NumericMatrix n2, NumericVector uetavec1, double eta1, double eta2)
{
  int S1=n1(_,0).size();
  int S2=n2(_,0).size();
  int K=n1(0,_).size();
  NumericVector p1(K), p2(K);
  NumericVector u1(K-1), u2(K-1);
  double sum_expu1=1, sum_expu2=1;
  
  for (int h=0;h<K-1;h++) {
    u1(h)=uetavec1(h);
    u2(h)=uetavec1(h+K-1);
    sum_expu1=sum_expu1+exp(u1(h));
    sum_expu2=sum_expu2+exp(u2(h));
  }
  
  for (int h=0;h<K-1;h++)
  {
    p1(h)=exp(u1(h))/sum_expu1;
    p2(h)=exp(u2(h))/sum_expu2;
  }
  p1(K-1)=1/sum_expu1;
  p2(K-1)=1/sum_expu2;
  
  double B=0;
  for (int j=0;j<S1;j++) {
    B=B+mlbeta(n1(j,_)+eta1*p1)-mlbeta(eta1*p1);
  }
  for (int j=0;j<S2;j++) {
    B=B+mlbeta(n2(j,_)+eta2*p2)-mlbeta(eta2*p2);
  }
  B=B+sum(log(p1)+log(p2))/K;
  return (-B);
}


//find loglikelihood under H0 with fixed eta and baseline logit transform and with priors
// [[Rcpp::export]]
double logitnoeta0(NumericMatrix n1, NumericMatrix n2, NumericVector uetavec0, double eta1, double eta2)
{
  int S1=n1(_,0).size();
  int S2=n2(_,0).size();
  int K=n1(0,_).size();
  double sum_expu=1;
  
  NumericVector p(K), u0(K-1);
  for (int h=0;h<K-1;h++) {
    u0(h)=uetavec0(h);
    sum_expu=sum_expu+exp(u0(h));
  }
  
  for (int h=0;h<K-1;h++)
  {
    p(h)=exp(u0(h))/sum_expu;
  }
  p(K-1)=1/sum_expu;
  
  double B=0;
  for (int j=0;j<S1;j++) {
    B=B+mlbeta(n1(j,_)+eta1*p)-mlbeta(eta1*p);
  }
  for (int j=0;j<S2;j++) {
    B=B+mlbeta(n2(j,_)+eta2*p)-mlbeta(eta2*p);
  }
  B=B+sum(log(p))/K;
  return (-B);
}
