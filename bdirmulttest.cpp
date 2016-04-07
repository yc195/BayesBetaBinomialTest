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


//define a function that maps a K-1 dimensional Euclidean space to a K-1 dimensional probability simplex
// [[Rcpp::export]]
NumericVector simplexMap(NumericVector u)
{
  int K=u.size()+1;
  NumericVector p(K);
  p(K-1)=(1-h(u(0)))/2;
  if (K==2) {
    p(0)=(1+h(u(0)))/2;
  } else if (K==3) {
    p(0)=(1+h(u(0)))*(1-h(u(1)))/4;
    p(1)=(1+h(u(0)))*(1+h(u(1)))/4;
  } else {
    p(0)=(1+h(u(0)))*(1-h(u(1)))/4;
    for (int i=1;i<K-2;i++) {
      p(i)=p(i-1)*(1+h(u(i)))*(1-h(u(i+1)))/(2*(1-h(u(i))));
    }
    p(K-2)=p(K-3)*(1+h(u(K-2)))/(1-h(u(K-2)));
  }
  return p;
}


//find the absolute value of determinant of the Jacobian matrix of mapping in simplexMap
// [[Rcpp::export]]
double calcJacobian(NumericVector u)
{
  int K=u.size()+1;
  NumericVector p(K);
  p=simplexMap(u);
  double absdetJ;
  NumericMatrix J(K-1,K-1);
  if (K==2)
  {
    absdetJ=hprime(u(0))/2;
  } else {
    for (int i=0; i<K-1;i++) {
      for (int j=0; j<K-1;j++) {
        if (j<=i) {
          J(i,j)=p(i)*hprime(u(j))/(1+h(u(j)));
        } else if (j==i+1) {
          J(i,j)=-p(i)*hprime(u(j))/(1-h(u(j)));
        } else {
          J(i,j)=0;
        }
      }
    }
  absdetJ=fabs(determ(J));
  }
  return absdetJ;
}


//find the integrand of the integral under H1 as function of vector u and eta 
// [[Rcpp::export]]
double funcH1(NumericMatrix n1, NumericMatrix n2, NumericVector uetavec1, double a1, double b1)
{
  int S1=n1(_,0).size();
  int S2=n2(_,0).size();
  int K=n1(0,_).size();
  NumericVector p1(K), p2(K);
  NumericVector u1(K-1), u2(K-1);
  
  for (int h=0;h<K-1;h++) {
   u1(h)=uetavec1(h);
   u2(h)=uetavec1(h+K-1);
  }
  double eta1=exp(uetavec1(2*K-2));
  double eta2=exp(uetavec1(2*K-1));
  
  p1=simplexMap(u1);
  p2=simplexMap(u2);
  
  double B=0;
  for (int j=0;j<S1;j++) {
    B=B+mlbeta(n1(j,_)+eta1*p1)-mlbeta(eta1*p1);
  }
  for (int j=0;j<S2;j++) {
    B=B+mlbeta(n2(j,_)+eta2*p2)-mlbeta(eta2*p2);
  }
  B=B+(1/K-1)*sum(log(p1)+log(p2))+log(calcJacobian(u1))+log(calcJacobian(u2));
  B=B+a1*(log(eta1)+log(eta2))-b1*(eta1+eta2);
  return (-B);
}


//find the integrand of the integral under H0 as function of vector u and eta 
// [[Rcpp::export]]
double funcH0(NumericMatrix n1, NumericMatrix n2, NumericVector uetavec0, double a0, double b0)
{
  int S1=n1(_,0).size();
  int S2=n2(_,0).size();
  int K=n1(0,_).size();
  NumericVector p(K), u0(K-1);
  for (int i=0;i<K-1;i++) {
    u0(i)=uetavec0(i);
  }
  double eta1=exp(uetavec0(K-1));
  double eta2=exp(uetavec0(K));
  
  p=simplexMap(u0);
  
  double B=0;
  for (int j=0;j<S1;j++) {
    B=B+mlbeta(n1(j,_)+eta1*p)-mlbeta(eta1*p);
  }
  for (int j=0;j<S2;j++) {
    B=B+mlbeta(n2(j,_)+eta2*p)-mlbeta(eta2*p);
  }
  B=B+(1/K-1)*sum(log(p))+log(calcJacobian(u0));
   B=B+a0*(log(eta1)+log(eta2))-b0*(eta1+eta2);
  return (-B);
}