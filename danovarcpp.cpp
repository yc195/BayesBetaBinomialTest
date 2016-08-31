#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//function for a_i(x) = \sum_j \digamma(x+m_ij) - r_i\digamma(x)
// [[Rcpp::export]]
double af(double x, arma::vec mi) {
  int K=mi.size();
  double a=0;
  
  for (int j=0;j<K;j++) {
    if (x<0.00001) {
      a=a+100000;
    } else {
      a=a+1/x;
    }
    for (int l=1;l<mi(j);l++) {
      a=a+1/(x+l);
    }
  }
  return a;
}

//function for b_i(x) = \sum_j \trigamma(x+m_ij) - r_i\trigamma(x)
// [[Rcpp::export]]
double bf(double x, arma::vec mi) {
  int K=mi.size();
  double a=0;
  
  for (int j=0;j<K;j++) {
    if (x<0.00001) {
      a=a-10000000000;
    } else {
      a=a-1/(x*x);
    }
    for (int l=1;l<mi(j);l++) {
      a=a-1/pow(x+l,2);
    }
  }
  return a;
}


//finds the loglikelihood of s \in \{0,1\}^p with fixed nu and with priors
// [[Rcpp::export]]
double calcll(arma::mat y, arma::mat fullx, arma::vec s, int nrun, arma::mat ntest, double nu, double sigma2) {
  int pf=fullx.n_cols;
  int p=sum(s)+1;
  int n=y.n_rows;
  int K=y.n_cols;
  
  arma::mat x=arma::zeros(n,p);
  
  double count=1;
  x.col(0) = fullx.col(0);
  for(int j=0; j<pf-1; j++) {
    if(s(j)==1) {
      x.col(count) = fullx.col(j+1);
      count++;
    }
  }
  
  arma::mat betadraw=arma::zeros(p,nrun);
  arma::mat W1=arma::zeros(n,n);
  arma::mat W2=arma::zeros(n,n);
  arma::mat diags=arma::eye(p,p);
  arma::vec theta(n),  z(n);
  
  for (int t=1; t<nrun; t++) {
    theta=1/(1+exp(-x*betadraw.col(t-1)));

    for (int i=0; i<n; i++) {
      z(i)=theta(i)*(1-theta(i));
      W1(i, i)=af(theta(i)*nu, (y.row(i)).t())-af((1-theta(i))*nu, (ntest.row(i)).t()-(y.row(i)).t());
      W2(i, i)=-nu*(bf(theta(i)*nu, (y.row(i)).t())+bf((1-theta(i))*nu, (ntest.row(i)).t()-(y.row(i)).t()))*pow(z(i), 2);
      W2(i, i)=W2(i,i)-W1(i,i)*z(i)*(1-2*theta(i));
    }
    betadraw.col(t)=betadraw.col(t-1)+(x.t()*W2*x+diags/(nu*sigma2)).i()*(x.t()*W1*z-betadraw.col(t-1)/(nu*sigma2));
  }
  double detHess=det(nu*x.t()*W2*x+diags/sigma2);
  double loglh=0;
  for (int i=0;i<n;i++) 
  {
    for (int j=0; j<K; j++) 
    {
      loglh=loglh+lgamma(theta(i)*nu+y(i, j))+lgamma((1-theta(i))*nu+ntest(i, j)-y(i, j))-lgamma(nu+ntest(i,j))-lgamma(theta(i)*nu)-lgamma((1-theta(i))*nu)+lgamma(nu);
    }
  }
  
  double sumdnorm=-0.5*p*log(sigma2);
  for (int u=0;u<p;u++) 
  {
    sumdnorm=sumdnorm-0.5*pow(betadraw(u, nrun-1), 2)/sigma2;
  }

  double loglaplace=loglh+sumdnorm-0.5*log(abs(detHess));

  return loglaplace;
}

