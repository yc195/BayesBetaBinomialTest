#include <RcppArmadillo.h>
#include <math.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;



// [[Rcpp::export]]
double aij(double x, double y) {   //digamma(x+y)-digamma(x)
  double a=0;
  if (x<0.00001) {
    a=100000;
  } else {
    a=1/x;
  }
  for (int l=1;l<y;l++) {
    a=a+1/(x+l);
  }
  return a;
}


// [[Rcpp::export]]
double bij(double x, double y) {    //trigamma( function of (x+y))x+y)-trigamma(x)
  double a=0;
  if (x<0.00001) {
    a=10000000000;
  } else {
    a=-1/(x*x);
  }
  for (int l=1;l<y;l++) {
    a=a-1/((x+l)*(x+1));
  }
  return a;
}


//calculates the log determinant
// [[Rcpp::export]]
double logdethess(arma::mat H) {  
  arma::mat G=chol(H);
  arma::vec d=diagvec(G);
  double detH=2*sum(log(d));
  return detH;
}


//finds the log likelihood under H1
// [[Rcpp::export]]
List calcloglhH1(arma::vec y, arma::mat fullx, arma::vec s, arma::vec r, int nrun, arma::vec ntest, arma::vec nu, double sigma2beta, double sigma2gamma) {
  int pf=fullx.n_cols;
  int m=sum(r);
  int n=r.size();
  int p=sum(s)+1;   
  
  arma::mat xidraw=arma::zeros(p+n-1,nrun);
  arma::vec beta=arma::zeros(p);
  arma::vec gamma=arma::zeros(n-1);
  arma::vec V=arma::zeros(n+p-1);
  arma::mat W=arma::zeros(n+p-1,n+p-1);  
  double zij;
  
  arma::mat x=arma::zeros(m, p);
  
  double count=1;
  x.col(0) = fullx.col(0);
  for(int j=0; j<pf-1; j++) {
    if(s(j)==1) {
      x.col(count) = fullx.col(j+1);
      count++;
    }
  }
  
  //(i, j)<-->mi(i)+j
  arma::vec mi(n);  //mi(i)=0+r(0)+...+r(i-1)
  mi(0)=0;
  for (int i=0;i<n-1;i++) {
    mi(i+1)=mi(i)+r(i);
  }
  
  arma::mat D=arma::zeros(n+p-1,n+p-1);  
  
  for (int q=0;q<p;q++) {     //diagonal matrix D
    D(q,q)=1/sigma2beta;
  }
  for (int i=0;i<n-1;i++) {
    D(p+i,p+i)=1/sigma2gamma;
  }
  
  double eta;
  double theta;
  for (int t=1; t<nrun; t++) {
    V=arma::zeros(n+p-1);
    W=arma::zeros(n+p-1,n+p-1);
    
    for (int i=0; i<n; i++) {
      arma::mat Ai=arma::zeros(r(i),r(i));
      arma::mat Ci=arma::zeros(r(i),r(i));  //Ci=Bi+Ai diga(1-2theta)
      arma::mat xEi(r(i), n+p-1);      //matrix [x Ei']
      arma::vec Ai1=arma::zeros(r(i));    //A_i 1_r(i)
      
      for (int j=0;j<r(i);j++){
        if (i==0) {
          eta=0;
        } else {
          eta=gamma(i-1);
        }
        
        for (int q=0;q<p;q++) {
          eta=eta+x(mi(i)+j, q)*beta(q);
        }
        theta=1/(1+exp(-eta));
        zij=nu(i)*theta*(1-theta);
        
        Ai(j, j)=(aij(theta*nu(i), y(mi(i)+j))-aij((1-theta)*nu(i), ntest(mi(i)+j)-y(mi(i)+j)))*zij;
        Ci(j, j)=(bij(theta*nu(i), y(mi(i)+j))+bij((1-theta)*nu(i), ntest(mi(i)+j)-y(mi(i)+j)))*zij*zij+Ai(j, j)*(1-2*theta);
        Ai1(j)=Ai(j, j);
        
        xEi=arma::zeros(r(i),n+p-1);
        for (int j=0;j<r(i);j++) {
          for (int q=0; q<p;q++) {
            xEi(j, q)=x(mi(i)+j, q);
          }
          if (i!=0) {
            xEi(j,p+i-1)=1;
          }
        }
      }
      V=V+xEi.t()*Ai1;
      W=W+xEi.t()*Ci*xEi;
    }
    xidraw.col(t)=xidraw.col(t-1)-(W-D).i()*(V-D*xidraw.col(t-1));
    
    for (int q=0;q<p;q++) {
      beta(q)=xidraw(q, t);
    }
    for (int i=0;i<n-1;i++) {
      gamma(i)=xidraw(p+i, t);
    }
  }
  
  double logdetHess=logdethess(D-W);
  double loglhH1=0;

  for (int i=0;i<n;i++) {
    for (int j=0; j<r(i); j++) {
      if(i==0) {
        eta=0;
      } else {
        eta=gamma(i-1);
      }
      
      for (int q=0;q<p;q++) {
        eta=eta+x(mi(i)+j, q)*beta(q);
      }
      theta=1/(1+exp(-eta));
      loglhH1=loglhH1+lgamma(theta*nu(i)+y(mi(i)+j))+lgamma((1-theta)*nu(i)+ntest(mi(i)+j)-y(mi(i)+j))
              -lgamma(nu(i)+ntest(mi(i)+j))-lgamma(theta*nu(i))-lgamma((1-theta)*nu(i))+lgamma(nu(i));
    }
  }
  double sumdnormbeta=-0.5*p*log(sigma2beta);
  for (int u=0;u<p;u++) {
    sumdnormbeta=sumdnormbeta-0.5*pow(beta(u), 2)/sigma2beta;
  }
  double sumdnormgamma=-0.5*(n-1)*log(sigma2gamma);
  for (int u=0;u<n-1;u++) {
    sumdnormgamma=sumdnormgamma-0.5*pow(gamma(u), 2)/sigma2gamma;
  }
  double loglaplaceH1=loglhH1+sumdnormbeta+sumdnormgamma-0.5*logdetHess;
  
  return Rcpp::List::create(
    Rcpp::Named("xidraw") = xidraw,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("logdetHess") = logdetHess,    
    Rcpp::Named("loglaplaceH1") = loglaplaceH1
    
  );
}


//finds the log likelihood under H0
// [[Rcpp::export]]
List calcloglhH0(arma::vec y, arma::mat fullx, arma::vec s, arma::vec r, int nrun, arma::vec ntest, arma::vec nu, double sigma2beta) {
  int pf=fullx.n_cols;
  int m=sum(r);     //m=fullx.n_rows
  int n=r.size();
  int p=sum(s)+1;   
  
  arma::mat betadraw=arma::zeros(p,nrun);
  arma::vec V=arma::zeros(p);
  arma::mat W=arma::zeros(p,p);  
  double zij;
  
  arma::mat x=arma::zeros(m, p);
  
  double count=1;
  x.col(0) = fullx.col(0);
  for(int j=0; j<pf-1; j++) {
    if(s(j)==1) {
      x.col(count) = fullx.col(j+1);
      count++;
    }
  }
  
  //(i, j)<-->mi(i)+j
  arma::vec mi(n);  //mi(i)=0+r(0)+...+r(i-1)
  mi(0)=0;
  for (int i=0;i<n-1;i++) {
    mi(i+1)=mi(i)+r(i);
  }
  
  arma::mat D=arma::eye(p,p)/sigma2beta;  

  double eta;
  double theta;
  for (int t=1; t<nrun; t++) {
    V=arma::zeros(p);
    W=arma::zeros(p,p);
    
    for (int i=0; i<n; i++) {
      arma::mat Ai=arma::zeros(r(i),r(i));
      arma::mat Ci=arma::zeros(r(i),r(i));  //Ci=Bi+Ai diga(1-2theta)
      arma::mat xi(r(i), p);      //matrix [x Ei']
      arma::vec Ai1=arma::zeros(r(i));    //A_i 1_r(i)
      
      for (int j=0;j<r(i);j++) {
        eta=0;
        for (int q=0;q<p;q++) {
          eta=eta+x(mi(i)+j, q)*betadraw(q, t-1);
        }
        theta=1/(1+exp(-eta));
        zij=nu(i)*theta*(1-theta);
        
        Ai(j, j)=(aij(theta*nu(i), y(mi(i)+j))-aij((1-theta)*nu(i), ntest(mi(i)+j)-y(mi(i)+j)))*zij;
        Ci(j, j)=(bij(theta*nu(i), y(mi(i)+j))+bij((1-theta)*nu(i), ntest(mi(i)+j)-y(mi(i)+j)))*zij*zij+Ai(j, j)*(1-2*theta);
        Ai1(j)=Ai(j, j);
        
        xi=arma::zeros(r(i),p);
        for (int j=0;j<r(i);j++) {
          for (int q=0; q<p;q++) {
            xi(j,q)=x(mi(i)+j, q);
          }
        }
      }

      V=V+xi.t()*Ai1;
      W=W+xi.t()*Ci*xi;
    }
    betadraw.col(t)=betadraw.col(t-1)-(W-D).i()*(V-D*betadraw.col(t-1));
  }

  
  double logdetHess=logdethess(D-W);
  double loglhH0=0;
  
  for (int i=0;i<n;i++) {
    for (int j=0; j<r(i); j++) {
      eta=0;
      for (int q=0;q<p;q++) {
        eta=eta+x(mi(i)+j, q)*betadraw(q, nrun-1);
      }
      theta=1/(1+exp(-eta));
      loglhH0=loglhH0+lgamma(theta*nu(i)+y(mi(i)+j))+lgamma((1-theta)*nu(i)+ntest(mi(i)+j)-y(mi(i)+j))
        -lgamma(nu(i)+ntest(mi(i)+j))-lgamma(theta*nu(i))-lgamma((1-theta)*nu(i))+lgamma(nu(i));
    }
  }
  double sumdnormbeta=-0.5*p*log(sigma2beta);
  for (int u=0;u<p;u++) {
    sumdnormbeta=sumdnormbeta-0.5*pow(betadraw(u, nrun-1), 2)/sigma2beta;
  }

  double loglaplaceH0=loglhH0+sumdnormbeta-0.5*logdetHess;
  
  return Rcpp::List::create(
    Rcpp::Named("betadraw") = betadraw,
    Rcpp::Named("logdetHess") = logdetHess,    
    Rcpp::Named("loglaplaceH0") = loglaplaceH0
 );
}