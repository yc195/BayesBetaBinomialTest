##This function implements a Bayesian Beta-binomial regression test for multivariate count

##y is the matrix representing the number of successes
##ntest is the matrix representing the number of trials
##nrun is the number of iterations run for MLE
##nq is the number of nodes for Gaussian-Laguerre quadrature

require(statmod)

##calculates the leave one out Bayes factor in favor of the full model for the variable at position pos
bbtestLOO = function(y,ntest,x,pos,sigma=5,a=10,b=0.1,nq=15)
{
  points=gauss.quad(nq,"laguerre",alpha=a-1)$nodes
  weights=gauss.quad(nq,"laguerre",alpha=a-1)$weights
  priors=exp((1-b)*points)
  H0=rep(0,nq)
  H1=rep(0,nq)
  s=rep(1,p)
  s[pos]=0
  
  for(t in 1:nq)
  {
    H0[t]=calcll(y, x, rep(1,p), nrun, ntest, points[t], sigma2)
    H1[t]=calcll(y, x, s, nrun, ntest, points[t], sigma2)
  }
  
  logBF = log(sum(weights*exp(H1-max(c(H0,H1)))*priors))-log(sum(weights*exp(H0-max(c(H0,H1)))*priors))
  
  return(exp(logBF))
}