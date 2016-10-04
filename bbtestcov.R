##This function implements a Bayesian Beta-binomial test with covariates and overdispersion

##y is the matrix representing the number of successes
##ntest is the matrix representing the number of trials
##r is a vector that gives the total number of samples in each group
##x is the matrix representing the covariates
##nrun is the number of iterations run for MLE
##nq is the number of nodes for Gaussian-Laguerre quadrature

##Example:
##y = c(1,4,2,4,7)
##ntest = c(3,6,3,5,10)
##x = matrix(c(1,1,1,1,1,7,5,3,2,4),nrow=5,ncol=2)
##r = c(2,3)
##BF = bbtestL(y,ntest,x,r)

require(statmod)

bbtestL = function(y,ntest,x,r,sigma2beta=5,sigma2gamma=5,a=10,b=0.1,nq=10)
{
  p=dim(x)[[2]]-1
  points=gauss.quad(nq,"laguerre",alpha=a-1)$nodes
  weights=gauss.quad(nq,"laguerre",alpha=a-1)$weights
  priors=exp((1-b)*points)
  H0=rep(0,nq)
  H1=rep(0,nq)

  for(t in 1:nq)
  {
    H0[t]=calcloglhH0(y, x, rep(1,p), r, nrun, ntest, rep(points[t],length(r)), sigma2beta)$loglaplaceH0
    H1[t]=calcloglhH1(y, x, rep(1,p), r, nrun, ntest, rep(points[t],length(r)), sigma2beta, sigma2gamma)$loglaplaceH1
  }
  
  logBF = log(sum(weights*exp(H1-max(c(H0,H1)))*priors))-log(sum(weights*exp(H0-max(c(H0,H1)))*priors))
  
  return(exp(logBF))
}