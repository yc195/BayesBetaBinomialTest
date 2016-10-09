##This function implements a Bayesian Beta-binomial test with covariates and overdispersion

##y is the matrix representing the number of successes
##ntest is the matrix representing the number of trials
##r is a vector that gives the total number of samples in each group
##x is the matrix representing the covariates
##nrun is the number of iterations run for MLE
##nq is the number of nodes for Riemannian sum

##Example:
##y = c(1,4,2,4,7)
##ntest = c(3,6,3,5,10)
##x = matrix(c(1,1,1,1,1,7,5,3,2,4),nrow=5,ncol=2)
##r = c(2,3)
##BF = bbtestL(y,ntest,x,r)

bbtestL = function(y,ntest,x,r,sigma2beta=5,nq=10)
{
  p=dim(x)[[2]]-1
  sq = seq(1,nq,1)
  points=exp(-1+5*sq/nq)
  
  H0st=matrix(NA,nrow=nq,ncol=nq)
  H1st=matrix(NA,nrow=nq,ncol=nq)
  
  for(t in 1:nq)
  {
    for(s in 1:nq)
    {
      H0st[s,t]=calcloglhH0(y, x, rep(1,p), r, nrun, ntest, c(points[s],points[s]), sigma2beta)$loglaplaceH0
      H1st[s,t]=calcloglhH1(y, x, rep(1,p), r, nrun, ntest, c(points[s],points[s]), sigma2beta, points[t])$loglaplaceH1
    }
  }
  logBF = log(mean(exp((H1st-max(c(H0st,H1st))))))-log(mean(exp(H0st-max(c(H0st,H1st)))))
  return(exp(logBF))
}