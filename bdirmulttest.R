##This function implements a Bayesian Dirichlet-multinomial test for two groups of multivariate count data

##n1 and n2 are two groups of samples
##a and b are hyperparameters for the within-group variations eta1, eta2 \sim Ga(a,b)
##nq is the number of nodes for Gaussian-Laguerre quadrature
bdirmulttest = function(n1,n2,a,b,nq)
{
  K=dim(n1)[2]
  points=gauss.quad(nq,"laguerre",alpha=a-1)$nodes
  weights=gauss.quad(nq,"laguerre",alpha=a-1)$weights
  priors=as.matrix(exp((1-b)*points))%*%exp((1-b)*points)
  H0st=matrix(NA,nrow=nq,ncol=nq)
  H1st=matrix(NA,nrow=nq,ncol=nq)

  for(t in 1:nq)
  {
    for(s in 1:nq)
    {
      flambda_min1=optim(par=rep(1, (2*K-2)), fn=logitnoeta1,  n1=n1,n2=n2, eta1=points[s], eta2=points[t], method="BFGS", hessian=TRUE)  ## under H1
      flambda_min0=optim(par=rep(1, (K-1)), fn=logitnoeta0,  n1=n1,n2=n2, eta1=points[s], eta2=points[t], method="BFGS", hessian=TRUE)  ## under H0
      
      flambda1=flambda_min1$value
      H1=flambda_min1$hessian
      flambda0=flambda_min0$value
      H0=flambda_min0$hessian
      
      H1st[s,t]=-flambda1-unlist(determinant(H1,logarithm=TRUE))[1]/2
      H0st[s,t]=-flambda0-unlist(determinant(H0,logarithm=TRUE))[1]/2
    }
  }
  logBF = (K-1)*log(2*pi)/2-mlbeta(rep(1/K, K))+log(sum(exp(H1st-max(c(H0st,H1st)))*priors))-log(sum(exp(H0st-max(c(H0st,H1st)))*priors))
  return(exp(logBF))
}
