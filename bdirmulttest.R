##This function implements a Bayesian Dirichlet-multinomial test for two groups of multivariate count data

##n1 and n2 are two groups of samples
##a and b are hyperparameters for the within-group variations eta1, eta2 \sim Ga(a,b)
bdirmulttest = function(n1,n2,a,b)
{
  K=dim(n1)[2]
  flambda_min1=optim(par=rep(1, (2*K)), fn=funcH1,  n1=n1,n2=n2, a1=a, b1=b, method="BFGS", hessian=TRUE)  ## under H1
  flambda1=flambda_min1$value
  H1=flambda_min1$hessian
  
  flambda_min0=optim(par=rep(1, (K+1)), fn=funcH0,  n1=n1,n2=n2, a0=a, b0=b, method="BFGS", hessian=TRUE)  ## under H0
  flambda0=flambda_min0$value
  H0=flambda_min0$hessian
  
  logBF=(-flambda1+flambda0+(K-1)*log(2*pi)/2+(unlist(determinant(H0,logarithm=TRUE))[1]-unlist(determinant(H1,logarithm=TRUE))[1])/2-mlbeta(rep(1/K, K)))[[1]]
  print(exp(logBF))
}