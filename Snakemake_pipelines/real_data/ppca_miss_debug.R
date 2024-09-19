
X=data
nComp=20
PPCAEM = function(X, nComp=3, tol=.00001, maxits=100, showits=T){
  # Arguments X: numeric data, nComp: number of components
  # tol = tolerance level, maxits: maximum iterations, showits: show iterations
  require(pracma) # for orthonormal basis of W; pcaMethods package has also
  require(psych)  # for tr
  
  # starting points and other initializations
  Xorig = X
  X = X
  N = nrow(Xorig)
  D = ncol(Xorig)
  L = nComp
  NAs = is.na(Xorig)
  
  X[NAs] = 0
  
  mu = rowMeans(X)
  X = X - mu
  
  S = (1/N) * t(X)%*%X
  evals = eigen(S)$values
  evecs = eigen(S)$vectors
  
  V = evecs[,1:L]
  Lambda = diag(evals[1:L])
  
  Z = t(replicate(L, rnorm(N)))                                    # latent variables
  sigma2 = 1/(D-L) * sum(evals[(L+1):D])                           # variance; average variance associated with discarded dimensions
  W = V %*% chol(Lambda-sigma2*diag(L)) %*% diag(L)                # loadings
  
  it = 0
  converged = FALSE
  ll = 0
  
  if (showits)                                                     # Show iterations
    cat(paste("Iterations of EM:", "\n"))
  while ((!converged) & (it < maxits)) {
    if(exists('W.new')){
      W.old = W.new
    }
    else {
      W.old = W
    }
    
    ll.old = ll
    
    proj = t(W.old%*%Z)
    Xnew = Xorig
    Xnew[NAs] = proj[NAs]
    X = Xnew
    
    Psi = sigma2*diag(L)
    M = t(W.old) %*% W.old + Psi
    
    W.new = S%*%W.old%*%solve(Psi + solve(M)%*%t(W.old)%*%S%*%W.old)   # E and M
    sigma2 = 1/D * tr(S - S%*%W.old%*%solve(M)%*%t(W.new))
    
    Z = solve(M)%*%t(W.new)%*%t(X)
    
    
    # log likelihood as in paper
    #     ZZ = sigma2*solve(M) + Z%*%t(Z)
    #     ll = .5*sigma2*D + .5*tr(ZZ) + .5*sigma2 * X%*%t(X) -
    #          1/sigma2 * t(Z)%*%t(W.new)%*%t(X) + .5*sigma2 * tr(t(W.new)%*%W.new%*%ZZ)
    #     ll = -sum(ll)
    
    # more straightforward
    ll = dnorm(X, mean=t(W.new%*%Z), sd=sqrt(sigma2), log=T)
    ll = -sum(ll)
    
    it = it + 1
    if (showits & (it == 1 | it%%5 == 0))                          # if showits, show first and every 5th iteration
      cat(paste(format(it), "...", "\n", sep = ""))
    converged = max(abs(ll.old-ll)) <= tol
  }
  
  W = orth(W.new)
  evs = eigen(cov(X %*% W))
  evecs = evs$vectors
  
  W = W %*% evecs
  Z = X %*% W
  Xrecon = Z %*% t(W)
  reconerr = sum((Xrecon-X)^2)
  
  if (showits)                                                     # Show last iteration
    cat(paste0(format(it), "...", "\n"))
  
  W1 = cbind(W, matrix(0, ncol=D-nComp, nrow=D))
  return(list(scores=Z, loadings=W1, Xrecon=Xrecon, reconerr=reconerr, ll=ll, sigma2=sigma2))
}



library(pcaMethods)

## Load a sample metabolite dataset with 5\% missing values (metaboliteData)
data(metaboliteData)
## Perform probabilistic PCA using the 3 largest components
result <- pca(t(metaboliteData), method="ppca", nPcs=3, seed=123)
## Get the estimated complete observations
cObs <- completeObs(result)
## Plot the scores
plotPcs(result, type = "scores")



