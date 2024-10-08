init_sigma <- function() {1}
init_cov = function(n_dims,n_pcs)diag(1, nrow=n_dims, ncol=n_pcs)

run_pca <- function(data){
  z = prcomp(data)
  return(unname(t(t(z$rotation) * z$sdev)))
}

ppca_direct <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
  n_dims = ncol(data)
  #eq 31
  mu = rowMeans(data,na.rm=TRUE)
  centered_data = data - mu

  S = t(centered_data) %*% centered_data / nrow(data)
  SVD = svd(S)

  #eq 8
  sigma = sum(SVD$d[(n_pcs+1):n_dims]) / (n_dims - n_pcs)

  #eq 7
  W = SVD$u[,1:n_pcs] %*% sqrt(diag(SVD$d[1:n_pcs] - sigma))
  e = svd(W)
  pcs = e$u %*% diag(e$d)
  pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))
  W = cbind(W, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))


  return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u))
}

ppca_EM <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
  n_dims = ncol(data)
  #eq 31

  S = t(data) %*% data / nrow(data)

  sigma = init_sigma()
  W = init_cov(n_dims, n_pcs)

  score = Inf

  for(i in 1:n_iter){
    #M is defined after eq 6
    M = t(W) %*% W + diag(sigma, n_pcs)
    M_inv = solve(M)

    #eq 29
    aux = diag(sigma, n_pcs) + M_inv %*% t(W) %*% S %*% W
    W_new = S %*% W %*% solve(aux)

    #eq 30
    sigma_aux = S - S %*% W %*% M_inv %*% t(W_new)
    sigma = 1 / n_dims * sum(diag(sigma_aux))

    W = W_new


    #simple convergence test
    new_score = sum(abs(W)) + sigma
    if(verbose)print(sprintf("iter: %d, score : |%.5f|",i, score-new_score))
    if(abs(score - new_score) < tol){
      break
    } else {
      score = new_score
    }


  }

  e = eigen(W %*% t(W))
  pcs = t(sqrt(abs(e$values)) * t(e$vectors))

  #posterior sanity checks, eq 9
  M = t(W) %*% W + diag(sigma, n_pcs)
  M_inv = solve(M)

  posterior_mean = M_inv %*% t(W) %*% t(data )
  approx_data = W %*% solve(t(W) %*% W) %*% M %*% posterior_mean + mu

  return(list(W=W, sigma=sigma, pcs=pcs, eval=e$values, evec=e$vectors,
              data=approx_data, est_data=posterior_mean + mu))
}


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

  e = svd(W)
  pcs = e$u %*% diag(e$d)
  pcs = cbind(pcs, matrix(0, ncol=D-nComp, nrow=D))
  #W = cbind(W, matrix(0, ncol=n_dims-nComp, nrow=D))

  if (showits)                                                     # Show last iteration
    cat(paste0(format(it), "...", "\n"))


  return(list(loadings=pcs, W=W, ll=ll, sigma2=sigma2))
}

ppca_r <- function(data, npc){
  library(pcaMethods)
  D=ncol(data)
  result <- pca(data, method="ppca", nPcs=20)
  res1=loadings(result)
  res2=cbind(res1, matrix(0, ncol=D-npc, nrow=D))
  return(res2)
}

double_center <- function(dist)
  t(dist - rowMeans(dist)) - rowMeans(dist) + mean(dist)

PCA1 <- function(data, n_pcs=4, tol=1e-4, verbose=T, normalized=T){
    n_dims = ncol(data)
    n_snps = nrow(data)

    if(normalized){ #such that data is between 0 and 1
        dhat = colMeans( data * (1-data))
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat)
    } else{ #unnormalized version
        dhat = colMeans( data * (2-data))
        gram_mat <- t(data) %*% data / n_snps - diag(dhat)
        cm = double_center(gram_mat)
    }

    e = svd(cm)
    pcs = e$u %*% sqrt(diag(e$d)  )


    return(list(pcs=pcs, eval=e$d))
}

ppca_cor <- function(data, n_pcs=4, normalized=F){
  n_dims = ncol(data)
  n_snps = nrow(data)

  if(normalized){ #such that data is between 0 and 1
    dhat = colMeans( data * (1-data))
    gram_mat <- t(data) %*% data / n_snps - diag(dhat)
    cm = double_center(gram_mat)
  } else{ #unnormalized version
    dhat = colMeans( data * (2-data))
    gram_mat <- t(data) %*% data / n_snps - diag(dhat)
    cm = double_center(gram_mat)
  }

  SVD = svd(cm)

  #eq 8
  sigma = sum(SVD$d[(n_pcs+1):n_dims]) / (n_dims - n_pcs)

  #eq 7
  W = SVD$u[,1:n_pcs] %*% sqrt(diag(SVD$d[1:n_pcs] - sigma))
  e = svd(W)
  pcs = e$u %*% diag(e$d)
  pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))


  return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u))
}



fstats_pca <- function(dataf, outfile, method, npcs){

  data=read.csv(file=dataf, header=F, sep=',')
  is.na(data) <- data == "9"

  data=read.csv(file=dataf, header=F, sep=',')
  data=data/2
  mu = rowMeans(data,na.rm=TRUE)

  if (method=='pca'){
    data = data - mu
    pc=t(run_pca(data))
  }else if (method=='ppca_direct'){
    data = data - mu
    pc=t(ppca_direct(as.matrix(data), n_pcs=npcs)$pcs)
  }else if (method=='ppca_EM'){
    data = data - mu
    pc=t(ppca_EM(as.matrix(data),n_pcs=npcs)$pcs)
  }else if (method=='ppca_miss'){
    data = data - mu
    pc=t(PPCAEM(X=as.matrix(data),nComp=npcs)$loadings)
  }else if (method=='ppca_r'){
    data = data - mu
    pc=t(ppca_r(data=as.matrix(data), npc=npcs))
  }else if (method=='PCA1'){
    pc=t(PCA1(data=as.matrix(data, n_pcs=npcs, normalized=T))$pcs)
  }else if (method=='ppca_cor'){
    pc=t(ppca_cor(data=as.matrix(data), n_pcs=npcs, normalized=T)$pcs)
  }

    write.table(pc, file = outfile, row.names=FALSE, col.names=FALSE, sep=',')

  }

  fstats_pca(dataf=snakemake@input[["dataf"]], outfile=snakemake@output[["outfile"]],
             method=snakemake@wildcards[["method"]], npcs=as.integer(snakemake@wildcards[["npcs"]]))
