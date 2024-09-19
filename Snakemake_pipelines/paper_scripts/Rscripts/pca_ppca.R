library(MASS)

simple_sim <- function(N_SAMPLES,TRUE_RANK,N_SNPS,SIGMA){

  points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
  points <- rbind(points)

  C = points %*% t(points)
  P = (outer(diag(C), diag(C)) + C^2)
  Q = outer(diag(C), diag(C), `+`)  #useful matrix for later


  X1 = mvrnorm(N_SNPS, rep(0, N_SAMPLES), C)
  N1 = nrow(X1)

  #data with added noise
  X3 = X1 + rnorm(length(X1), 0, sqrt(SIGMA))
  return(X3)
}

run_pca <- function(data){
  z = prcomp(data, center=F, scale=F)
  pcs = unname(t(z$rotation) * z$sdev)
  scores = z$x
  data_rec = z$x %*% (pcs/z$sdev)

  return(list(pcs = pcs, scores = scores, data_rec = data_rec, eval=z$sdev))
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

  #reconstructing data and score mat

  T1= W %*% solve(t(W) %*% W)
  M= t(W) %*% W + (sigma * diag(n_pcs))
  xntn=solve(M) %*% t(W) %*% t(data)

  T_rec= t(T1 %*% M %*% xntn)

  scor=t(xntn)

  e = svd(W)
  pcs = e$u %*% diag(e$d)
  pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))
  W = cbind(W, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))


  return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u,data_rec=T_rec, scores=scor))
}


fstats_pca <- function(X3,method,npcs){

  data=X3
  mu = rowMeans(data,na.rm=TRUE)
  data = data - mu

  if (method=='pca'){
    xxx=run_pca(data)
  }else if (method=='ppca_direct'){
    xxx = ppca_direct(as.matrix(data), n_pcs=npcs)
  }

  return(xxx)

}

plotf <- function(outf){

  N_SNPS = 100000
  N_SAMPLES = 20
  TRUE_RANK = 4
  N_PCS=4
  SIGMA=0.5

  X3= simple_sim(N_SAMPLES=N_SAMPLES , TRUE_RANK=TRUE_RANK , N_SNPS=N_SNPS, SIGMA=SIGMA)
  xxx= fstats_pca(X3=X3, method='pca', npcs=N_PCS)
  xxx1= fstats_pca(X3=X3, method='ppca_direct', npcs=N_PCS)


  eigenvalues=xxx$eval[1:19]
  PC=1:19
  #plot(PC, eigenvalues, main="PCA", ylim=c(0,1.7))
  eigenvalues1=c(sqrt(xxx1$eval),rep(0,15))
  #plot(PC, eigenvalues1, main="Probabilistic PCA", ylim=c(0,1.7))

  png(file=outf, width=600, height=350)
  par(mfrow=c(1,2))
  plot(PC, eigenvalues, main="PCA", ylim=c(0,1.7))
  plot(PC, eigenvalues1, main="Probabilistic PCA", ylim=c(0,1.7))

  dev.off()

}

plotf(outf=snakemake@output[["outf"]])
