library(MASS)

simple_sim <- function(N_SAMPLES,TRUE_RANK,N_SNPS,SIGMA){

  points = matrix(runif(N_SAMPLES * TRUE_RANK), nrow=N_SAMPLES)
  points <- rbind(points)

  C = points %*% t(points)
  P = (outer(diag(C), diag(C)) + C^2)
  Q = outer(diag(C), diag(C), `+`)  #useful matrix for later


  X1 = mvrnorm(N_SNPS, rep(0.5, N_SAMPLES), C)
  N1 = nrow(X1)

  #data with added noise
  X3 = X1 #+ rnorm(length(X1), 0.5, sqrt(SIGMA))
  X4 = pmax(pmin(X3, 1), 0)
  #d2=data + mvrnorm(N_SNPS, mu=rep(0, N_SAMPLES), Sigma=diag(NOISE, N_SAMPLES))
  X5 = t(apply(X4, 1, function(p)rbinom(ncol(X4), 2, p)))
  return(X5)
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

  e = eigen(cm)
  pcs = e$vectors %*% diag(sqrt(abs(e$values)) * sign(e$values)  )


  return(list(pcs=pcs, eval=e$values))
}


ppca_cor <- function(data, n_pcs=4, normalized=T){
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




fstats_pca <- function(X3,method,npcs){

  data11=X3
  data12=data11/2
  mu = rowMeans(data11,na.rm=TRUE)
  data = data11 - mu
  data = data/2

  if (method=='pca'){
    xxx=run_pca(data)
  }else if (method=='ppca_direct'){
    xxx = ppca_direct(as.matrix(data), n_pcs=npcs)
  }else if (method=='PCA1'){
    xxx = PCA1(as.matrix(data12), n_pcs=npcs, normalized=TRUE)
  }else if (method=='ppca_cor'){
    xxx = ppca_cor(as.matrix(data12), n_pcs=npcs, normalized=TRUE)
  }

  return(xxx)

}

#genf = "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/test_samplesize100/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.geno_pc"


plotf <- function(genf,outf){

  X=read.csv(file=genf, header=F, sep=',')

  N_PCS=8
  xxx= fstats_pca(X3=X, method='pca', npcs=N_PCS)
  xxx1= fstats_pca(X3=X, method='ppca_direct', npcs=N_PCS)
  xxx2= fstats_pca(X3=X, method='PCA1', npcs=N_PCS)
  xxx3= fstats_pca(X3=X, method='ppca_cor', npcs=N_PCS)

  eigenvalues=xxx$eval[1:19]
  PC=1:19
  #plot(PC, eigenvalues, main="PCA", ylim=c(0,1.7))
  eigenvalues1=c(sqrt(xxx1$eval),rep(0,11))
  #plot(PC, eigenvalues1, main="Probabilistic PCA", ylim=c(0,1.7))
  eigenvalues2=xxx2$eval[1:19]
  eigenvalues3=c(sqrt(xxx3$eval),rep(0,11))

  png(file=outf, width=1024, height=600, res=200)
  par(mfrow=c(1,3))
  plot(PC, eigenvalues, main="Classical PCA \n(No noise)",ylim=c(0,1.3), xlab="", ylab="Eigenvalues", pch=19, cex.main=0.9)
  abline(h=0, col="red", lwd=2, lty=2)
  plot(PC, eigenvalues1, main="Probabilistic PCA \n(Normal homoskedastic noise)",ylim=c(0,1.3), ylab="", pch=19, cex.main=0.9)
  abline(h=0, col="red", lwd=2, lty=2)
  plot(PC, eigenvalues2, main="LSE \n(Binomial noise)",ylim=c(0,1.3), ylab="", xlab="", pch=19, cex.main=0.9)
  abline(h=0, col="red", lwd=2, lty=2)
  #plot(PC, eigenvalues3, main="PPCA with binomial error")
  dev.off()

}

plotf(genf=snakemake@input[["genf"]],outf=snakemake@output[["outf"]])
