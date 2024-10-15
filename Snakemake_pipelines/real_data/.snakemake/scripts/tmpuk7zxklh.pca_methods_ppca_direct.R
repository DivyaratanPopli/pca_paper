
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('toy_example/all_ind.geno_pc', 'to_delete/n_blocks.csv', "dataf" = 'toy_example/all_ind.geno_pc', "filter" = 'to_delete/n_blocks.csv'),
    output = list('to_delete/filter1276/pcs_ppca_miss/pcs_npcs2.csv', 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2.csv_scores', 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_data_rec.csv', 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_sigma.txt', 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_eval.csv', "outfile" = 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2.csv', "outfile_s" = 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2.csv_scores', "outfile_d" = 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_data_rec.csv', "sigmaf" = 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_sigma.txt', "evalf" = 'to_delete/filter1276/pcs_ppca_miss/pcs_npcs2_eval.csv'),
    params = list(),
    wildcards = list('1276', 'ppca_miss', '2', "f" = '1276', "method" = 'ppca_miss', "npcs" = '2'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list(),
    rule = 'ppcaMethods_sub',
    bench_iteration = as.numeric(NA),
    scriptdir = '/mnt/diversity/divyaratan_popli/pca_paper/Snakemake_pipelines/real_data/Rscripts',
    source = function(...){
        wd <- getwd()
        setwd(snakemake@scriptdir)
        source(...)
        setwd(wd)
    }
)


######## snakemake preamble end #########
init_sigma <- function() {1}
init_cov = function(n_dims,n_pcs)diag(1, nrow=n_dims, ncol=n_pcs)

run_pca <- function(data){
  z = prcomp(data)
  pcs = unname(t(z$rotation) * z$sdev)
  scores = z$x
  data_rec = z$x %*% (pcs/z$sdev)

  return(list(pcs = pcs, scores = scores, data_rec = data_rec))
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


ppca_miss <- function(data, n_pcs=4, n_iter=100, tol=1e-5, verbose=T){

  data_orig = data
  NAs = is.na(data_orig)
  data[NAs] = 0
  W=0
  i=0
  while(i<n_iter){

    res = ppca_direct(data=as.matrix(data), n_pcs=n_pcs)
    data[NAs] = res$data_rec[NAs]

    diff=sum((res$W - W)^2)
    print(diff)
    W=res$W
    i=i+1

    if(diff < tol){
      return(list(W=res$W, sigma=res$sigma, pcs=res$pcs, eval=res$eval,
                  evec=res$evec,data_rec=res$data_rec, scores=res$scores))
    }
  }
  return(list(W=res$W, sigma=res$sigma, pcs=res$pcs, eval=res$eval,
              evec=res$evec,data_rec=res$data_rec, scores=res$scores))

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

fstats_pca <- function(dataf, filter, f, outfile, outfile_s, outfile_d, method, npcs, sigmaf, evalf){

  data_un=read.csv(file=dataf, header=F, sep=',')
  fi=read.csv(file=filter, header=F, sep=',')

  rm_start = sum(fi[1:f-1,1])+1
  rm_end = sum(fi[1:f,1])
  data = data_un[-seq(rm_start, rm_end),]

  is.na(data) <- data == "9"
  data=data/2
  mu = rowMeans(data,na.rm=TRUE)
  data = data - mu

  if (method=='pca'){
    pc=run_pca(data)$pcs
    scores=run_pca(data)$scores
    data_rec=run_pca(data)$data_rec

  }else if (method=='ppca_direct'){
    xxxx=ppca_direct(as.matrix(data), n_pcs=npcs)
    pc=t(xxxx$pcs)
    scores=xxxx$scores
    data_rec=xxxx$data_rec
    sigma=xxxx$sigma
  }else if (method=='ppca_EM'){
    pc=t(ppca_EM(as.matrix(data),n_pcs=npcs)$pcs)
  }else if (method=='ppca_miss_old'){
    pc=t(PPCAEM(X=as.matrix(data),nComp=npcs)$loadings)
  }else if (method=='ppca_r'){
    pc=t(ppca_r(data=as.matrix(data), npc=npcs))
  }else if (method=='ppca_miss'){
    xxxx=ppca_miss(as.matrix(data), n_pcs=npcs)
    pc=t(xxxx$pcs)
    scores=xxxx$scores
    data_rec=xxxx$data_rec
    sigma=xxxx$sigma
    eval=xxxx$eval
  }



  write.table(pc, file = outfile, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(scores, file = outfile_s, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(data_rec, file = outfile_d, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(eval, file = evalf, row.names=FALSE, col.names=FALSE, sep=',')
  file <- file(sigmaf, "w")
  write(sigma,file)
  close(file)

}

fstats_pca(dataf=snakemake@input[["dataf"]], filter=snakemake@input[["filter"]],outfile=snakemake@output[["outfile"]],
           method=snakemake@wildcards[["method"]], npcs=as.integer(snakemake@wildcards[["npcs"]]), f=as.integer(snakemake@wildcards[["f"]]),
           outfile_s=snakemake@output[["outfile_s"]],outfile_d=snakemake@output[["outfile_d"]], sigmaf=snakemake@output[["sigmaf"]], evalf=snakemake@output[["evalf"]])

#dataf = '/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/missing/simfiles/mu0/run1/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/eigen_sub.geno_pc'
