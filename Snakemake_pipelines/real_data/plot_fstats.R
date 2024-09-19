
#flist=c("/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/noisy_val/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv",
 #       "/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/true_val/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv",
  #      "/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/pca_val/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv",
   #     "/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/ppca_direct_val/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv",
    #    "/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/ppca_EM_val/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv",
     #   "/mnt/diversity/divyaratan_popli/fstats/preliminary/pipeline/f_statistics/nsnp100000_nind6_rank3_noise0.05_missing0/admixtools2/fstats_scale4_p1pop1_p2pop2_p3pop3_p4pop4.csv")
plot_f <- function(flist, n_pcs, outplot, S1,S2,S3,S4, n_pcs1, n_pcs2, fscale2){

  noisy=read.csv(file=flist[1], header=F, sep=',')
  true=read.csv(file=flist[2], header=F, sep=',')
  pca=read.csv(file=flist[3], header=F, sep=',')
  ppca_direct=read.csv(file=flist[4], header=F, sep=',')
  ppca_EM=read.csv(file=fscale2, header=F, sep=',')
  admix=read.csv(file=flist[5], header=F, sep=',')
  adjusted=read.csv(file=flist[6], header=F, sep=',')
  N_SAMPLES=dim(pca)[2]

  png(file=outplot,
      width=1000, height=500)

  trx=N_SAMPLES 
  nosx=N_SAMPLES -(N_SAMPLES/50)
  adjx=N_SAMPLES -(2*N_SAMPLES/50)
  admx=N_SAMPLES -(3*N_SAMPLES/50)

  par(mfrow=c(1,3), cex=1.2)
  require(glue)
  stat="f2"
  pos=1
  plot(NA,
       xlab="# PCS", ylab=stat,
       xlim=c(0, N_SAMPLES), ylim=c(min(noisy[pos,], true[pos,], min(pca[pos,]), min(ppca_direct[pos,]), min(ppca_EM[pos,]), admix[pos,], adjusted[pos,]), max(noisy[pos,], true[pos,], max(pca[pos,]), max(ppca_direct[pos,]), max(ppca_EM[pos,]), admix[pos,],adjusted[pos,])))
  title(sprintf("f2(%s, %s)", S1, S2))
  #abline(h=c(true[pos,], noisy[pos,], admix[pos,], adjusted[pos,]), col=1:4, lty=2)
  points(x=trx, y=true[pos,], col="black", cex=1, pch=10)
  points(x=nosx, y=noisy[pos,], col="red", cex=1, pch=3)
  points(x=adjx, y=adjusted[pos,], col="green", cex=1, pch=4)
  points(x=admx, y=admix[pos,], col="blue", cex=1, pch=7)
  lines(t(pca[pos,]),  col=5,lwd=2)
  lines(t(ppca_direct[pos,]),  col=6)
  lines(t(ppca_EM[pos,]),  col=7)

  stat="f3"
  pos=2
  plot(NA,
       xlab="# PCS", ylab=stat,
       xlim=c(0, N_SAMPLES), ylim=c(min(noisy[pos,], true[pos,], min(pca[pos,]), min(ppca_direct[pos,]), min(ppca_EM[pos,]), admix[pos,], adjusted[pos,]), max(noisy[pos,], true[pos,], max(pca[pos,]), max(ppca_direct[pos,]), max(ppca_EM[pos,]), admix[pos,],adjusted[pos,])))
  title(sprintf("f3(%s, %s, %s)", S1, S2, S3))
  #abline(h=c(true[pos,], noisy[pos,], admix[pos,], adjusted[pos,]), col=1:4, lty=2)
  points(x=trx, y=true[pos,], col="black", cex=1, pch=10)
  points(x=nosx, y=noisy[pos,], col="red", cex=1, pch=3)
  points(x=adjx, y=adjusted[pos,], col="green", cex=1, pch=4)
  points(x=admx, y=admix[pos,], col="blue", cex=1, pch=7)
  lines(t(pca[pos,]),  col=5,lwd=2)
  lines(t(ppca_direct[pos,]),  col=6)
  lines(t(ppca_EM[pos,]),  col=7)


  stat="f4"
  pos=3
  plot(NA,
       xlab="# PCS", ylab=stat,
       xlim=c(0, N_SAMPLES), ylim=c(min(noisy[pos,], true[pos,], min(pca[pos,]), min(ppca_direct[pos,]), min(ppca_EM[pos,]), admix[pos,], adjusted[pos,]), max(noisy[pos,], true[pos,], max(pca[pos,]), max(ppca_direct[pos,]), max(ppca_EM[pos,]), admix[pos,],adjusted[pos,])))
  title(sprintf("f4(%s, %s; %s, %s)", S1, S2, S3, S4))
  #abline(h=c(true[pos,], noisy[pos,], admix[pos,], adjusted[pos,]), col=1:4, lty=2)
  points(x=trx, y=true[pos,], col="black", cex=1, pch=10)
  points(x=nosx, y=noisy[pos,], col="red", cex=1, pch=3)
  points(x=adjx, y=adjusted[pos,], col="green", cex=1, pch=4)
  points(x=admx, y=admix[pos,], col="blue", cex=1, pch=7)
  lines(t(pca[pos,]),  col=5,lwd=2)
  lines(t(ppca_direct[pos,]),  col=6)
  lines(t(ppca_EM[pos,]),  col=7)

  legend("bottomleft", col=c(5,6,7,1,2,3,4),
         legend=c("pca", glue("pPCA rank {n_pcs1}"),
                  glue("pPCA rank {n_pcs2}"),
                  "truth", "full noisy data", "adjusted", "admixtools2"),
         lty=c(1,1,1,NA,NA,NA,NA), pch = c(NA, NA, NA, 10, 3, 4, 7))

  dev.off()
}

plot_f(flist=snakemake@input[["flist"]], n_pcs1=snakemake@wildcards[["npcs"]], outplot=snakemake@output[["outplot"]], S1=snakemake@wildcards[["pop1"]], S2=snakemake@wildcards[["pop2"]], S3=snakemake@wildcards[["pop3"]], S4=snakemake@wildcards[["pop4"]],
  n_pcs2 = snakemake@wildcards[["npcs2"]], fscale2=snakemake@input[["fscale2"]])
