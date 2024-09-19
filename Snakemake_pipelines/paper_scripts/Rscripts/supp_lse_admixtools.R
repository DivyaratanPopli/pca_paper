library(ggplot2)


#admixf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/mu0.05/run1/npop10_nind100/admixtools2Norm/f2mat102"
#lsef=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/mu0.05/run1/npop10_nind100/PCA1_val/f2mat104","/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/mu0.05/run1/npop10_nind100/PCA1_val/f2mat103","/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/mu0.05/run1/npop10_nind100/PCA1_val/f2mat102")

plotf <- function(admixf, lsef, truef, outplotf){

  admix = read.csv(file=admixf, header=F, sep=',')[1,4]
  tr = read.csv(file=truef, header=F, sep=',')[1,4]

  xx1=c()
  for(fi in seq(1,length(lsef))){

    xx = read.csv(file=lsef[fi], header=F, sep=',')[1,4]
    xx1=append(xx1,xx)

  }
  xx1=xx1/116762
  admix=admix/116762
  tr=tr/116762
  png(file=outplotf, res=200, height=5, width=8, units = "in")
  plot(xx1, ylab=expression('f2(X'[1]*',X'[4]*')'), xlab="Number of PCs")
  abline(h=admix, col="red")
  abline(h=tr, col="black", lty=2)
  legend("bottomright", legend=c("LSE", "ADMIXTOOLS 2", "True value"),
         col=c("black", "red", "black"), lty=c(NA,1,2), pch=c(1,NA, NA))
  dev.off()

}

plotf(admixf=snakemake@input[["admixf"]],
        truef=snakemake@input[["truef"]],
       lsef=snakemake@input[["lsef"]],
       outplotf=snakemake@output[["outplotf"]])
