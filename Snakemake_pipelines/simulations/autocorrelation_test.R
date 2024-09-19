library(sns)
library(dplyr)

#fname="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/test_hypothesis_admixTools_x/simfiles/mu0/run1/npop100_nind100/eigen.geno_pc"
#outfile="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/admixpops_split_time/test"
ess_r <- function(fname, outfile){
  geno=read.csv(file=fname, header=F, sep=',')

  geno1=geno%>%mutate(freq_0=rowSums(.==0))%>%mutate(freq_1=rowSums(.==1))%>%mutate(freq_2=rowSums(.==2))
  geno2=geno1[!((geno1["freq_0"]==0 & geno1["freq_1"]==0) | (geno1["freq_2"]==0 & geno1["freq_1"]==0)),] #removing monomorphic sites

  geno_s=geno2[!((geno2["freq_2"]==103 & geno2["freq_0"]==0) | (geno2["freq_0"]==103 & geno2["freq_2"]==0)),]
  geno_d=geno_s[!((geno_s["freq_2"]==103 & geno_s["freq_0"]==1) | (geno_s["freq_0"]==103 & geno_s["freq_2"]==1) | (geno_s["freq_0"]==102 & geno_s["freq_1"]==2) | (geno_s["freq_2"]==102 & geno_s["freq_1"]==2) ),]



  ncol=mean(ess(geno_d))
  #ncol=mean(ess(data))
  write(ncol,file=outfile)
}

ess_r(fname=snakemake@input[["fname"]], outfile=snakemake@output[["outfile"]])
