library(admixtools)
library(tidyverse)

get_fstats <- function(prefix, outfile, outf2, f2_dir, pop1, pop2, pop3, pop4){

  prefix = strsplit(prefix, split=".geno")[[1]]

  if ( !dir.exists(f2_dir) ) {
    admixtools::extract_f2(pref=prefix, outdir=f2_dir, maxmiss=1)
  }

  f2s = admixtools::read_f2(f2_dir)


  f2_mean=mean(f2s[pop1,pop2,])
  f3_mean=admixtools::f3(f2s, pop2, pop1, pop3)$est
  f4_mean=admixtools::f4(f2s, pop1, pop2, pop3, pop4)$est
  f4_std=admixtools::f4(f2s, pop1, pop2, pop3, pop4)$se
  f2mat1=apply(f2s, 1:2, mean)

  f2mat=matrix(c(f2mat1[pop1,pop1],f2mat1[pop2,pop1],f2mat1[pop3,pop1],f2mat1[pop4,pop1],
  f2mat1[pop1,pop2],f2mat1[pop2,pop2],f2mat1[pop3,pop2],f2mat1[pop4,pop2],
  f2mat1[pop1,pop3],f2mat1[pop2,pop3],f2mat1[pop3,pop3],f2mat1[pop4,pop3],
  f2mat1[pop1,pop4],f2mat1[pop2,pop4],f2mat1[pop3,pop4],f2mat1[pop4,pop4]), nrow = 4, ncol = 4, byrow = FALSE)

  #f2_mean=f2mat[pop1,pop2]
  #f3_mean=(f2mat[pop2, pop3] + f2mat[pop1, pop3] - f2mat[pop1, pop2])/2
  #f4_mean=(f2mat[pop1, pop3] + f2mat[pop2, pop4] - f2mat[pop1, pop2] - f2mat[pop3, pop4])/2


  write.table(c(f2_mean, f3_mean, f4_mean, f4_std), file = outfile, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(round(f2mat,4), file = outf2, row.names=FALSE, col.names=FALSE, sep=',')

}

#prefix = '/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations_comparison/simfiles/mu0/run1/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/eigen_sub.geno'
#f2_dir = '/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations_comparison/simfiles/mu0/run1/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/admixtools2_fmat'

get_fstats(prefix=snakemake@input[["prefix"]], outfile=snakemake@output[["outfile"]],
           f2_dir=snakemake@params[["f2_dir"]],outf2=snakemake@output[["outf2"]],
           pop1=snakemake@params[["pop1"]], pop2=snakemake@params[["pop2"]],
           pop3=snakemake@params[["pop3"]], pop4=snakemake@params[["pop4"]])
