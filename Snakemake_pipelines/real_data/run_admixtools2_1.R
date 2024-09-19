library(admixtools)
library(tidyverse)

get_fstats <- function(prefix, outfile, f2_dir){

  prefix = strsplit(prefix, split=".geno")[[1]]

  if ( !dir.exists(f2_dir) ) {
    admixtools::extract_f2(prefix, f2_dir)
  }

  f2s = admixtools::read_f2(f2_dir)

  pop1="Altai"
  pop2="Vindija33.19"
  pop3=c("Les_Cottes_L35MQ25","Goyet_L35MQ25","Mezmaiskaya1_L35MQ25","Mezmaiskaya2_L35MQ25","VindijaG1_L35MQ25","Spy_L35MQ25")

  f3stats = qp3pop(f2s, pop1, pop2, pop3)
  print(f3stats[,3])
  f3stats[,3] = c("Goyet Q56-1","Les Cottes Z4-1514","Mezmaiskaya1","Mezmaiskaya2","Spy 94a","Vindija 87")
  write.table(f3stats, file = outfile, row.names=FALSE, col.names=FALSE, sep=',')

  print(admixtools::f3(f2s, pop1, pop2, "Goyet_L35MQ25")$est)

}

#prefix = '/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations_comparison/simfiles/mu0/run1/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/eigen_sub.geno'
#f2_dir = '/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations_comparison/simfiles/mu0/run1/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/admixtools2_fmat'

get_fstats(prefix=snakemake@input[["prefix"]], outfile=snakemake@output[["outfile"]],
           f2_dir=snakemake@params[["f2_dir"]])
