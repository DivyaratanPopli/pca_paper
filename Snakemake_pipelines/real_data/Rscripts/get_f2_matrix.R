library(admixtools)
library(tidyverse)

get_fstats <- function(prefix, indf, outf2, outstd, f2_dir){
  
  prefix = strsplit(prefix, split=".geno")[[1]]
  
  if ( !dir.exists(f2_dir) ) {
    admixtools::extract_f2(pref=prefix, outdir=f2_dir)
  }
  
  f2s = admixtools::read_f2(f2_dir)
  
  allpop <- read.delim(indf, sep="\t", header = FALSE)
  pop=allpop$V3
  
  n <- length(pop)
  f2mat_mean <- matrix(NA, nrow = n, ncol = n)
  f2mat_std <- matrix(NA, nrow = n, ncol = n)
  
  # Step 4: Fill the matrix using a loop
  for (i in 1:n) {
    for (j in 1:n) {
      if (i==j){
        f2mat_mean[i, j] <- 0
        f2mat_std[i, j] <- 0
      }
      else{
        f2mat_mean[i, j] <- mean(f2s[pop[i], pop[j],])
        f2mat_std[i, j] <- sd(f2s[pop[i], pop[j],])
      }
      
    }
  }
  
  rownames(f2mat_mean) <- pop
  rownames(f2mat_std) <- pop
  
  write.table(f2mat_mean, file = outf2, row.names=TRUE, col.names=TRUE, sep=',')
  write.table(f2mat_std, file = outstd, row.names=TRUE, col.names=TRUE, sep=',')
  
}

#prefix = '/mnt/diversity/divyaratan_popli/pca_paper/Snakemake_pipelines/real_data/toy_example/all_ind.geno'
#f2_dir = '/mnt/diversity/divyaratan_popli/pca_paper/Snakemake_pipelines/real_data/to_delete/admixtools2_fmat'
#indf= '/mnt/diversity/divyaratan_popli/pca_paper/Snakemake_pipelines/real_data/toy_example/all_ind.ind'
get_fstats(prefix=snakemake@input[["prefix"]], indf=snakemake@input[["indf"]], outstd=snakemake@output[["outstd"]],
           f2_dir=snakemake@params[["f2_dir"]],outf2=snakemake@output[["outf2"]], )
