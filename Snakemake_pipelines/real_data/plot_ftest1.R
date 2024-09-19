library(ggplot2)
library(dplyr)

#fname="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations/simfiles/ftest_npop4_nind80/ppca_direct_val/ll_20_p1pop1_p2pop2_p3pop3_p4pop4.csv"

plot_f <- function(flist,fscale2,slist,s_scale2){

  data=read.csv(file=fname, header=F, sep=',')
  mat=as.data.frame(data)
  colnames(mat) <- c("migration_rate", "run", "loglikelihood_ratio")
  mat$migration_rate <- as.factor(mat$migration_rate)
  
  mat.summary <- mat %>%
    group_by(migration_rate) %>%
    summarise(
      sd = sd(loglikelihood_ratio, na.rm = TRUE),
      ll = mean(loglikelihood_ratio)
    )
  mat.summary
 
  ggplot(data=mat, aes(x=migration_rate, y=loglikelihood_ratio)) +
    geom_jitter(position = position_jitter(0.2), color="darkgrey") + 
    geom_errorbar(data=mat.summary,aes(ymax = loglikelihood_ratio + sd, ymin=loglikelihood_ratio - sd))
  
 
}