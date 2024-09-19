library(ggplot2)
library(tidyr)

#fname="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/migrations_comparison/simfiles/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/ppca_direct_val_scale20/ll.csv"

plot_comp <- function(fname,fout){

  data=read.csv(file=fname, header=F, sep=',')

  mat=as.data.frame(data)
  colnames(mat) <- c("migration_rate", "run", "ppca", "admixtools")
  mat$migration_rate <- as.factor(mat$migration_rate)


  gg_mat <- mat %>% pivot_longer(
    cols = c("ppca","admixtools")
  )


  ggplot(data=gg_mat, aes(x=migration_rate, y=value, col=name)) +
    geom_jitter(position = position_jitter(0.2)) +
    #ylim(c(-10,50)) +
    geom_hline(yintercept=c(3,-3), col="black", linetype="dashed") +
    geom_hline(yintercept=c(2,-2), col="black", linetype="dashed") +

    xlab("Migration rate") + ylab("Standard deviations away from 0") +
    scale_y_continuous(breaks=c(-3,-2,0,2,3))+
    scale_colour_discrete("Method",labels=c('admixtools2', 'ppca_scale8'))  +
    theme_classic()

  ggsave(fout,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
}

plot_comp(fname=snakemake@input[["fname"]],
         fout=snakemake@output[["fout"]])
