library(dplyr)
library(ggplot2)

#admf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/admixtools2/f3stats.csv"
#ppca_muf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/mu.csv"
#ppca_stdf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/std.csv"

plotf3 <- function(admf, ppca_muf, ppca_stdf, plotf){

  adm=read.csv(file=admf, header=F, sep=',')
  ppca_mu=read.csv(file=ppca_muf, header=F, sep=',')

  ppca_std=read.csv(file=ppca_stdf, header=F, sep=',')

  ppca1 = as.data.frame(cbind(ppca_mu, ppca_std))
  colnames(ppca1)=c('mu','std')
  ppca1$pop1=rep('Altai',6)
  ppca1$pop2=rep('Vindija33.19',6)
  ppca1$pop3=c('Les Cottes Z4-1514','Goyet Q56-1','Mezmaiskaya1','Mezmaiskaya2','Vindija 87','Spy 94a')
  ppca1$method='ppca'


  adm1=adm[,1:5]
  colnames(adm1)=c('pop1','pop2','pop3','mu','std')
  adm1$method="admixtools"

  vec=full_join(ppca1, adm1)

  f3_plot=ggplot(vec, aes(x=mu, y=pop3, color=method)) +
    geom_point() +
    geom_errorbar(aes(xmin=mu-2*std, xmax=mu+2*std), width=0.1) +
    scale_color_discrete(labels=c('ADMIXTOOLS 2', 'PPCA')) +
    theme_classic() + ylab('Poplation X') + xlab("F3(Altai ; Vindija33.19 , X)") +
    guides(color=guide_legend(title="Method")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title=element_text(size=16))

  ggsave(plotf, f3_plot,
         width = 10, height = 3, dpi = 300, units = "in", device='png')
}

plotf3(admf=snakemake@input[["admf"]],
         ppca_muf=snakemake@input[["ppca_muf"]],
         ppca_stdf=snakemake@input[["ppca_stdf"]],
         plotf=snakemake@output[["plotf"]])
