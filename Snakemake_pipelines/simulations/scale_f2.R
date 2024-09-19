library(tidyverse)
library(ggplot2)

#flist=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_ppca_direct_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_pca_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_PCA1_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_ppca_direct_val_scale20",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_pca_val_scale20",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/avgAccuracy_PCA1_val_scale20")

#slist=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_ppca_direct_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_pca_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_PCA1_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_ppca_direct_val_scale20",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_pca_val_scale20",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/average_run/npop10_nind100/stdDev_PCA1_val_scale20")

#truef="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0/avgrun/npop10_nind100/true_val/f2mat8"



plot_f2s <- function(flist,slist,f2plot, pc1, pc2, sp, truef){

  df1 <- data.frame(matrix(ncol = 5, nrow = 0))
  x <- c("Method", "Pair", "Scale", "F2", "SE")
  colnames(df1) <- x

  tr=read.csv(file=truef, header=F, sep=',')

  for(i in seq(1:length(flist))){
    print(i)
    xx=read.csv(file=flist[i], header=F, sep=',')
    ss=read.csv(file=slist[i], header=F, sep=',')
    temp=strsplit(flist[i], split="Accuracy_")[[1]][2]
    Method=rep(strsplit(temp, split="_val_")[[1]][1], 4)
    Scale=rep(as.integer(strsplit(strsplit(temp, split="_val_")[[1]][2], split="scale")[[1]][2]), 4)
    Pair=c("1_2", "3_4", "2_3", "1_4")
    F2=c(xx[1,2],xx[3,4],xx[2,3],xx[1,4])
    SE=c(ss[1,2],ss[3,4],ss[2,3],ss[1,4])

    df2=data.frame(Method,Pair,Scale,F2,SE)
    df1=rbind(df1,df2)
  }


  df1$Pair_f = factor(df1$Pair, levels=c('1_2','3_4','2_3','1_4'))
  df1[df1$Method == "pca","Method"]="PCA"
  df1[df1$Method == "PCA1","Method"]="LSE"
  df1[df1$Method == "ppca_direct","Method"]="PPCA"

  tru=data.frame("Pair_f"=c("1_2", "3_4", "2_3", "1_4"), "true_val"=c(tr[1,2],tr[3,4],tr[2,3],tr[1,4]))
  tru$Pair_f = factor(tru$Pair_f, levels=c('1_2','3_4','2_3','1_4'))

  ggplot(df1, aes(x=Scale, y=F2, col=Method), alpha=0.1) +
    geom_point() +
    facet_grid(cols=vars(Pair_f)) +
    geom_errorbar(aes(ymin=F2-(2*SE), ymax=F2+(2*SE)), width=.2) +
    scale_color_manual(values=c(PCA="#999999", LSE="#E69F00", PPCA="#56B4E9")) +
    geom_hline(data = tru, aes(yintercept = true_val, linetype="True F2")) +
    scale_linetype_manual(name = "", values = 2,
              guide = guide_legend(override.aes = list(linetype = "dashed"))) +
    xlab("Number of PC's used") + ylab("f2") +
    theme_bw()

  ggsave(f2plot,
         width = 8, height = 5, dpi = 150, units = "in", device='png')


}

plot_f2s(flist=snakemake@input[["flist"]],
         slist=snakemake@input[["slist"]],
         f2plot=snakemake@output[["f2plot"]],
         pc1=snakemake@wildcards[["npcs"]],
         pc2=snakemake@wildcards[["npcs2"]], sp=snakemake@wildcards[["sp"]],
         truef=snakemake@input[["true"]])
