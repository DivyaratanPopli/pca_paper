library(tidyverse)
library(ggplot2)


f_avg=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_PCA1_val_scale8",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_PCA1_val_scale12",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_ppca_direct_val_scale8",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_ppca_direct_val_scale12",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_admixtools2Norm_scale8")

f_std=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_PCA1_val_scale8",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_PCA1_val_scale12",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_ppca_direct_val_scale8",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_ppca_direct_val_scale12",
       "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_admixtools2Norm_scale8")


ftrue="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_true_val"
       


plot_f2s <- function(f_avg, f_std, ftrue, fplot){
  
  df = data.frame("Method"=factor(),"std"=numeric(),"mean_F"=numeric(), "Fstatistic"=numeric())
  
  
  for (i in seq(1,length(f_avg))){
    
    xx_avg = read.csv(file=f_avg[i], header=F, sep=',')
    xx_std = read.csv(file=f_std[i], header=F, sep=',')
    
    meth = strsplit(f_avg[i],'fstats_avg_')[[1]][2]

    
    xx1 = data.frame("Method"=meth, "std"=xx_std$V1[1], "mean_F"=xx_avg$V1[1], "Fstatistic"='1;3,4')
    xx2 = data.frame("Method"=meth, "std"=xx_std$V1[2], "mean_F"=xx_avg$V1[2], "Fstatistic"='4;1,2')
    xx3 = data.frame("Method"=meth, "std"=xx_std$V1[3], "mean_F"=xx_avg$V1[3], "Fstatistic"='1,3;2,4')
    
    df=rbind(df,xx1,xx2,xx3)
    
  }
  
  df1= df[df$Fstatistic=="1;3,4" | df$Fstatistic=="4;1,2",]
  df2= df[df$Fstatistic=="1,3;2,4",]
  
  
  true = read.csv(file=ftrue, header=F, sep=',')
  true1=true$V1[1]
  true2=true$V1[2]
  true3=true$V1[3]

  #gg$mean_f2 = abs(gg$mean_f2)
  fac=c("1;3,4","4;1,2","1,3;2,4")
  df$Fstatistic=factor(df$Fstatistic,levels=fac)
  
  ggplot(df1, aes(x=Fstatistic, y=mean_F, fill=Method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_F-(2*std), ymax=mean_F+(2*std)), width=.2,
                  position=position_dodge(.9)) +
    geom_segment(aes(x = 0.52, y = true1, xend = 1.5, yend = true1, linetype="True F-statistic"), colour= 'red') +
    geom_segment(aes(x = 1.52, y = true2, xend = 2.5, yend = true2, linetype="True F-statistic"), colour= 'red') +
    
    scale_linetype_manual(name ="", values = 2) +
    xlab("Pair of individuals") + ylab("F3 estimates") + scale_x_discrete(labels=c('1;3,4', '4;1,2')) +
    theme_classic()
  
  ggplot(df2, aes(x=Fstatistic, y=mean_F, fill=Method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_F-(2*std), ymax=mean_F+(2*std)), width=.2,
                  position=position_dodge(.9)) +
    geom_segment(aes(x = 0.52, y = true3, xend = 1.5, yend = true3, linetype="True F-statistic"), colour= 'red') +
    
    scale_linetype_manual(name ="", values = 2) +
    xlab("Pair of individuals") + ylab("F4 estimates") + scale_x_discrete(labels=c('1,3;2,4')) +
    theme_classic()
  
  
  #geom_hline(yintercept=c(true1,true2,true3,true3m), linetype="dashed", color = "black") + theme_classic()
  #coord_cartesian(ylim = c(0, 2000))
  
  ggsave(f2plot,
         width = 8, height = 5, dpi = 150, units = "in", device='png')
  
  
}

plot_f2s(flist=snakemake@input[["flist"]],
         fscale2=snakemake@input[["fscale2"]],
         pscale2=snakemake@input[["pscale2"]],
         slist=snakemake@input[["slist"]],
         s_scale2=snakemake@input[["s_scale2"]],
         p_scale2=snakemake@input[["p_scale2"]],
         f2plot=snakemake@output[["f2plot"]],
         pc1=snakemake@wildcards[["npcs"]],
         pc2=snakemake@wildcards[["npcs2"]], sp=snakemake@wildcards[["sp"]],
         truef=snakemake@input[["true"]])
