library(tidyverse)
library(ggplot2)


#f_avg=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_PCA1_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_PCA1_val_scale12",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_ppca_direct_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_ppca_direct_val_scale12",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_admixtools2Norm_scale8")

#f_std=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_PCA1_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_PCA1_val_scale12",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_ppca_direct_val_scale8",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_ppca_direct_val_scale12",
#        "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_std_admixtools2Norm_scale8")


#ftrue="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_true_val"


plot_f2s <- function(f_avg, f_std, ftrue, f3out, f4out){
  
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
  
  write.table(df1, file = f3out, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(df2, file = f4out, row.names=FALSE, col.names=FALSE, sep=',')
  
  
}

plot_f2s(f_avg=snakemake@input[["f_avg"]],
                     f_std=snakemake@input[["f_std"]],
                     ftrue=snakemake@input[["f_true"]],
                     f3out=snakemake@output[["f3out"]],
                     f4out=snakemake@output[["f4out"]])
  
