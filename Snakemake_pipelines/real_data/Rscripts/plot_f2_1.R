library(ggplot2)

#flist=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/noisy_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#"/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/pca_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#"/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/ppca_direct_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#"/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/adjusted_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#"/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/admixtools2Norm_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv")
#fscale2="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/avgtables_npop4/ppca_direct_val_scale20_p1pop1_p2pop2_p3pop3_p4pop4.csv"

#slist=c("/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/noisy_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#          "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/pca_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#          "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/ppca_direct_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#          "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/adjusted_val_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv",
#          "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/admixtools2Norm_scale3_p1pop1_p2pop2_p3pop3_p4pop4.csv")
#s_scale2="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/all_pops_f2/simfiles/stdtables_npop4/ppca_direct_val_scale20_p1pop1_p2pop2_p3pop3_p4pop4.csv"


avg <-function(noisy,pca,ppca_direct,ppca_s2,adjusted,admix){
  
  
  xx=noisy
  tab_n=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  xx=pca
  tab_p=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  xx=ppca_direct
  tab_pd=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  xx=ppca_s2
  tab_p2=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  xx=adjusted
  tab_adj=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  xx=admix
  tab_adm=c(xx[1,2] , xx[3,4], (xx[1,3]+xx[1,4]+xx[2,3]+xx[2,4])/4)
  
  allt=as.data.frame(rbind(tab_n,tab_p,tab_pd,tab_p2,tab_adj,tab_adm))
  
  allt$V4 <- rownames(allt)
  return(allt)
}



plot_f2s <- function(flist,fscale2,slist,s_scale2){
  

  noisy=read.csv(file=flist[1], header=F, sep=',')
  pca=read.csv(file=flist[2], header=F, sep=',')
  ppca_direct=read.csv(file=flist[3], header=F, sep=',')
  ppca_s2=read.csv(file=fscale2, header=F, sep=',')
  adjusted=read.csv(file=flist[4], header=F, sep=',')
  admix=read.csv(file=flist[5], header=F, sep=',')
  
  
  std_noisy=read.csv(file=slist[1], header=F, sep=',')
  std_pca=read.csv(file=slist[2], header=F, sep=',')
  std_ppca_direct=read.csv(file=slist[3], header=F, sep=',')
  std_ppca_s2=read.csv(file=s_scale2, header=F, sep=',')
  std_adjusted=read.csv(file=slist[4], header=F, sep=',')
  std_admix=read.csv(file=slist[5], header=F, sep=',')
  
  
  
  meant=avg(noisy=noisy, pca=pca, ppca_direct=ppca_direct , ppca_s2=ppca_s2, adjusted=adjusted, admix=admix)
  stdt=avg(noisy=std_noisy, pca=std_pca, ppca_direct=std_ppca_direct , ppca_s2=std_ppca_s2, adjusted=std_adjusted, admix=std_admix)
  
  
  gg_mean <- meant %>% pivot_longer(
    cols = c("V1","V2","V3")
  )
  colnames(gg_mean) = c("method","split_time","mean_f2")
  
  
  gg_std <- stdt %>% pivot_longer(
    cols = c("V1","V2","V3")
  )
  colnames(gg_std) = c("method","split_time","std_error")
  
  
  gg=merge(gg_mean, gg_std, by=c("method","split_time"))
  
  gg[gg$method=="tab_n","method"] = "noisy"
  gg[gg$method=="tab_p","method"] = "pca"
  gg[gg$method=="tab_pd","method"] = "ppca_scale3"
  gg[gg$method=="tab_p2","method"] = "ppca_scale20"
  gg[gg$method=="tab_adj","method"] = "adjusted"
  gg[gg$method=="tab_adm","method"] = "admixtools"
  
  gg[gg$split_time=="V1","split_time"] = "1000"
  gg[gg$split_time=="V2","split_time"] = "2000"
  gg[gg$split_time=="V3","split_time"] = "5000"
  
   
  
  ggplot(gg, aes(x=split_time, y=mean_f2, fill=method)) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_f2-std_error, ymax=mean_f2+std_error), width=.2,
                  position=position_dodge(.9)) 
    
  ggsave(outplot,
         width = 8, height = 5, dpi = 150, units = "in", device='png')


}

plot_f2s(flist=snakemake@input[["flist"]], 
         fscale2=snakemake@input[["fscale2"]],
         slist=snakemake@input[["slist"]],
         s_scale2==snakemake@input[["s_scale2"]],
         f2plot==snakemake@output[["f2plot"]])
