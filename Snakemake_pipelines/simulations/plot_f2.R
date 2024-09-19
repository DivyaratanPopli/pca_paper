library(tidyverse)
library(ggplot2)

#flist=c("simfiles/mu0/average_run/npop4_nind400/missing0/subsetInds40/p1pop1_p2pop2_p3pop3_p4pop4/avgAccuracy_ppca_direct_val_scale20",
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


avg <-function(pca,ppca_direct,admix,ppca_s2,pca_s2){


  xx=pca
  tab_p=c(xx[1,2] , xx[3,4], xx[2,3], xx[1,4])

  xx=ppca_direct
  tab_pp=c(xx[1,2] , xx[3,4], xx[2,3], xx[1,4])

  xx=admix
  tab_adm=c(xx[1,2] , xx[3,4], xx[2,3], xx[1,4])

  xx=ppca_s2
  tab_pp2=c(xx[1,2] , xx[3,4], xx[2,3], xx[1,4])

  xx=pca_s2
  tab_p2=c(xx[1,2] , xx[3,4], xx[2,3], xx[1,4])





  allt=as.data.frame(rbind(tab_p,tab_pp,tab_adm,tab_pp2,tab_p2))

  allt$V5 <- rownames(allt)
  return(allt)
}



plot_f2s <- function(flist,fscale2,pscale2,slist,s_scale2,p_scale2,f2plot,pc1,pc2, sp,truef){


  pca=read.csv(file=flist[1], header=F, sep=',')
  ppca=read.csv(file=flist[2], header=F, sep=',')
  admix=read.csv(file=flist[3], header=F, sep=',')
  ppca_s2=read.csv(file=fscale2, header=F, sep=',')
  pca_s2=read.csv(file=pscale2, header=F, sep=',')
  true=read.csv(file=truef, header=F, sep=',')

  std_pca=read.csv(file=slist[1], header=F, sep=',')
  std_ppca=read.csv(file=slist[2], header=F, sep=',')
  std_admix=read.csv(file=slist[3], header=F, sep=',')
  std_ppca_s2=read.csv(file=s_scale2, header=F, sep=',')
  std_pca_s2=read.csv(file=p_scale2, header=F, sep=',')



  meant=avg(pca=pca, ppca_direct=ppca, admix=admix, ppca_s2=ppca_s2, pca_s2=pca_s2)
  stdt=avg(pca=std_pca, ppca_direct=std_ppca, admix=std_admix, ppca_s2=std_ppca_s2, pca_s2=std_pca_s2)

  true1=true[1,2]
  true2=true[3,4]
  true3m=true[2,3]
  true3=true[1,4]

  gg_mean <- meant %>% pivot_longer(
    cols = c("V1","V2","V3","V4")
  )
  colnames(gg_mean) = c("method","split_time","mean_f2")


  gg_std <- stdt %>% pivot_longer(
    cols = c("V1","V2","V3","V4")
  )
  colnames(gg_std) = c("method","split_time","std_error")

  label1=paste("PPCA_scale", pc1, sep="")
  label2=paste("PPCA_scale", pc2, sep="")
  label3=paste("LSE_scale", pc1, sep="")
  label4=paste("LSE_scale", pc2, sep="")

  gg=merge(gg_mean, gg_std, by=c("method","split_time"))

  gg[gg$method=="tab_p","method"] = label3
  gg[gg$method=="tab_pp","method"] = label1
  gg[gg$method=="tab_adm","method"] = "ADMIXTOOLS 2"
  gg[gg$method=="tab_pp2","method"] = label2
  gg[gg$method=="tab_p2","method"] = label4


  gg[gg$split_time=="V1","split_time"] = 1
  gg[gg$split_time=="V2","split_time"] = 2
  gg[gg$split_time=="V3","split_time"] = 3
  gg[gg$split_time=="V4","split_time"] = 4

  #gg$mean_f2 = abs(gg$mean_f2)
  fac=c("Pop1_Pop2","Pop3_Pop4","Pop2_Pop3","Pop1_Pop4")
  #gg$split_time=factor(gg$split_time,levels=fac)

  ggplot(gg, aes(x=split_time, y=mean_f2, fill=method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_f2-(2*std_error), ymax=mean_f2+(2*std_error)), width=.2,
                  position=position_dodge(.9)) +
    #geom_hline(aes(yintercept= true1, linetype = "Pop1_Pop2"), colour= 'red') +
    geom_segment(aes(x = 0.52, y = true1, xend = 1.5, yend = true1, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 1.52, y = true2, xend = 2.5, yend = true2, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 2.52, y = true3m, xend = 3.5, yend = true3m, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 3.52, y = true3, xend = 4.5, yend = true3, linetype="True F2"), colour= 'red') +
    scale_linetype_manual(name ="", values = 2) +
    xlab("Pair of individuals") + ylab("F2 estimates") + scale_x_discrete(labels=c('Pop1_Pop2', 'Pop3_Pop4', 'Pop2_Pop3', 'Pop1_Pop4')) +
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
