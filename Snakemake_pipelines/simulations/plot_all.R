library(gridExtra)
library(ggplot2)



#f2f="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/npop10_nind100/plots_8_12/mu0.05_f2.csv"
#f3f="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/npop10_nind100/plots_8_12/mu0.05_f3.csv"
#f4f="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/npop10_nind100/plots_8_12/mu0.05_f4.csv"

#ftrue="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/fstats_avg_true_val"
#true_f2="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/avgrun/npop10_nind100/true_val/f2mat15"

plot_all <- function(f2f, f3f, f4f, ftrue, true_f2, plotf){

  gg=read.csv(file=f2f, header=T, sep=',')
  df1=read.csv(file=f3f, header=T, sep=',')
  df2=read.csv(file=f4f, header=T, sep=',')

  true = read.csv(file=ftrue, header=F, sep=',')
  true1=true$V1[1]
  true2=true$V1[2]
  true3=true$V1[3]

  truef2=read.csv(file=true_f2, header=F, sep=',')
  tr1=truef2[1,2]
  tr2=truef2[3,4]
  tr3m=truef2[2,3]
  tr3=truef2[1,4]
  
  gg$split_time=as.factor(gg$split_time)
  
  g_f2=ggplot(gg, aes(x=split_time, y=mean_f2, fill=method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_f2-(2*std_error), ymax=mean_f2+(2*std_error)), width=.2,
                  position=position_dodge(.9)) +
    geom_segment(aes(x = 0.52, y = tr1, xend = 1.5, yend = tr1, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 1.52, y = tr2, xend = 2.5, yend = tr2, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 2.52, y = tr3m, xend = 3.5, yend = tr3m, linetype="True F2"), colour= 'red') +
    geom_segment(aes(x = 3.52, y = tr3, xend = 4.5, yend = tr3, linetype="True F2"), colour= 'red') +
    scale_linetype_manual(name ="", values = 2) + scale_x_discrete(labels=c(expression(F[2](X[1],X[2])), expression(F[2](X[3],X[4])),expression(F[2](X[2],X[3])),expression(F[2](X[1],X[4])))) +
    xlab("Statistic") + ylab("Estimates") + 
    theme_classic()


  g_f3=ggplot(df1, aes(x=Fstatistic, y=mean_F, fill=Method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_F-(2*std), ymax=mean_F+(2*std)), width=.2,
                  position=position_dodge(.9)) +
    geom_segment(aes(x = 0.52, y = true1, xend = 1.5, yend = true1, linetype="True F-statistic"), colour= 'red') +
    geom_segment(aes(x = 1.52, y = true2, xend = 2.5, yend = true2, linetype="True F-statistic"), colour= 'red') +

    scale_linetype_manual(name ="", values = 2) + scale_x_discrete(labels=c(expression(F[3](X[1]*";"*X[3],X[4])), expression(F[3](X[4]*";"*X[1],X[2])))) +
    xlab("Statistic") + ylab("Estimates") +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank())

  g_f4=ggplot(df2, aes(x=Fstatistic, y=mean_F, fill=Method)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean_F-(2*std), ymax=mean_F+(2*std)), width=.2,
                  position=position_dodge(.9)) +
    geom_segment(aes(x = 0.52, y = true3, xend = 1.5, yend = true3, linetype="True F-statistic"), colour= 'red') +

    scale_linetype_manual(name ="", values = 2) +
    xlab("Statistic") + ylab("Estimates") + scale_x_discrete(labels=c(expression(F[4](X[1],X[3]*";"*X[2],X[4])))) +
    theme_classic() +
    theme(legend.position = "none",
          legend.title = element_blank())



  combined_plot = grid.arrange(
    g_f2,
    arrangeGrob(g_f3, g_f4, ncol = 2, widths = c(1, 0.5)),  # Adjust the widths here
    heights = c(2, 1)
  )

  ggsave(plotf, combined_plot,
         width = 8, height = 5, dpi = 300, units = "in", device='png')
}


plot_all(f2f=snakemake@input[["f2f"]],
         f3f=snakemake@input[["f3f"]],
         f4f=snakemake@input[["f4f"]],
         ftrue=snakemake@input[["ftrue"]],
         true_f2=snakemake@input[["true_f2"]],
         plotf=snakemake@output[["plotf"]])
