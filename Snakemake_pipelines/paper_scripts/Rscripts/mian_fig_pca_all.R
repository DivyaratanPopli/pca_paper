library(ggplot2)

x1="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/npop10_nind100/mu0.05_f2_fig_all.csv"
x2="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/npop10_nind100/mu0.05_f2_fig_all.csv"
x3="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind_missing/simfiles/Ne1000/split_times1000/npop10_nind100/missing0.5/mu0.05_f2_fig_all.csv"
truef="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/avgrun/npop10_nind100/true_val/f2mat8"
nfile_avg1="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/avgAccuracy_noisy_val_scale8"
nfile_avg2="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/avgAccuracy_noisy_val_scale8"
nfile_avg3="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Fig_comparison_PCA_PPCA_LSE_1ind_missing/simfiles/Ne1000/split_times1000/mu0.05/average_run/npop10_nind100/missing0.5/avgAccuracy_noisy_val_scale8"


df1=read.csv(file=x1, header=T, sep=',')
df2=read.csv(file=x2, header=T, sep=',')
df3=read.csv(file=x3, header=T, sep=',')

tr=read.csv(file=truef, header=F, sep=',')
noisy1=read.csv(file=nfile_avg1, header=F, sep=',')
noisy2=read.csv(file=nfile_avg2, header=F, sep=',')
noisy3=read.csv(file=nfile_avg3, header=F, sep=',')



tru=data.frame("Pair"="1_4", "true_val"=tr[1,4])

noi1=data.frame("Pair"="1_4", "noisy_val"=noisy1[1,4])
noi2=data.frame("Pair"="1_4", "noisy_val"=noisy2[1,4])
noi3=data.frame("Pair"="1_4", "noisy_val"=noisy3[1,4])
noi=rbind(noi1,noi2,noi3)
noi$data=c("Population-based", "Individual-based","With missing data")

df=rbind(df1,df2,df3)

df[df$Method=="emu","Method"] = "PCA"
df[df$Method=="ppca_miss","Method"] = "PPCA"
df[df$data=="Indiviadual-based-missing","data"] = "With missing data"

ggplot(df, aes(x=Scale, y=F2, col=Method)) +
 geom_point(alpha=0.5) + facet_grid(~factor(data,levels=c("Population-based", "Individual-based","With missing data"))) +
  geom_errorbar(aes(ymin=F2-(2*SE), ymax=F2+(2*SE)), width=.2, alpha=0.8) +
  scale_color_manual(values=c(PCA="#E69F00", LSE="#009E73", PPCA="#CC79A7")) +
  geom_hline(data = tru, aes(yintercept = true_val, linetype="True F2"), alpha=0.2) +
  geom_hline(data = noi, aes(yintercept = noisy_val, linetype="Uncorrected f2"), alpha=0.5) +
  scale_linetype_manual(name = "", values = c(2,3),
                        guide = guide_legend(override.aes = list(linetype = "dashed"))) +
  xlab("Number of PC's used") + ylab("f2") +
  theme_bw() 

ggsave(f2plot, xx, 
       width = 8, height = 5, dpi = 150, units = "in", device='png')
