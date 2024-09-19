library(ggplot2)


ppca_miss <- function(data, n_pcs=4, n_iter=100, tol=1e-5, verbose=T){

  data_orig = data
  NAs = is.na(data_orig)
  data[NAs] = 0
  W=0
  i=0
  while(i<n_iter){

    res = ppca_direct(data=as.matrix(data), n_pcs=n_pcs)
    data[NAs] = res$data_rec[NAs]

    diff=sum((res$W - W)^2)
    print(diff)
    W=res$W
    i=i+1

    if(diff < tol){
      return(list(W=res$W, sigma=res$sigma, pcs=res$pcs, eval=res$eval,
                  evec=res$evec,data_rec=res$data_rec, scores=res$scores))
    }
  }
  return(list(W=res$W, sigma=res$sigma, pcs=res$pcs, eval=res$eval,
              evec=res$evec,data_rec=res$data_rec, scores=res$scores))

}

plotf <- function(n_pcs, ppcaf, ppcaf1, indf, outplot1, outplot2){

  #ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/nea.evec"


  pc_1=t(read.csv(file=ppcaf, header=F, sep=''))

  pc_2= matrix(as.numeric(pc_1[2:(n_pcs+1),2:ncol(pc_1)]),n_pcs,ncol(pc_1)-1)
  eval=as.numeric(pc_1[2:(n_pcs+1),1])
  pc=pc_2 * sqrt(eval)

  ind1=pc_1[nrow(pc_1),2:ncol(pc_1)]

  df=as.data.frame(cbind(ind1, pc[1,], pc[2,]*(-1)))
  df$pop=ind1
  colnames(df) = c("ind", "pc1", "pc2", "pop")

  #colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet Q56-1" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")

  df[df$pop=="Les_Cottes_L35MQ25", "pop"] = "Les Cottes Z4-1514"
  df[df$pop=="Goyet_L35MQ25", "pop"] = "Goyet Q56-1"
  df[df$pop=="Mezmaiskaya1_L35MQ25", "pop"] = "Mezmaiskaya1"
  df[df$pop=="Mezmaiskaya2_L35MQ25", "pop"] = "Mezmaiskaya2"
  df[df$pop=="VindijaG1_L35MQ25", "pop"] = "Vindija 87"
  df[df$pop=="Spy_L35MQ25", "pop"] = "Spy 94a"


  colors1 <- c("Les Cottes Z4-1514" = "#88CCEE", "Goyet Q56-1" = "#CC6677", "Mezmaiskaya1" = "#DDCC77", "Mezmaiskaya2"= "#117733",
               "Vindija 87" = "#332288", "Spy 94a" = "#888888", "Vindija33.19" = "#44AA99", "Altai" = "#999933",
               "Denisova" = "#882255")

  df$pc1=as.numeric(df$pc1)
  df$pc2=as.numeric(df$pc2)


  p1 = ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
    #geom_point(position = position_jitter(w = 0.02, h = 0.02)) +
    geom_point() +
    scale_color_manual(values = colors1, name="pop") + theme_classic() + xlab("PC 1 (82.78%)")+ ylab("PC 2 (17.21%)") + coord_fixed(ratio = 1) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_hline(yintercept = 0, color = "grey60") +
    guides(color=guide_legend(title="")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title=element_text(size=16))


  ggsave(outplot1, p1,
         width = 8, height = 5, dpi = 300, units = "in", device='png')



  #ppcaf1="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/filter3/pcs_ppca_miss/pcs_npcs2.csv"
  pc=read.csv(file=ppcaf1, header=F, sep=',')[1:2,]
  #indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.ind"
  inds=read.csv(file=indf, header=F, sep='\t')

  df=as.data.frame(cbind(inds$V1, t(pc[1,]), t(pc[2,])))
  df$pop=df$V1
  colnames(df) = c("ind", "pc1", "pc2", "pop")

  #colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")


  df$pc1=as.numeric(df$pc1) *(-1)
  df$pc2=as.numeric(df$pc2) *(-1)

  df[df$pop=="Les_Cottes_L35MQ25", "pop"] = "Les Cottes Z4-1514"
  df[df$pop=="Goyet_L35MQ25", "pop"] = "Goyet Q56-1"
  df[df$pop=="Mezmaiskaya1_L35MQ25", "pop"] = "Mezmaiskaya1"
  df[df$pop=="Mezmaiskaya2_L35MQ25", "pop"] = "Mezmaiskaya2"
  df[df$pop=="VindijaG1_L35MQ25", "pop"] = "Vindija 87"
  df[df$pop=="Spy_L35MQ25", "pop"] = "Spy 94a"

  p2=ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
    geom_point() +
    scale_color_manual(values = colors1, name="pop") + theme_classic()+ xlab("PC 1 (71.59%)")+ ylab("PC 2 (28.40%)") + coord_fixed(ratio = 1) +
    geom_vline(xintercept = 0, color = "grey60") +
    geom_hline(yintercept = 0, color = "grey60") +
    guides(color=guide_legend(title="")) +
    theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16), legend.text=element_text(size=14), legend.title=element_text(size=16))

  ggsave(outplot2, p2,
         width = 8, height = 5, dpi = 300, units = "in", device='png')

}

plotf(n_pcs=snakemake@params[["n_pcs"]],
         ppcaf=snakemake@input[["ppcaf"]],
         ppcaf1=snakemake@input[["ppcaf1"]],
         indf=snakemake@input[["indf"]],
         outplot1=snakemake@output[["outplot1"]],
         outplot2=snakemake@output[["outplot2"]])
