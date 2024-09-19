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


ppca_direct <- function(data, n_pcs=4, n_iter=1000, tol=1e-4, verbose=T){
  n_dims = ncol(data)
  #eq 31
  mu = rowMeans(data,na.rm=TRUE)
  centered_data = data - mu
  
  S = t(centered_data) %*% centered_data / nrow(data)
  SVD = svd(S)
  
  #eq 8
  sigma = sum(SVD$d[(n_pcs+1):n_dims]) / (n_dims - n_pcs)
  
  #eq 7
  W = SVD$u[,1:n_pcs] %*% sqrt(diag(SVD$d[1:n_pcs] - sigma))
  
  #reconstructing data and score mat
  
  T1= W %*% solve(t(W) %*% W)
  M= t(W) %*% W + (sigma * diag(n_pcs))
  xntn=solve(M) %*% t(W) %*% t(data)
  
  T_rec= t(T1 %*% M %*% xntn)
  
  scor=t(xntn)
  
  e = svd(W)
  pcs = e$u %*% diag(e$d)
  pcs = cbind(pcs, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))
  W = cbind(W, matrix(0, ncol=n_dims-n_pcs, nrow=n_dims))
  
  
  return(list(W=W, sigma=sigma, pcs=pcs, eval=e$d, evec=e$u,data_rec=T_rec, scores=scor))
}

run_pca <- function(data){
  z = prcomp(data)
  pcs = unname(t(z$rotation) * z$sdev)
  scores = z$x
  data_rec = z$x %*% (pcs/z$sdev)
  
  return(list(pcs = pcs, scores = scores, data_rec = data_rec))
}



genof="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/eigen.miss.geno_pc"
geno=read.csv(file=genof, header=F, sep=',')

indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/eigen.miss.ind"
ind=read.csv(file=indf, header=F, sep='\t')


ind[ind$V1=="pop4_1",1] = "anc4_pop4"
ind[ind$V1=="pop3_1",1] = "anc3_pop3"
ind[ind$V1=="pop2_1",1] = "anc2_pop2"
ind[ind$V1=="pop1_1",1] = "anc1_pop1"

ind$pop=unlist(lapply(strsplit(as.character(ind$V1), "_"), '[[', 1))

index1=ind$pop!="anc4" & ind$pop!="anc3" & ind$pop!="anc2" & ind$pop!="anc1"
index2=ind$pop=="anc4" | ind$pop=="anc3" | ind$pop=="anc2" | ind$pop=="anc1"
geno1=geno[,index1]
ind1=ind[index1,]

geno2=geno[,index2]
ind2=ind[index2,]

genom=as.matrix(geno2)
size <- nrow(genom)*ncol(genom)
indm <- sample(1:size, floor(0*size))
genom[indm] <- NA
as.data.frame(genom)


genom1=as.matrix(geno1)
size1 <- nrow(genom1)*ncol(genom1)
indm1 <- sample(1:size1, floor(0*size1))
genom1[indm1] <- NA
as.data.frame(genom1)


###ppca on the whole dataset

geno_all=cbind(genom1,genom)
ind_all=rbind(ind1,ind2)


#write it down for smartpca
geno_w <-geno_all
geno_w[is.na(geno_w)] <-9
outgeno="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/eigen.smartpca.geno"
write.table(geno_w, file = outgeno, row.names=FALSE, col.names=FALSE, sep='')

ind_w=ind_all[,c(1,2,4)]
outind="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/eigen.smartpca.ind"
write.table(ind_w, file = outind, row.names=FALSE, col.names=FALSE, sep='\t',quote=FALSE)


#############################


geno_all=geno_all/2
mu = rowMeans(geno_all,na.rm=TRUE)
geno_norm = geno_all - mu


xxxx=ppca_miss(as.matrix(geno_norm), n_pcs=8)
pc=t(xxxx$pcs)
#ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/pcs_method_ppca_direct/pcs_npcs8.csv"
#pc=read.csv(file=ppcaf, header=F, sep=',')

#indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/eigen.miss.ind"
#ind=read.csv(file=indf, header=F, sep='\t')
#ind[ind$V1=="pop4_1",1] = "anc4_pop4"
#ind[ind$V1=="pop3_1",1] = "anc3_pop3"
#ind[ind$V1=="pop2_1",1] = "anc2_pop2"
#ind[ind$V1=="pop1_1",1] = "anc1_pop1"

df=as.data.frame(cbind(ind_all$V1, pc[1,]*(-1), pc[2,]))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("anc4" = "pink","anc3" = "pink","anc2" = "pink","anc1" = "pink","pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 


#plot with smart pca##############################################
ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/shrinkage/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/missing0/smartpca.evec"
pc_1=t(read.csv(file=ppcaf, header=F, sep=''))
pc_2= matrix(as.numeric(pc_1[2:(n_pcs+1),2:ncol(pc_1)]),n_pcs,ncol(pc_1)-1)
eval=as.numeric(pc_1[2:(n_pcs+1),1])
pc=pc_2 * sqrt(eval)

ind1=pc_1[nrow(pc_1),2:ncol(pc_1)]
#indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/missing_ppca/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
#inds=read.csv(file=indf, header=F, sep='\t')
#ind_miss1 = ind_all[ind_all$V3!="pop5_7",]
df=as.data.frame(cbind(ind1, pc[1,], pc[2,]*(-1)))
df$pop=ind1
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("anc4" = "pink","anc3" = "pink","anc2" = "pink","anc1" = "pink","pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  #geom_point(position = position_jitter(w = 0.02, h = 0.02)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 


#plot neandertals with smart pca##############################################
n_pcs=2
ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/nea.evec"
pc_1=t(read.csv(file=ppcaf, header=F, sep=''))
pc_2= matrix(as.numeric(pc_1[2:(n_pcs+1),2:ncol(pc_1)]),n_pcs,ncol(pc_1)-1)
eval=as.numeric(pc_1[2:(n_pcs+1),1])
pc=pc_2 * sqrt(eval)

ind1=pc_1[nrow(pc_1),2:ncol(pc_1)]
#indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/missing_ppca/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
#inds=read.csv(file=indf, header=F, sep='\t')
#ind_miss1 = ind_all[ind_all$V3!="pop5_7",]
df=as.data.frame(cbind(ind1, pc[1,], pc[2,]*(-1)))
df$pop=ind1
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  #geom_point(position = position_jitter(w = 0.02, h = 0.02)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") + theme_classic()


####ppca for the same

genof="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.geno_pc"
geno=read.csv(file=genof, header=F, sep=',')

indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.ind"
ind=read.csv(file=indf, header=F, sep='\t')

geno_all=geno/2
mu = rowMeans(geno_all,na.rm=TRUE)
geno_norm = geno_all - mu

n_pcs=2
xxxx=ppca_miss(as.matrix(geno_norm), n_pcs=n_pcs)
pc=t(xxxx$pcs)
pc=pc[1:n_pcs,]
df=as.data.frame(cbind(ind$V1, pc[1,]*(-1), pc[2,]*(-1)))
df$pop=df$V1
  #unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")


df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 


##f stats
##admixtools
library(admixtools)
library(tidyverse)

f2dir = "/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/admixtools2_fmat/"
f2s = admixtools::read_f2(f2dir)

pop1="Altai"
pop2="Vindija33.19"
pop3=c("Les_Cottes_L35MQ25","Goyet_L35MQ25","Mezmaiskaya1_L35MQ25","Mezmaiskaya2_L35MQ25","VindijaG1_L35MQ25","Spy_L35MQ25")

f3stats = qp3pop(f2s, pop1, pop2, pop3)

dataf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.geno_pc"

filter="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/n_blocks.csv"

data_un=read.csv(file=dataf, header=F, sep=',')
fi=read.csv(file=filter, header=F, sep=',')

rm_start = sum(fi[1:f-1,1])+1 
rm_end = sum(fi[1:f,1])
data = data_un[-seq(rm_start, rm_end),]

is.na(data) <- data == "9"
data=data/2
mu = rowMeans(data,na.rm=TRUE)
data = data - mu

#plot a ppca 
ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/filter3/pcs_ppca_miss/pcs_npcs2.csv"
pc=read.csv(file=ppcaf, header=F, sep=',')[1:2,]
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[1,]), t(pc[2,])))
df$pop=df$V1
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")


df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 



#plot f3's

admf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/admixtools2/f3stats.csv"
adm=read.csv(file=admf, header=F, sep=',')

ppca_muf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/mu.csv"
ppca_mu=read.csv(file=ppca_muf, header=F, sep=',')

ppca_stdf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/std.csv"
ppca_std=read.csv(file=ppca_stdf, header=F, sep=',')

ppca1 = as.data.frame(cbind(ppca_mu, ppca_std))
colnames(ppca1)=c('mu','std')
ppca1$pop1=rep('Altai',6)
ppca1$pop2=rep('Vindija33.19',6)
ppca1$pop3=c('Les_Cottes_L35MQ25','Goyet_L35MQ25','Mezmaiskaya1_L35MQ25','Mezmaiskaya2_L35MQ25','VindijaG1_L35MQ25','Spy_L35MQ25')
ppca1$method='ppca'


adm1=adm[,1:5]
colnames(adm1)=c('pop1','pop2','pop3','mu','std')
adm1$method="admixtools"

vec=full_join(ppca1, adm1) 

ggplot(vec, aes(x=mu, y=pop3, color=method)) +
  geom_point() +
  geom_errorbar(aes(xmin=mu-2*std, xmax=mu+2*std), width=0.1) +
  theme_classic() + ylab('Poplation X') + xlab("F3(Altai ; Vindija33.19 , X)")
###########################ppca

p_muf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/ppca_mu.csv"
p_mu=read.csv(file=p_muf, header=F, sep=',')

p_stdf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/ppca_miss_val_scale2/ppca_std.csv"
p_std=read.csv(file=p_stdf, header=F, sep=',')

pc=p_mu
std=p_std
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/neandertal_shrinkage/all_ind.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t((-1)*pc[1,]), t((-1)*pc[2,]),t(std[1,]),t(std[2,])))
df$pop=df$V1
colnames(df) = c("ind", "pc1", "pc2", "std1", "std2", "pop")

colors1 <- c("Les_Cottes_L35MQ25" = "dodgerblue2","Goyet_L35MQ25" = "green4", "Mezmaiskaya1_L35MQ25" = "#6A3D9A","Mezmaiskaya2_L35MQ25" = "#FF7F00", "VindijaG1_L35MQ25" = "gold1", "Spy_L35MQ25"= "skyblue2", "Altai" = "steelblue4", "Vindija33.19" = "maroon", "Denisova" = "darkorange1")
"dodgerblue2", "#E31A1C", # red
"green4",
"#6A3D9A", # purple
"#FF7F00", # orange
"black", "gold1",
"skyblue2", "#FB9A99", # lt pink
"palegreen2",
"#CAB2D6", # lt purple
"#FDBF6F", # lt orange
"gray70", "khaki2",
"maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
"darkturquoise", "green1", "yellow4", "yellow3",
"darkorange4", "brown"

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)
df$std1=as.numeric(df$std1)
df$std2=as.numeric(df$std2)

df$xmin= df$pc1 - df$std1
df$xmax= df$pc1 + df$std1
df$ymin= df$pc2 - df$std2
df$xmin= df$pc1 - df$std1

ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point() +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") +
  geom_errorbar(aes(xmin=pc1-(2*std1), xmax=pc1+(std1*2))) +
  geom_errorbar(aes(ymin=pc2-(2*std2), ymax=pc2+(std2*2))) + theme_classic()


#Alba's neandertals

#plot neandertals with smart pca##############################################
library(ggrepel)

n_pcs=2
ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/alba_sub_Chag_Vin_Goyet/nea.evec"
pc_1=t(read.csv(file=ppcaf, header=F, sep=''))
pc_2= matrix(as.numeric(pc_1[2:(n_pcs+1),2:ncol(pc_1)]),n_pcs,ncol(pc_1)-1)
eval=as.numeric(pc_1[2:(n_pcs+1),1])
pc=pc_2 * sqrt(eval)

ind1=pc_1[nrow(pc_1),2:ncol(pc_1)]
#indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/missing_ppca/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
#inds=read.csv(file=indf, header=F, sep='\t')
#ind_miss1 = ind_all[ind_all$V3!="pop5_7",]
df=as.data.frame(cbind(ind1, pc[1,], pc[2,]*(-1)))
df$pop=ind1
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors2 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  #geom_point(position = position_jitter(w = 0.02, h = 0.02)) +
  geom_point() +
  geom_label_repel(aes(label = ind),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.overlaps=30) +
  #geom_point()
  #scale_color_manual(values = colors1, name="pop") + 
  theme_classic()# +xlim(c(-1.7,0.9))

#same for ppca

genof="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/alba_sub_Chag_Vin_Goyet/allsamples_Ref.sub.geno_pc"
geno=read.csv(file=genof, header=F, sep=',')

indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/alba_sub_Chag_Vin_Goyet/allsamples_Ref.sub.ind"
ind=read.csv(file=indf, header=F, sep='\t')

is.na(geno) <- geno == "9"
geno_all=geno/2
mu = rowMeans(geno_all,na.rm=TRUE)
geno_norm = geno_all - mu

n_pcs=2
xxxx=ppca_miss(as.matrix(geno_norm), n_pcs=n_pcs)
pc=t(xxxx$pcs)
pc=pc[1:n_pcs,]
df=as.data.frame(cbind(ind$V1, pc[1,]*(-1), pc[2,]*(-1)))
df$pop=df$V1
#unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("Les_Cottes_L35MQ25" = "pink","Goyet_L35MQ25" = "black", "Mezmaiskaya1_L35MQ25" = "grey","Mezmaiskaya2_L35MQ25" = "cyan", "VindijaG1_L35MQ25" = "blue", "Spy_L35MQ25"= "green", "Altai" = "yellow4", "Vindija33.19" = "red", "Denisova" = "darkorange1")


df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  #geom_point(position = position_jitter(w = 0.02, h = 0.02)) +
  geom_point() +
  geom_label_repel(aes(label = ind),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',max.overlaps=50) +
  #geom_point()
  #scale_color_manual(values = colors1, name="pop") + 
  theme_classic() +guides(color=guide_legend(ncol =1)) #+xlim(c(-1.7,0.9))





