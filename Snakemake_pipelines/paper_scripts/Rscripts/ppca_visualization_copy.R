library(dplyr)
library(ggplot2)
library(magrittr)

# load slendr itself

library(slendr)
init_env()

sp=100
Ne=1000
gentime=30
mgen_array = c(1.5, 0.4, 0.1, 2.4, 0.2)
madm_array = c(0.8, 0.5, 0.4, 0.6, 0.6)
nanc_array = c(0.001, 0.001, 0.001, 0.001)
st_array = c(1,2,5,5.5)
mutation_rate=1e-8
recombination_rate=1e-8
slength=1e8
mgen=0.5
mprop=0
s1=10
genfile="test.geno"

sim_func <- function(sp, Ne, gentime, st_array, mgen_array, madm_array, nanc_array, mgen, mprop,
                     mutation_rate, recombination_rate, slength, s1,
                     genfile){
  
  
  
  mgen=mgen*sp*gentime    #gentime is multiplied to convert to years.
  
  prefix = strsplit(genfile, split=".geno")[[1]]
  
  mgen_array=mgen_array*sp*gentime   #
  mgenF = mgen_array[1]
  mgenG = mgen_array[2]
  mgenH = mgen_array[3]
  mgenI = mgen_array[4]
  mgenJ = mgen_array[5]
  
  
  madmF = madm_array[1]
  madmG = madm_array[2]
  madmH = madm_array[3]
  madmI = madm_array[4]
  madmJ = madm_array[5]
  
  st_array=st_array*sp*gentime
  st01 = st_array[1]
  st02 = st_array[2]
  st03 = st_array[3]
  st04 = st_array[4]
  
  nanc_array=nanc_array*sp*gentime
  
  n1anc = nanc_array[1]
  n2anc = nanc_array[2]
  n3anc = nanc_array[3]
  n4anc = nanc_array[4]
  
  #creating first 5 pops
  pop5 <- population("pop5", N = Ne, time = st04+1)
  pop1 <- population("pop1", N = Ne, time = st04, parent = pop5)
  pop4 <- population("pop4", N = Ne, time = st03, parent = pop1)
  pop3 <- population("pop3", N = Ne, time = st02, parent = pop4)
  pop2 <- population("pop2", N = Ne, time = st01, parent = pop1)
  
  #admixing first five populations to get five admixed populations
  
  pop6 <- population("pop6", N = Ne, time = mgenF, parent = pop4)
  pop7 <- population("pop7", N = Ne, time = mgenG, parent = pop3)
  pop8 <- population("pop8", N = Ne, time = mgenH, parent = pop2)
  pop9 <- population("pop9", N = Ne, time = mgenI, parent = pop1)
  pop10 <- population("pop10", N = Ne, time = mgenJ, parent = pop1)
  
  gf <- list(
    gene_flow(from = pop5, to = pop6, start = mgenF, end = mgenF-gentime, rate = 1-madmF),
    gene_flow(from = pop4, to = pop7, start = mgenG, end = mgenG-gentime, rate = 1-madmG),
    gene_flow(from = pop3, to = pop8, start = mgenH, end = mgenH-gentime, rate = 1-madmH),
    gene_flow(from = pop4, to = pop9, start = mgenI, end = mgenI-gentime, rate = 1-madmI),
    gene_flow(from = pop2, to = pop10, start = mgenJ, end = mgenJ-gentime, rate = 1-madmJ),
    
    gene_flow(from = pop3, to = pop2, start = mgen, end = mgen-gentime, rate = mprop)
  )
  
  model <- compile_model(
    populations = list(pop5,pop1,pop4,pop3,pop2,pop6,pop7,pop8,pop9,pop10),
    generation_time = gentime,
    direction = "backward",
    gene_flow = gf
    
    
  )
  plot_model(model)
  
  s_present <- schedule_sampling(model, times = 0, list(pop1, s1),list(pop2, s1),list(pop3, s1),list(pop4, s1),list(pop5, s1),list(pop6, s1),list(pop7, s1),list(pop8, s1),list(pop9, s1),list(pop10, s1))
  s_ancient1 <- schedule_sampling(model, times = n1anc, list(pop1, 1))
  s_ancient2 <- schedule_sampling(model, times = n2anc, list(pop2, 1))
  s_ancient3 <- schedule_sampling(model, times = n3anc, list(pop3, 1))
  s_ancient4 <- schedule_sampling(model, times = n4anc, list(pop4, 1))
  samples <- rbind(s_present,s_ancient1,s_ancient2,s_ancient3,s_ancient4)
  
  ts <- msprime(model, samples=samples, sequence_length = slength, recombination_rate = recombination_rate) %>%
    ts_mutate(mutation_rate = mutation_rate)
  
  ts_eigenstrat(ts, prefix = prefix)
  
  
  
}




#ppca visualization

ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/pcs_method_ppca_direct/pcs_npcs10.csv"
pc=read.csv(file=ppcaf, header=F, sep=',')
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[5,]), t(pc[6,])))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point(position = position_jitter(w = 0.01, h = 0.01)) +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 



#with sp=100 

ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times100/mu0/run1/npop10_nind100/pcs_method_ppca_direct/pcs_npcs10.csv"
pc=read.csv(file=ppcaf, header=F, sep=',')
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times100/mu0/run1/npop10_nind100/eigen.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[5,]), t(pc[6,])))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  #geom_point(pprint(sp, Ne, gentime, st_array, mgen_array, madm_array, nanc_array, mgen, mprop,
  mutation_rate, recombination_rate, slength, s1,
genfile, outfile, outf2)osition = position_jitter(w = 0.02, h = 0.02)) +
  geom_point() +
  scale_color_manual(values = colors1, name="pop") 

cat(sp, Ne, gentime, st_array, mgen_array, madm_array, nanc_array, mgen, mprop,
    mutation_rate, recombination_rate, slength, s1,
    genfile, outfile, outf2)                



######################




ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Ne5000sp1000/simfiles/Ne5000/split_times1000/mu0/run1/npop10_nind100/pcs_method_ppca_direct/pcs_npcs10.csv"
pc=read.csv(file=ppcaf, header=F, sep=',')
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/Ne5000sp1000/simfiles/Ne5000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[1,]), t(pc[2,])))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point(position = position_jitter(w = 0.03, h = 0.03)) +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 




###########################


#ppca visualization

ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/sp200Ne1000/simfiles/Ne1000/split_times100/mu0/run1/npop10_nind100/pcs_method_ppca_direct/pcs_npcs10.csv"
pc=read.csv(file=ppcaf, header=F, sep=',')
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/sp200Ne1000/simfiles/Ne1000/split_times100/mu0/run1/npop10_nind100/eigen.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[1,]), t(pc[2,])))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point(position = position_jitter(w = 0.01, h = 0.01)) +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 


#ppca visualization with bjk

ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/partial_resampling/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/block1/method_ppca_direct/pcs_npcs10.csv"
ppcaf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/partial_resampling/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/pcs_method_ppca_direct/pcs_npcs10.csv"

pc=read.csv(file=ppcaf, header=F, sep=',')
indf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/block_resampling/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.ind"
inds=read.csv(file=indf, header=F, sep='\t')

df=as.data.frame(cbind(inds$V1, t(pc[5,]), t(pc[6,])))
df$pop=unlist(lapply(strsplit(as.character(df$V1), "_"), '[[', 1))
colnames(df) = c("ind", "pc1", "pc2", "pop")

colors1 <- c("pop5" = "black", "pop1" = "red", "pop2" = "blue", "pop3"= "green", "pop4" = "yellow2", "pop6" = "yellow4", "pop7" = "yellowgreen", "pop8" = "cyan3", "pop9" = "darkorange1", "pop10" = "purple4")

df$pc1=as.numeric(df$pc1)
df$pc2=as.numeric(df$pc2)


ggplot(df, aes(x=pc1, y=pc2, color=pop)) +
  geom_point(position = position_jitter(w = 0, h = 0)) +
  #geom_point()
  scale_color_manual(values = colors1, name="pop") 






