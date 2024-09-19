genof="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times1000/mu0/run1/npop10_nind100/eigen.geno_pc"
data=read.csv(file=genof, header=F, sep=',')
acf(data[,1])
dim(data)
mean(ess(data))

