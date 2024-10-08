library(dplyr)
library(ggplot2)
library(magrittr)

# load slendr itself
library(slendr)
init_env()

sim_func <- function(sp, Ne, gentime, st_array, mgen_array, madm_array, nanc_array, mgen, mprop,
                     mutation_rate, recombination_rate, slength, s1,
                     genfile, outfile, outf2){
  
  mutation_rate=as.numeric(mutation_rate)
  recombination_rate=as.numeric(recombination_rate)
  slength=as.numeric(slength)
  
  s1=as.integer(as.integer(s1)/10)
  Ne=as.integer(Ne)
  sp=as.integer(sp)
  mgen=mgen*sp*gentime    #gentime is multiplied bcoz slendr works with yrs. and not generation
  
  prefix = strsplit(genfile, split=".geno")[[1]]
  
  mgen_array=mgen_array*sp*gentime
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
  
  #plot_model(model)
  
  s_present <- schedule_sampling(model, times = 0, list(pop1, s1),list(pop2, s1),list(pop3, s1),list(pop4, s1),list(pop5, s1),list(pop6, s1),list(pop7, s1),list(pop8, s1),list(pop9, s1),list(pop10, s1))
  s_ancient1 <- schedule_sampling(model, times = n1anc, list(pop1, 1))
  s_ancient2 <- schedule_sampling(model, times = n2anc, list(pop2, 1))
  s_ancient3 <- schedule_sampling(model, times = n3anc, list(pop3, 1))
  s_ancient4 <- schedule_sampling(model, times = n4anc, list(pop4, 1))
  samples <- rbind(s_present,s_ancient1,s_ancient2,s_ancient3,s_ancient4)
  
  ts <- msprime(model, samples=samples, sequence_length = slength, recombination_rate = recombination_rate) %>%
    ts_mutate(mutation_rate = mutation_rate)
  
  ts_eigenstrat(ts, prefix = prefix)
  
  #f2 for the ancient samples
  t1 <- schedule_sampling(model, times = n1anc, list(pop1, 10000))
  t2 <- schedule_sampling(model, times = n2anc, list(pop2, 10000))
  t3 <- schedule_sampling(model, times = n3anc, list(pop3, 10000))
  t4 <- schedule_sampling(model, times = n4anc, list(pop4, 10000))
  
  true_samples <- rbind(t1,t2,t3,t4)
  ts_f2 <- msprime(model, samples=true_samples, sequence_length = slength, recombination_rate = recombination_rate)
  
  s_ts <- ts_samples(ts_f2) %>%
    split(., .$pop) %>%
    lapply(pull, "name")
  
  f12 = ts_f2(ts_f2, A=s_ts$pop1, B=s_ts$pop2,mode = "branch")$f2
  f13 = ts_f2(ts_f2, A=s_ts$pop1, B=s_ts$pop3,mode = "branch")$f2
  f14 = ts_f2(ts_f2, A=s_ts$pop1, B=s_ts$pop4,mode = "branch")$f2
  f23 = ts_f2(ts_f2, A=s_ts$pop2, B=s_ts$pop3,mode = "branch")$f2
  f24 = ts_f2(ts_f2, A=s_ts$pop2, B=s_ts$pop4,mode = "branch")$f2
  f34 = ts_f2(ts_f2, A=s_ts$pop3, B=s_ts$pop4,mode = "branch")$f2
  
  f3_213 = ts_f3(ts_f2, A=s_ts$pop2, B=s_ts$pop1, C=s_ts$pop3, mode = "branch")$f3
  f4_1234 = ts_f4(ts_f2,W = s_ts$pop1,X = s_ts$pop2,Y = s_ts$pop3, Z = s_ts$pop4)$f4
  
  f2mat = matrix(c(0,f12,f13,f14,f12,0,f23,f24,f13,f23,0,f34,f14,f24,f34,0),
                 nrow = 4, ncol = 4, byrow = FALSE)
  #ts_samples(ts)
  
  write.table(c(f12, f3_213, f4_1234), file = outfile, row.names=FALSE, col.names=FALSE, sep=',')
  write.table(round(f2mat,4), file = outf2, row.names=FALSE, col.names=FALSE, sep=',')
  
}
