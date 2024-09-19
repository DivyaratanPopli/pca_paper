
sim_func <- function(sp, Ne, gentime, st_array, mgen_array, madm_array, mgen, mprop, mutation_rate, recombination_rate, slength, s1, genfile, truef2){
  
  s1=as.integer(as.integer(s1)/10)
  Ne=as.integer(Ne)
  sp=as.integer(sp)
  mgen=mgen*sp*gentime
  
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
  
  samples <- schedule_sampling(model, times = 0, list(pop1, s1),list(pop2, s1),list(pop3, s1),list(pop4, s1),list(pop5, s1),list(pop6, s1),list(pop7, s1),list(pop8, s1),list(pop9, s1),list(pop10, s1))
  
  ts <- msprime(model, samples=samples, sequence_length = 1e8, recombination_rate = 1e-8) %>%
    ts_mutate(mutation_rate = 1e-8)
  
  ts_eigenstrat(ts, prefix = prefix)
  
  #get the tree
  ts_samples(ts)
  
  s_ts <- ts_samples(ts1) %>%
    split(., .$pop) %>%
    lapply(pull, "name")
  
  ts_f2(ts1, A=s_ts$pop6, B=s_ts$pop5,mode = "branch")
  
  tree <- ts_tree(ts, i = 1)
  ts_draw(tree, time_scale = "rank")
}



ts_1000 <- msprime(model, sequence_length = 1e8, recombination_rate = 1e-8) 

s_ts <- ts_samples(ts_1000) %>%
  split(., .$pop) %>%
  lapply(pull, "name")

ts_f2(ts_1000, A=s_ts$pop1, B=s_ts$pop2,mode = "branch")
