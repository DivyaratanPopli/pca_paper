library(admixtools)

convert_geno <- function(genof,outf){

  genof1=strsplit(genof, split=".geno")[[1]]
  xx=read_eigenstrat(genof1)
  xx$geno[is.na(xx$geno)] =9

  write.table(xx$geno, file = outf, row.names=FALSE, col.names=FALSE, sep=',')

}
convert_geno(genof=snakemake@input[["genof"]], outf=snakemake@output[["outf"]])

#genof="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times1000/mu0/run1/nind100_eigen"
#outf="/mnt/diversity/divyaratan_popli/fstats/genetic_simulations/slendr_anc/simfiles/Ne1000/split_times1000/mu0/run1/nind100_eigen.geno_pc"
