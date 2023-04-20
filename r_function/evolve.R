## ---------------------------
#   calc_score()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)


evolve <- function(genomes, mut_rate, conjug_rate) {
  for (i in 1:nrow(genomes)) {
    genomes[i, ] <- mutation(genome = genomes[i, ], mu = mut_rate)
  }
  return(genomes)
}

mutation <- function(genome, mu) {
  genome <- as.logical(genome)
  to_mutate <- 
    as.logical(rbinom(n = length(genome), size = 1, prob = mu))
  genome <- xor(to_mutate, genome)
  return(genome)
}

crossing_over <- function(genome, mu) {}