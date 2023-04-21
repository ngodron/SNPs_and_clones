## ---------------------------
#   evolve()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

mutation <- function(genome, mu) {
  genome <- as.logical(genome)
  to_mutate <- 
    as.logical(rbinom(n = length(genome), size = 1, prob = mu))
  genome <- xor(to_mutate, genome)
  return(genome)
}

crossing_over <- function(genomes, cr) {
  # Randomised order pairwise single crossing over.
  # To do: Location of crossing over is picked in a Gaussian centered around |SNPs|/2
  genomes <- genomes[sample(1:nrow(genomes)), ]
  crossed <- genomes
  for (i in 1:(nrow(genomes))/2) {
    # print("Before")
    # print(genomes[(2*i-1):(2*i),])
    if (rbinom(1, 1, cr)) {
     halfway <- floor(ncol(genomes)/2)
     crossed[(2*i-1),] <- c(genomes[(2*i-1), 1:halfway],
                                 genomes[(2*i), (halfway+1):ncol(genomes)])
     crossed[(2*i),] <- c(genomes[(2*i), 1:halfway],
                               genomes[(2*i-1), (halfway+1):ncol(genomes)])
     
   } else {
     crossed[(2*i-1):(2*i), ] <- genomes[(2*i-1):(2*i), ]
   }
    # print("After")
    # print(crossed[(2*i-1):(2*i),])
  }
  return(crossed)
}

evolve <- function(genomes, mut_rate, conjug_rate) {
  for (i in 1:nrow(genomes)) {
    genomes[i, ] <- mutation(genome = genomes[i, ], mu = mut_rate)
  }
  genomes <- crossing_over(genomes <- genomes, cr = conjug_rate)
  return(genomes)
}