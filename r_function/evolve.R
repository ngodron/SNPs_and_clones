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

mutation_balanced <- function(genome, mu, mut_max = 5, min_mut = 1) {
  # Mutation function to have balanced 0->1 and 1->0 mutations between generations,
  # but with diversity in mutation rates and gain/loss balance at the scale of a genome.
  # This function is implemented to avoid the overall inflation of SNPs used.
  
  # print(genome) # Debugging
  # Binomial sampling of count of mutation events for each genome:
  # min_mut is the minimum count of mutations acquired by a genome.
  mut_count <- max(sum(rbinom(mut_max, 1, 0.5)), min_mut)
    
  # Gain/Loss mutations counts allocation (from mut_count):
  gain_count <- sum(rbinom(mut_count, 1, 0.5))
  loss_count <- mut_count - gain_count
  # print(c(gain_count, loss_count)) # Debugging
  # print(c(sum(genome == 0)-1, sum(genome == 1)-1)) # Debugging
    
  # Making Gain and Loss mutations:
  zero_bits_to_swap <- sample(which(genome == 0), min(gain_count, sum(genome == 0)-1))
  one_bits_to_swap <- sample(which(genome == 1), min(loss_count, sum(genome == 1)-1))
  # print(c(zero_bits_to_swap, one_bits_to_swap)) # Debugging
  genome[c(zero_bits_to_swap, one_bits_to_swap)] <- 
    ! genome[c(zero_bits_to_swap, one_bits_to_swap)]
  
  # hist(mut_count, breaks = unique(mut_count))
  # counts <- data.frame(mut_count, gain_count, loss_count)
  return(genome)
}

crossing_over <- function(genomes, cr) {
  # Randomised order pairwise single crossing over.
  # To do: Location of crossing over is picked in a Gaussian centered around |SNPs|/2
  genomes <- genomes[sample(1:nrow(genomes)), ]
  crossed <- genomes
  for (i in 1:(nrow(genomes)/2)) {
    # print("Before")
    # print(genomes[(2*i-1):(2*i),])
    if (rbinom(1, 1, cr)) {
     halfway <- floor(ncol(genomes)/2)
     gen_1 <- 2*i - 1
     gen_2 <- 2*i
     crossed_1 <- 
       c(genomes[gen_1, 1:halfway], genomes[gen_2, (halfway+1):ncol(genomes)])
     crossed_2 <- 
       c(genomes[gen_2, 1:halfway], genomes[gen_1, (halfway+1):ncol(genomes)])
     crossed[c(gen_1, gen_2), ] <- rbind(crossed_1, crossed_2)
    } 
   #  else {
   #   crossed[(2*i-1):(2*i), ] <- genomes[(2*i-1):(2*i), ]
   # }
    # print("After")
    # print(crossed[(2*i-1):(2*i),])
  }
  return(crossed)
}

evolve <- function(genomes, mut_rate, conjug_rate, min_mut) {
  maximum_mutations <- round((ncol(genomes) * mut_rate * 2))
  for (i in 1:nrow(genomes)) {
    genomes[i, ] <- mutation_balanced(genome = genomes[i, ],
                                      mu = mut_rate,
                                      mut_max = maximum_mutations,
                                      min_mut = 1)
  }
  
  # for (i in 1:nrow(genomes)) {
  #   genomes[i, ] <- mutation(genome = genomes[i, ], mu = mut_rate)
  # }
  
  genomes <- 
    crossing_over(genomes <- genomes, cr = conjug_rate)
  return(genomes)
}