## ---------------------------
#   evolve()
## ---------------------------

# INPUT:
#   genomes of gen i: 0/1 matrix with individuals in rows and SNPs in cols
#   mut_rate: between 0 and 1, probability of swapping each feature's (SNP) bit.
#   conjug_rate: between 0 and 1, probability for each pair of genomes to be crossed over.
#   min_mut: minimum count of mutations to be realised. (To implement)

# OUTPUT:
#   genomes of gen i+1, after mutation and crossing-over.

mutation_balanced <- function(genome, mu, max_tot) {
  # Mutation function to have balanced 0->1 and 1->0 mutations between generations,
  # but with diversity in mutation rates and gain/loss balance at the scale of a genome.
  # This function avoids the overall inflation of features used by max_tot threshold.
  # NB: max_tot can be surpassed after crossing-over is realised.
  
  genome_length <- length(genome)
  n_set_mut <- sum(genome) # Count of features used
  n_unset_mut <- genome_length - n_set_mut # Count of features NOT used
  
  # Binomial sampling of count of mutation events for each genome:
  n_mut <- 
    rbinom(n = 1, size = genome_length, prob = mu)
  
  # Gain/Loss mutations counts allocation (at random from n_mut):
  mut_gain <- rbinom(n = 1, size = n_mut, prob = 0.5)
  mut_loss <- n_mut - mut_gain
  mut_gain <- min(mut_gain, n_unset_mut)
  mut_loss <- min(mut_loss, n_set_mut)
  
  total_mut <- n_set_mut + mut_gain - mut_loss 
  # cat('n_mut: ', n_mut)
  # cat('\t\ttotal_mut: ', total_mut,'\n')
  # cat('mut gain \t', mut_gain, '\t\tmut_loss\t', mut_loss, '\n')
  
  if (total_mut > max_tot) {
    mut_loss <- mut_loss + mut_gain
    mut_loss <- min(mut_loss, n_set_mut)
    mut_gain <- 0
  }
  
  # Making Gain and Loss mutations:
  to_mutate <- c(sample(x = which(genome == 0), size = mut_gain, replace = FALSE),
    sample(x = which(genome == 1), size = mut_loss, replace = FALSE))
  # print(to_mutate)
  genome[to_mutate] <- ! genome[to_mutate]
  
  return(genome)
}

crossing_over <- function(genomes, cr) {
  # Randomised order pairwise single crossing over.
  genomes <- genomes[sample(1:length(genomes))] # randomizes genomes order
  crossed <- genomes
  for (i in 1:(length(genomes)/2)) {
    if (rbinom(1, 1, cr)) {
      pivot <- floor(runif(1, min = 2, max = length(genomes[[1]])-1))
      gen_1 <- 2*i - 1
      gen_2 <- 2*i
      crossed_1 <- 
         c(genomes[[gen_1]][1:pivot], genomes[[gen_2]][(pivot+1):length(genomes[[1]])])
      crossed_2 <- 
        c(genomes[[gen_2]][1:pivot], genomes[[gen_1]][(pivot+1):length(genomes[[1]])])
      
    if (rbinom(n = 1, size = 1, prob = cr / 100)) { # rare inversion event
      crossed[c(gen_1, gen_2)] <- list(crossed_2, crossed_1)
    } else {
      crossed[c(gen_1, gen_2)] <- list(crossed_1, crossed_2)
      }
    }
  }
  return(crossed)
}

evolve <- function(genomes, mut_rate, conjug_rate, min_mut) {
  out_genomes <- vector(mode = 'list', length = length(genomes))
  
  for (i in 1:length(genomes)) {
    out_genomes[[i]] <- mutation_balanced(genome = genomes[[i]],
                                      mu = mut_rate,
                                      max_tot = 20)
  }
  
  if (! length(genomes) > 2) {
    warning("Number of SNPs in genomes is not superior to 2, crossing-over was not performed.")
  }
  else {
    out_genomes <- 
      crossing_over(genomes = out_genomes, cr = conjug_rate)
  }
  return(out_genomes)
}
