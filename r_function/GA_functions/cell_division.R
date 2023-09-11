## ---------------------------
#   cell_division()
## ---------------------------

# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   scores = vector-like list of scores, one for each genome of this generation
#   

# OUTPUT: 
#   genomes of next generation
#

cell_division <- function(genomes, scores, n_best, n_child, n_elite, n_novel,
                          mu, cr) {
  if (length(genomes) != length(scores)) {
    stop('cell_division error: scores length should be equal to the number of genomes.')
  }
  # if ((n_best * n_child + n_novel) != length(genomes)) {
  #   stop('n_best * n_child - n_elite should be equal to the number of genomes')
  # }
  # Models are ordered by desc. scores values and number of features used
  sorted_models <- order(scores, sapply(genomes, sum))
  
  novel <- generate_G0(n_snps = length(genomes[[1]]),
                       n_indiv = n_novel,
                       p = mu)
  
  to_keep <- rep(sorted_models[1:n_best], times = n_child)

  is_elite <- unique(to_keep)[0:n_elite] # Works as long as n_elite is an integer.
  # is_elite <- to_keep[0:n_elite] # Identical elites do not occur in early gens.

  out_gen <- vector(mode = 'list', length = length(genomes))
  out_gen[1:(n_best * n_child)] <- genomes[to_keep]
  out_gen[((n_best * n_child)+1):((n_best * n_child)+n_novel)] <- novel

  out_gen[1:(length(out_gen) - n_elite)] <- 
    evolve(genomes = out_gen[1:(length(out_gen) - n_elite)], 
           mut_rate = mu, conjug_rate = cr)

    out_gen[((n_best * n_child)+n_novel+1):length(out_gen)] <- 
      genomes[is_elite]
  return(out_gen)
}

loss_to_BW <- function(scores, penalty_power) {
  # This function takes as input the vector of losses from the fitness function and transforms it 
  # into a vector whose sum is 1, each element is the probability for a genome to be selected by 
  # the Biased Wheel at each sampling.
  # scores is expected to be of same length as genomes with values between 0 and 1, minimal scores being best
  
  # A bigger penalty power favours better solutions more heavily
  penalized_scores <- (1+median(scores)-scores) ** penalty_power
  BW_probabilities <- penalized_scores/sum(penalized_scores)
  return(BW_probabilities)
}