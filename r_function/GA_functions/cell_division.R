## ---------------------------
#   cell_division()
## ---------------------------

# INPUT:
#   
#   
#   

# OUTPUT: 
#   
#

cell_division <- function(genomes, scores, n_best, n_child, n_elite, n_novel,
                          mu, cr) {
  if (length(genomes) != length(scores)) {
    stop('scores length should be equal to the number of genomes')
  }
  # if ((n_best * n_child + n_novel) != length(genomes)) {
  #   stop('n_best * n_child - n_elite should be equal to the number of genomes')
  # }
  # Models are ordered by desc. scores values and number of features used
  sorted_models <- order(scores, sapply(genomes, sum))
  to_keep <- rep(sorted_models[1:n_best], times = n_child)
  
  novel <- generate_G0(n_snps = length(genomes[[1]]), 
                       n_indiv = n_novel, 
                       p = mu)
  
  is_elite <- unique(to_keep)[0:n_elite] # Works on integers.
  # is_elite <- to_keep[0:n_elite] # Identical elites do not occur in early gens.

  out_gen <- vector(mode = 'list', length = length(genomes))
  out_gen[1:(n_best * n_child)] <- genomes[to_keep]
  out_gen[((n_best * n_child)+1):((n_best * n_child)+n_novel)] <- 
    novel
  ####### LAST THING TO DO HERE IS TO CORRECT EVOLVE
  out_gen[1:(length(out_gen) - n_elite)] <- 
    evolve(genomes = out_gen[1:(length(out_gen) - n_elite)], 
           mut_rate = mu, conjug_rate = cr)
  #######
    out_gen[((n_best * n_child)+n_novel+1):length(out_gen)] <- 
      genomes[is_elite]
  return(out_gen)
}

