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

cell_division <- function(genomes, scores, n_best, n_child, n_elite, mu, cr) {
  if (nrow(genomes) != length(scores)) {
    stop('scores length should be equal to genomes number of rows')
  }
  if ((n_best * n_child + n_elite) != nrow(genomes)) {
    stop('n_best * n_child - n_elite should be equal to genomes number of rows')
  }
  # deserves a comment
  to_keep <- rep(order(scores, rowSums(genomes))[1:n_best], n_child)
  is_elite <- to_keep[0:n_elite]
  out_gen <- genomes[to_keep, ]
  out_gen <- evolve(genomes = out_gen, mut_rate = mu, conjug_rate = cr)
  out_gen <- rbind(out_gen, genomes[is_elite, ])
  return(out_gen)
}

