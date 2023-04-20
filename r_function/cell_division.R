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
  
  ranking <- 
    rank(x = scores, ties.method = 'random') 
  to_keep <- rep(which(ranking <= n_best), each = n_child)
  out_gen <- genomes[to_keep, ]
  out_gen <- evolve(genomes = out_gen, mut_rate = mu, conjug_rate = cr)
  out_gen <- rbind(out_gen, genomes[scores[which(ranking <= n_elite)]])
  return(out_gen)
}

