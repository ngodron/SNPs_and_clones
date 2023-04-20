## ---------------------------
#   generate_GO()
## ---------------------------

# INPUT:
#   n_snps,
#   n_indiv,
#   p, vector of probability

# OUTPUT: 
#   out_gen0, 0/1 matrix with indiv
#

generate_G0 <- function(n_snps, n_indiv, p) {
  # make p a matrix
  if (length(p) != n_snps) {
    if (length(p) > 1) {
      stop("Wrong number of elements in vector p:
    Should be either a single value or have a length equal to n_snps")
    }
    if (any(p < 0 | p > 1)) {
      stop("Wrong value(s) in vector p:
      Should be included in [0-1].")
    }
  } 
  out_gen0 <- 
    t(replicate(n = n_indiv, expr = rbinom(n = n_snps, size = 1, prob = p)))
  
  # each genome should include 1 snp at least
  for (i in 1:nrow(out_gen0)) {
    if(sum(out_gen0[i,] == 0)) {
      out_gen0[i, sample(x = 1:ncol(out_gen0), size = 1)] <- 1
    }
  }
  return(out_gen0)
}
