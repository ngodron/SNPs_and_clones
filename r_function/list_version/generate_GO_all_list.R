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
    lapply(X = 1:n_indiv, 
           FUN = function(x) rbinom(n = n_snps, size = 1, prob = p))
  
  # each genome should include 1 snp at least
  for (i in 1:n_indiv) {
    if(sum(out_gen0[[i]] == 0)) {
      out_gen0[[i]][sample(x = 1:n_snps, size = 1)] <- 1
    }
  }
  return(out_gen0)
}
