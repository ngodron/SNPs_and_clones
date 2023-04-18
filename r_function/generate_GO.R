## ---------------------------
# generate_GO()
## ---------------------------

#
#
#

generate_G0 <- function(n_snps, n_indiv, p) {
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
  which(rowSums(out_gen0) == 0)
  for (i in 1:nrow(out_gen0)) {
    out_gen0[i, sample(x = 1:ncol(out_gen0), size = 1)] <- 1  
  }
  return(out_gen0)
}
