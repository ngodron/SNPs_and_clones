## ---------------------------
# generate_GO()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

calc_score <- function(genomes, snps, phenotype, fitness) {
  all_scores <- rep(NA, nrow(genomes))
  for (i in 1:nrow(genomes)) {
    curr_genome <- genomes[i, ] == 1
    curr_snps <- snps[ , curr_genome]
    all_scores[i] <-  fitness(snps = curr_snps, pheno = phenotype) 
  }
  return(all_scores)
}

fish_fitness_mock <- function(snps, pheno) {
  out_score <- 0
  for (i in 1:ncol(snps)) {
    
    curr_snp <- as.matrix(snps[ , i])
    curr_p <- fisher.test(table(curr_snp, pheno))$p
    out_score <- 
      out_score + (1- curr_p)
  }
  out_score <- out_score / ncol(snps)
  return(out_score)
}

