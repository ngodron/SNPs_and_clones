## ---------------------------
#   calc_score()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

calc_score <- function(genomes, snps, phenotype, fitness) {
  all_scores <- rep(NA_real_, nrow(genomes))
  for (i in 1:nrow(genomes)) {
    curr_genome <- which(genomes[i, ] == 1)
    # drop = F allows to keep single SNP genomes as matrix
    curr_snps <- snps[ , curr_genome, drop = FALSE] 
    all_scores[i] <-  fitness(snps = curr_snps, pheno = phenotype) 
  }
  return(all_scores)
}

fish_fitness_mock <- function(snps, pheno) {
  out_score <- rep(NA_real_, ncol(snps))
  for (i in 1:ncol(snps)) {
    curr_snp <- as.matrix(snps[ , i])
    curr_p <- fisher.test(table(curr_snp, pheno))$p.value
    out_score[i] <- curr_p
  }
  out_score <- mean(out_score) #+ log10(ncol(snps))
  return(out_score)
}

glm_fitness_mock <- function(snps, pheno) {
   
  if (ncol(snps) == 0) return(1) 
  
  df <- data.frame(pheno, snps)
  
  model <- glm(formula = pheno ~ ., data = df, family = binomial(link = "logit"))
  predicted <- predict(object = model, type = 'response')
  predictions <- ifelse(test = predicted >= 0.5, yes = TRUE, no = FALSE)

  pred_TP <- sum(pheno & predictions)
  pred_TN <- sum(! pheno & ! predictions)
  pred_FP <- sum(! pheno & predictions)
  pred_FN <- sum(pheno & ! predictions)
  
  accu <- (pred_TP + pred_TN) / (pred_TP + pred_TN + pred_FP + pred_FN)
  
  out_score <- 1 - accu
  if (ncol(snps) > 20) out_score <- 1
  out_score <- out_score + 2^(ncol(snps)/1000)
   
  # out_score <- ifelse(test = ncol(snps) > 15, yes =  out_score * ncol(snps), no = out_score)
  return(out_score)
}

rand_fitness_mock <- function(snps, pheno) {
  return(runif(n = 1))
}
