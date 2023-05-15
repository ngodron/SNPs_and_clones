## ---------------------------
#   calc_score()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

calc_score_nopar <- function(genomes, snps, phenotype, fitness, covars) {
  all_scores <- vector(mode = 'list', length = nrow(genomes))
  for (i in 1:nrow(genomes)) {
    curr_genome <- which(genomes[i, ] == 1)
    # drop = F allows to keep single SNP genomes as matrix
    curr_snps <- snps[ , curr_genome, drop = FALSE] 
    all_scores[[i]] <-  
      fitness(snps = curr_snps, 
              pheno = phenotype, 
              genome_size = length(genomes), 
              covariables = covars) 
  }
  return(all_scores)
}

calc_score <- function(genomes, snps, phenotype, fitness, covars) {

  # all_scores <- rep(NA_real_, nrow(genomes))
  all_scores <- 
    foreach(i = 1:nrow(genomes), .combine = 'c') %dopar% {
      curr_genome <- which(genomes[i, ] == 1)
      # drop = F allows to keep single SNP genomes as matrix
      curr_snps <- snps[ , curr_genome, drop = FALSE] 
      fitness(snps = curr_snps, pheno = phenotype, genome_size = ncol(genomes), covariables = covars) 
  }
  return(all_scores)
}

glm_fitness_mock <- function(snps, pheno, genome_size, covariables) {
  
  if (ncol(snps) == 0) return(list(1, NA)) # NA for the model
  
  df <- data.frame(pheno, snps, covariables)
  
  model <- glm(formula = pheno ~ ., data = df, family = binomial(link = "logit"))
  predicted <- predict(object = model, type = 'response')
  predictions <- ifelse(test = predicted >= 0.5, yes = TRUE, no = FALSE)
  
  pred_TP <- sum(pheno & predictions)
  pred_TN <- sum(! pheno & ! predictions)
  pred_FP <- sum(! pheno & predictions)
  pred_FN <- sum(pheno & ! predictions)
  
  accu <- (pred_TP + pred_TN) / (pred_TP + pred_TN + pred_FP + pred_FN)
  
  # num_mcc <- pred_TN * pred_TP - pred_FN * pred_FP
  # den_mcc <- 
  #   sqrt((pred_TP + pred_FP) * 
  #          (pred_TP + pred_FN) * 
  #          (pred_TN + pred_FP) * 
  #          (pred_TN + pred_FN))
  # 
  # mcc <- num_mcc / den_mcc
  out_score <- 1 - accu
  # out_score <- 1 - mcc
  # if (ncol(snps) > 30) out_score <- 1
  # if(ncol(snps) > 20) out_score <- 1
  out_score <- out_score  #* (ncol(snps) / genome_size/2)
  
  # out_score <- ifelse(test = ncol(snps) > 15, yes =  out_score * ncol(snps), no = out_score)
  return(list(out_score, model))
}

decision_tree_fitness <- function(snps, pheno, genome_size, covariables) {
  
  if (ncol(snps) == 0) return(list(1, NA)) # NA for the model
  
  df <- data.frame(pheno, snps, covariables)
  #df <- as.data.frame(apply(X = df, MARGIN = 2, FUN = as.logical, simplify = TRUE))
  formu <- formula(pheno ~ .)
  model <- 
    rpart::rpart(formula = formu, 
                 data = df, 
                 minbucket = 10, 
                 method = 'class')
  predicted <- predict(object = model)
  predictions <- predicted[ , 1] <= 0.5 # predicted[ ,1] is smallest alphanum value
  # predictions <- predicted
  pred_TP <- sum(pheno & predictions)
  pred_TN <- sum(! pheno & ! predictions)
  pred_FP <- sum(! pheno & predictions)
  pred_FN <- sum(pheno & ! predictions)
  
  accu <- (pred_TP + pred_TN) / (pred_TP + pred_TN + pred_FP + pred_FN)
  
  out_score <- 1 - accu
  rm(formu)
  return(list(out_score, model))
}

rand_fitness_mock <- function(snps, pheno) {
  return(runif(n = 1))
}

genomes_diversity <- function(genomes) {
  n_by_cols <- colSums(genomes)
  genomes <- genomes[ ,n_by_cols >= 1, drop = FALSE]
  prop <- colSums(genomes) / nrow(genomes)
  out <- 1 - abs(mean(prop) - 0.5) * 2
  return(out)
}
