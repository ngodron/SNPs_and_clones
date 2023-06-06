## ---------------------------
#   calc_score()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

calc_score <- function(genomes, snps, phenotype, fitness, covars, weights = NULL) {
  if (is.null(weights)) {
    weights <- rep(1, ncol(snps))
    names(weights) <- colnames(snp_df)
  }
  all_scores <- vector(mode = 'list', length = length(genomes))
  for (i in 1:length(genomes)) {
    # drop = F allows to keep single SNP genomes as matrix
    curr_snps <- snps[ , genomes[[i]] == 1, drop = FALSE] 
    all_scores[[i]] <-  
      fitness(snps = curr_snps, 
              pheno = phenotype, 
              genome_size = ncol(genomes), 
              covariables = covars,
              wgts = weights) 
  }
  return(all_scores)
}


decision_tree_fitness <- function(snps, pheno, genome_size, covariables, wgts) {
  if (ncol(snps) == 0) return(list(1, NA)) # NA for the model
  df <- data.frame(pheno, snps, covariables)
  #df <- as.data.frame(apply(X = df, MARGIN = 2, FUN = as.logical, simplify = TRUE))
  wgts <- wgts[colnames(snps)] # filtering for wanted snps
  wgts <- c(wgts,rep(mean(wgts), ncol(covariables))) # adding mean weight for the covariables
  wgts <- 1/wgts
  formu <- formula(pheno ~ .)
  model <- 
    rpart::rpart(formula = formu, 
                 data = df, 
                 minbucket = 10, 
                 method = 'class',
                 cost = wgts,
                 parms = list(prior = prop_priors),
                 model = TRUE)
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

genomes_diversity <- function(genomes) {
  n_by_cols <- colSums(genomes)
  genomes <- genomes[ ,n_by_cols >= 1, drop = FALSE]
  prop <- colSums(genomes) / nrow(genomes)
  out <- 1 - abs(mean(prop) - 0.5) * 2
  return(out)
}