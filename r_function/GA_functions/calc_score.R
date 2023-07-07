## ---------------------------
#   calc_score()
## ---------------------------

# 
# INPUT:
#   genomes = 0/1 matrix with indiv in rows and snps in cols
#   snps = 0/1 matrix with indiv in rows and snps in cols
#   pheno = 0/1 vector of length = nrow(snps)

calc_score <- function(genomes, snps, phenotype, covars,
                       fitness, costs = NULL, met = "mcc") {
  if (is.null(costs)) {
    costs <- rep(1, (ncol(snps)+ ncol(covars)))
    names(costs) <- c(colnames(snp_df), colnames(covars))
    costs <- as.matrix(costs)
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
              costs_split = costs,
              metric = "mcc") 
  }
  return(all_scores)
}


decision_tree_fitness <- function(snps, pheno, genome_size, covariables, 
                                  costs_split, metric) {
  if (ncol(snps) == 0) return(list(1, NA)) # NA for the model
  df <- data.frame(pheno, snps, covariables)
  to_keep <- c(colnames(snps), colnames(covariables))
  #df <- as.data.frame(apply(X = df, MARGIN = 2, FUN = as.logical, simplify = TRUE))
  costs_split <- costs_split[to_keep,] # filtering for wanted snps
  # costs_split <- c(costs_split, rep(1, ncol(covariables))) # adding weight 1 for the covariables
  formu <- formula(pheno ~ .)
  model <- 
    rpart::rpart(formula = formu, 
                 data = df, 
                 minbucket = 10, 
                 method = 'class',
                 cost = costs_split,
                 model = TRUE)
  
  predicted <- predict(object = model)
  predictions <- predicted[ , 1] <= 0.5 # predicted[ ,1] is smallest alphanum value
  # predictions <- predicted
  TP <- sum(pheno & predictions)
  TN <- sum(! pheno & ! predictions)
  FP <- sum(! pheno & predictions)
  FN <- sum(pheno & ! predictions)
  
  if (metric == "loss") {
    accu <- (TP + TN) / (TP + TN + FP + FN)
    out_score <- 1 - accu
  }
  if (metric == "mcc") {
    margins <- ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    if (margins  == 0 | is.na(margins)){
      mcc <- 0
    } else {
      mcc <- (TP*TN - FP*FN) / ((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))**0.5
    }
    out_score <- 1 - mcc 
  }

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
