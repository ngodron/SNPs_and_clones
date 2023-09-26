

## GA iterator: gen_algo ----
gen_algo <- function(snp_matrix, 
                     pheno_matrix, 
                     covar_matrix = NULL, 
                     costs = NULL,
                     n_gen = 1e2,
                     n_iter = 1,
                     n_ind = 1e2,
                     n_eli = 2,
                     n_nov = 5,
                     n_chi = 4,
                     mu = 1e-3,
                     cr = 0.5,
                     curr_gen = NULL,
                     fitness_fun = decision_tree_fitness) {
  ## Sourcing functions ----
  sapply(X = list.files(path = './r_function/GA_functions', 
         full.names = TRUE), 
         FUN = source)
  cat("\n", "\n")

  
  # Variables
  all_gen <- vector(mode = 'list', length = n_iter)
  score_list <- vector(mode = 'list', length = n_iter)
  model_list <- vector(mode = 'list', length = n_iter)
  
  if (is.null(covar_matrix)) {
    covar_matrix <- matrix(data = 1, nrow = nrow(snp_matrix))
    colnames(covar_matrix) <- 'COVAR_1'
  }
  if (is.null(costs)) {
    costs <- rep(1, (ncol(snp_matrix)+ ncol(covar_matrix)))
    names(costs) <- c(colnames(snp_matrix), colnames(covar_matrix))
    costs <- as.matrix(costs)
  }
  
  # Population parameters ----
  n_top <- 
    (n_ind / n_chi) - 
    (n_eli / n_chi) -
    (n_nov / n_chi)
  n_top <- ceiling(n_top)
  
  n_ind <- (n_top * n_chi) +
    (n_eli + n_nov) # Total count of individuals
  
  # First generation ----
  if (is.null(curr_gen)) {
    curr_gen <- generate_G0(n_snps = ncol(snp_matrix), 
                n_indiv = n_ind, 
                p = mu)
  } else if (exists(x = curr_gen)) {
    #TODO
    cat('Loading ...')
    
  }
  
  # All other generations ----
  count_iter <- 0
  for (i in 1:n_gen) {
    cat('Generation ', count_iter + i, '/', n_gen * n_iter, '\n')
    all_gen[[i]] <- 
      sapply(curr_gen, function(x) which(x == 1, arr.ind = TRUE))
    all_gen[[i]] <- sapply(all_gen[[i]], function(x) colnames(snp_matrix)[x])
    
    print(summary(sapply(all_gen[[i]], length)))
    
    curr_scores_models <- 
      calc_score(genomes = curr_gen, 
                 snps = snp_matrix, 
                 phenotype = pheno_matrix, 
                 covars = covar_matrix,
                 fitness = fitness_fun,
                 costs = costs,
                 met = 'loss')
    
    curr_scores <- 
      unlist(lapply(curr_scores_models, function(x) x[[1]]))
    curr_models <- 
      (lapply(curr_scores_models, function(x) x[[2]]))
    model_list[[i]] <-
      curr_models[[which(curr_scores == min(curr_scores))[1]]] 
    rm(curr_scores_models)
    print(summary(curr_scores))
    score_list[[i]]<- c(curr_scores)
    
    next_gen <- 
      cell_division(genomes = curr_gen, 
                    scores = curr_scores, 
                    n_best = n_top, 
                    n_child = n_chi, 
                    n_elite = n_eli, 
                    n_novel = n_nov,
                    mu = mu, cr = cr)
    # print(length(next_gen))
    curr_gen <- next_gen
  }
  
  score_df <-
    data.frame(gen = rep(1:length(score_list), each = n_ind),
               score = unlist(score_list),
               n_snps = unlist(sapply(all_gen, FUN = function(x) {
                 sapply(x, length)}
                 , simplify = FALSE))
    )
  output <- list(all_gen, score_list, model_list, curr_gen, score_df)
  return(output)
}