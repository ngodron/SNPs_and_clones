## Sourcing functions ----
sapply(X = list.files(path = './r_function/GA_functions', 
                      full.names = TRUE), 
       FUN = source)

## Arguments parsing ----
config_file <- commandArgs(trailingOnly = TRUE)[1]
n_iter <- commandArgs(trailingOnly = TRUE)[2]
count_iter <- commandArgs(trailingOnly = TRUE)[3]


if (file.exists(config_file)) {
  source(config_file)
} else {
  stop("The path to config file is either wrong or unreadable")
}

# Unpacking arguments and parameters into global environment
unpacking <- function(list_args) {
  for (i in 1:length(list_args)) {
    assign(names(list_args)[i], list_args[[i]], envir = .GlobalEnv)
    cat(names(list_args)[i], ":", list_args[[i]], "\n")
  }
} 

unpacking(arguments)
unpacking(params_list)

## Inputs ----
source('./r_function/load_inputs.R')
matrices <- input_loading(arguments$snp,
                          arguments$pheno, 
                          arguments$pheno_index,
                          arguments$covar_index)

snp_df <- matrices[[1]]
pheno <- matrices[[2]]
if (length(matrices) == 3) {
  covar <- matrices[[3]]
}

# rm(matrices)

if (file.exists('./curr_gen.GAG')) {
  source('curr_gen.GAG')
} else {
  curr_gen <-
    generate_G0(n_snps = ncol(snp_df), n_indiv = n_ind, p = mutation_rate)
}

## Variables declaration ----
print_params <- 
  paste0('n_ind = ', n_ind, '\nmu = ', mutation_rate)

all_gen <- vector(mode = 'list', length = n_iter)
score_list <- vector(mode = 'list', length = n_iter)
model_list <- vector(mode = 'list', length = n_iter)


# Temporary /!\ ----
# debugSource("~/2023/All_list/SNPs_and_clones/r_function/list_version/calc_score_all_list.R")
pheno <- as.integer(pheno != "sputum")
prop_priors <- sapply(X = sort(unique(pheno)), function(x) sum(pheno == x)) / length(pheno)
save <- 1

## G.A. iterator: gen_algo ----
gen_algo <- function(snp_matrix, pheno_matrix, covar_matrix, parameters, fitness_fun) {
  unpacking(parameters)
  cat("\n", "\n")
  
  for (i in 1:n_iter) {
    cat('Generation ', i, '/', n_iter, '\n')
    all_gen[[i]] <- 
      sapply(curr_gen, function(x) which(x == 1, arr.ind = TRUE))
    all_gen[[i]] <- sapply(all_gen[[i]], function(x) colnames(snp_df)[x])
    
    print(summary(sapply(all_gen[[i]], length)))
    
    curr_scores_models <- 
      calc_score(genomes = curr_gen, 
                 snps = snp_df, 
                 phenotype = pheno, 
                 fitness = fitness_fun, 
                 covars = covar,
                 weights = NULL)
    
    curr_scores <- 
      unlist(lapply(curr_scores_models, function(x) x[[1]]))
    curr_models <- 
      (lapply(curr_scores_models, function(x) x[[2]]))
    model_list[[i]] <-
      curr_models[[which(curr_scores == min(curr_scores))[1]]] 
    rm(curr_scores_models)
    print(summary(curr_scores))
    score_list[[i]]<- c(curr_scores)
    #diversity[i] <- genomes_diversity(curr_gen)
    
    next_gen <- 
      cell_division(genomes = curr_gen, 
                    scores = curr_scores, 
                    n_best = n_top, 
                    n_child = n_chi, 
                    n_elite = n_eli, 
                    n_novel = n_nov,
                    mu = mutation_rate, cr = crossing_rate)
    # print(length(next_gen))
    curr_gen <- next_gen
  }
output <- list(all_gen, score_list, model_list, curr_gen)
  return(output)
}

output_algo <- gen_algo(snp_df, pheno, covar, params_list, decision_tree_fitness)
# 
all_gen <- output_algo[[1]]
score_list <- output_algo[[2]]
model_list <- output_algo[[3]]
curr_gen <- output_algo[[4]]

if (save == 1) {
  dump("curr_gen", file = "curr_gen.GAG")
}
