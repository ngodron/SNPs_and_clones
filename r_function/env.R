arguments$dir <- '~/Kevin_These/Current_File/'
# arguments$dir <- '../../Kevin_These/Current_File/'

arguments$snp <- '230516_matrix_SNP_gene_roary_594.csv'
arguments$pheno <- 'database_594.csv'
arguments$pheno_index <- 2
arguments$covar_index <- 7

arguments$verbose <- 2

params_list <- list()
# General parameters
## Number of generations
params_list$n_iter <- 1e1
## Number of individuals
params_list$n_ind <- 1e2 

# Population parameters
params_list$n_eli <- ceiling(params_list$n_ind / 20) # Count of elites
params_list$n_nov <- ceiling(params_list$n_ind / 20) # Count of novel (random) individuals
params_list$n_chi <- 4 # Number of children per "top genome"

params_list$n_top <- 
  (params_list$n_ind / params_list$n_chi) - 
  (params_list$n_eli / params_list$n_chi) -
  (params_list$n_nov / params_list$n_chi)
params_list$n_top <- ceiling(params_list$n_top)

params_list$n_ind <- (params_list$n_top * params_list$n_chi) +
  (params_list$n_eli + params_list$n_nov) # Total count of individuals

params_list$mutation_rate <- 1e-3
params_list$crossing_rate <- 0.5

fitness_fun <- decision_tree_fitness
