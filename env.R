# Meta arguments ----
arguments <- list()

# Directory in which are the SNP and phenotype files.
arguments$dir <- '~/Kevin_These/Current_File/'
# arguments$dir <- '../../Kevin_These/Current_File/'

# Path to SNP presence/absence tab-separated file (absolute, or relative to --dir).
arguments$snp <- '230615_SNP_roary_594.csv'
# Path to phenotype tab-separated file (absolute, or relative to --dir).
arguments$pheno <- 'database_594.csv'

# Index of phenotype and optional covariate(s) columns in pheno file.
arguments$pheno_index <- 2
arguments$covar_index <- 7

# Path to (optional) cost tab-separated file (absolute, or relative to --dir).
arguments$cost <- '2023_06_29-homoplasy.tsv'

# Index of optional weight column in weights file.
arguments$cost_index <- 4

# Optional verbose argument can take values 0, 1 or 2.
# To do: Make it functional, set default to 1 in gen_algo.R
# 0: Neither input nor genetic algorithm logs are printed to console.
# 1: Only genetic algorithm logs are printed to console.
# 2: Input and genetic algorithm are printed to console.
arguments$verbose <- 2

# Genetic algoritm parameters ----
params_list <- list()

## Number of individuals
params_list$n_ind <- 1e2

# Population parameters: number of elites, novel individuals and children per individual.
params_list$n_eli <- ceiling(params_list$n_ind / 20) # Count of elites
params_list$n_nov <- ceiling(params_list$n_ind / 20) # Count of novel (random) individuals
params_list$n_chi <- 4 # Number of children per "top genome"

# Determining number of "top" individuals that will reproduce
params_list$n_top <- 
  (params_list$n_ind / params_list$n_chi) - 
  (params_list$n_eli / params_list$n_chi) -
  (params_list$n_nov / params_list$n_chi)
params_list$n_top <- ceiling(params_list$n_top)

# Number of individuals can be increased
# to satisfy all conditions with an equal number of children 
params_list$n_ind <- (params_list$n_top * params_list$n_chi) +
  (params_list$n_eli + params_list$n_nov) # Total count of individuals

# Mutation and crossing-over influence evolution between Gen(n) and Gen(n+1)
params_list$mutation_rate <- 1e-3
params_list$crossing_rate <- 0.5

# Fitness function choice.
# Advanced users can add their custom fitness functions to calc_score.R
fitness_fun <- decision_tree_fitness
