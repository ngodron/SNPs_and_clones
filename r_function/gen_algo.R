library(profvis)
library(pryr)
library(tidyverse)

arguments <- list()

if (file.exists('./r_function/env.R')) {
  source('./r_function/env.R')
} else {
  # Directory in which are the SNP and phenotype files.
  arguments$dir <- commandArgs(trailingOnly = TRUE)[1]
  
  # Path to SNP presence/absence tab-separated file (absolute, or relative to --dir).
  arguments$snp <- commandArgs(trailingOnly = TRUE)[2]
  
  # Path to phenotype tab-separated file (absolute, or relative to --dir).
  arguments$pheno <- commandArgs(trailingOnly = TRUE)[3]
  
  # Index of phenotype and optional covariate(s) columns
  arguments$pheno_index <- commandArgs(trailingOnly = TRUE)[4]
  arguments$covar_index <- commandArgs(trailingOnly = TRUE)[5] # Can be "NULL"
  
  # List of genetic algorithm parameters
  params_list <- commandArgs(trailingOnly = TRUE)[6]
  names(params_list) <- c("n_iter","n_ind","n_eli","n_nov","n_chi","n_top",
                          "mutation_rate","crossing_rate")
  # Should have following format (no spacing character between arguments or in ""):
  # list(n_iter,n_ind,n_eli,n_nov,n_chi,n_top,mutation_rate,crossing_rate)
  #   Where all parameters are integers,
  #   except mutation_rate and crossing_rate which are decimal values.
}

source('./r_function/load_inputs.R')
matrices <- input_loading(arguments$snp,
                          arguments$pheno, 
                          arguments$pheno_index,
                          arguments$covar_index)

snp_df <- matrices[[1]]
pheno <- matrices[[2]]
covar <- matrices[[3]]

source('./r_function/generate_GO.R')
source('./r_function/calc_score.R')
source('./r_function/cell_division.R')
source('./r_function/evolve.R')

