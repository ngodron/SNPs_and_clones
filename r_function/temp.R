library(profvis)
# base_dir <- '../../THESE_KLA/'
base_dir <- '~/Kevin_These/'
setwd("/home/nicolas/2023/SNPs_and_clones")


snp_df <- 
  read.delim(file = paste0(base_dir,'Current_File/230331_matrix_snp_gene_613.csv'))

metadata <- 
  read.delim(file = paste0(base_dir,'/database/database_613.csv'))

pheno <- 
  merge(y = metadata, x = snp_df[,1:2], by.y = 'Reads', by.x = 'X') 


pheno <- pheno$Disease != 'sputum'

strain_names <- snp_df[,1]
snp_df <- snp_df[,-1]

all_lin <- sort(unique(metadata$Lineage))
n_lin <- length(all_lin)
covar <- matrix(nrow = nrow(metadata), ncol = n_lin)
colnames(covar) <- all_lin
for (i in 1:n_lin) {
  curr_lin <- all_lin[i]
  covar[,i] <- metadata$Lineage == curr_lin
}


source('./r_function/generate_GO.R')
source('./r_function/calc_score.R')
source('./r_function/cell_division.R')
source('./r_function/evolve.R')

genomes_diversity <- function(genomes) {
  n_by_cols <- colSums(genomes)
  genomes <- genomes[ ,n_by_cols >= 1, drop = FALSE]
  prop <- colSums(genomes) / nrow(genomes)
  out <- 1 - abs(mean(prop) - 0.5) * 2
  return(out)
}

# General parameters
n_iter <- 1e2 # Number of generations
n_ind <- 1e2 # Number of individuals

# Population parameters
n_eli <- ceiling(n_ind / 20) # Count of elites
n_nov <- ceiling(n_ind / 20) # Count of novel (random) individuals
n_chi <- 4 # Number of children per "top genome"

n_top <- n_ind/n_chi - (n_eli/n_chi) - (n_nov/n_chi)
n_top <- ceiling(n_top) 

n_ind <- n_top * n_chi + (n_eli + n_nov) # Total count of individuals

mutation_rate <- 1e-3

params <- 
  paste0('n_ind = ', n_ind, '\nmu = ', mutation_rate)

all_gen <- vector('list', n_iter)

curr_gen <- 
  generate_G0(n_snps = ncol(snp_df), 
              n_indiv = n_ind, 
              p = runif(n = ncol(snp_df), min = 0, max = mutation_rate))
rowSums(curr_gen)

score_list <- vector(mode = 'list', length = n_iter)
model_list <- vector(mode = 'list', length = n_iter)
diversity <- rep(NA_real_, n_iter)

# library(doParallel)
# doParallel::registerDoParallel(cores = 48)
# doParallel::registerDoParallel(cores = 1)

# library("profmem")
# memory_usage <- profmem({
# system.time({
profvis({
  
  for (i in 1:n_iter) {
    cat('Generation ', i, '/', n_iter, '\n')
    set_snps <- which(curr_gen == 1, arr.ind = TRUE)
    set_snps <- set_snps[order(set_snps[,1]), ]
    snp_list <- vector(mode = 'list', length = n_ind)
    for (j in 1:n_ind) {
      snp_list[[j]] <- 
        set_snps[set_snps[,1] == j, 2, drop = TRUE]
    }
    all_gen[[i]] <- snp_list
    print(summary(sapply(all_gen[[i]], length)))
    
    curr_scores_models <- 
      calc_score_nopar(genomes = curr_gen, 
                       snps = snp_df, 
                       phenotype = pheno, 
                       fitness = decision_tree_fitness, 
                       covars = covar)
    
    curr_scores <- 
      unlist(lapply(curr_scores_models, function(x) x[[1]]))
    curr_models <- 
      (lapply(curr_scores_models, function(x) x[[2]]))
    model_list[[i]] <-
      curr_models[[which(curr_scores == min(curr_scores))[1]]] 
    
    rm(curr_scores_models)
    print(summary(curr_scores))
    score_list[[i]]<- c(curr_scores)
    diversity[i] <- genomes_diversity(curr_gen)
    
    # if (min(curr_scores) < 0.0) {generate_G0
    #   print('tadaaaa')
    #   all_gen <- all_gen[1:i]
    #   score_list <- score_list[1:i]
    #   break()
    # }
    
    next_gen <- 
      cell_division(genomes = curr_gen, 
                    scores = curr_scores, 
                    n_best = n_top, 
                    n_child = n_chi, 
                    n_elite = n_eli, 
                    n_novel = n_nov,
                    mu = mutation_rate, cr = 0)
    rm(list = ls()[grep(pattern = 'curr_.*', x = ls())])
    curr_gen <- next_gen
  }
})

# Rprofmem(NULL)

score_df <- 
  data.frame(gen = rep(1:length(score_list), each = n_ind), 
             score = unlist(score_list),
             n_snps = unlist(sapply(all_gen, FUN = function(x) {
               sapply(x, length)}
               , simplify = FALSE)),
             diversity = rep(diversity, each= n_ind))


require(tidyverse)
conf_int <- 0.05
score_df |> 
  group_by(gen) |>
  mutate(gen_min = min(score, na.rm = TRUE)) |>
  mutate(is_min = score == gen_min) |>
  mutate(n_snps_min = min(n_snps[score == gen_min])[1]) |>
  mutate(conf_low = quantile(score, probs = conf_int)) |> 
  mutate(conf_high = quantile(score, probs = 1-conf_int)) |>
  ungroup() |> 
  filter(gen >= 0) |> 
  identity() -> score_df

if (n_iter > 200){
  score_df |> 
    filter(gen %% (n_iter/100) == 0 | gen < 100) |>
    identity() -> score_df
}

ggplot(score_df, mapping = aes(x = gen, y = diversity, colour = gen_min)) +
  geom_point() +
  geom_smooth() +
  theme_bw()

ggplot(score_df) +
  geom_point(aes(x = 1 - gen_min,  y = diversity)) +
  geom_smooth(aes(x = 1 - gen_min,  y = diversity)) +
  theme_bw()

# ggplot(score_df) +
#   geom_boxplot(aes(x = gen, y = n_snps, group = gen)) +
#   theme_bw()

legend_pos <- c(x = 3/4*n_iter,y = max(score_df$conf_high))
ggplot(score_df) +
  geom_ribbon(aes(x = gen, ymin = conf_low, ymax = conf_high), alpha = 0.05) +
  geom_smooth(aes(x = gen, y = score), level = 0.99) +
  geom_point(aes(x = gen, y = gen_min, group = gen, colour = n_snps_min), size = 1, shape = 3) +
  geom_text(aes(x = legend_pos[1], y = legend_pos[2]), label = params) +
  scale_colour_viridis_c() +
  theme_bw() +
  geom_blank()


