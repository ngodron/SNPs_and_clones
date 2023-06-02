library(profvis)
library(pryr)
library(tidyverse)
# set.seed(123456)
# base_dir <- '~/Kevin_These/Current_File/'

base_dir <- '../../THESE_KLA/Current_File/'

curr_unique_name <- 
  c(letters, LETTERS, 0:9)
curr_unique_name <- sample(curr_unique_name, 24)
curr_unique_name <- paste0(curr_unique_name, collapse = "")
curr_unique_name <- 
  paste0('../output/',format(Sys.time(), "%Y_%m_%d-%H_%M_%S-"),'-genet-', 
         curr_unique_name, '/', collapse = "")

dir.create(path = curr_unique_name)


snp_df <- 
  read.delim(file = paste0(base_dir, '230516_matrix_SNP_gene_roary_594.csv'))
strain_names <- snp_df[,1]

metadata <- 
  read.delim(file = paste0(base_dir, 'database_594.csv'))

pheno <- 
  merge(y = metadata, x = snp_df[,1:2], by.y = 'Reads', by.x = 'X') 
pheno <- pheno$Disease != 'sputum'


prop_priors <- sapply(X = unique(pheno), function(x) 1 - sum(pheno == x) / length(pheno))
# prop_priors <- sapply(X = unique(pheno), function(x) sum(pheno == x) / length(pheno))
prop_priors <- c(0.5,0.5)
  
weights_df <- 
  read.csv(file = '../../THESE_KLA/Current_File/2023_05_31-homoplasy.tsv', 
           sep = '\t')
w_in_snp <- weights_df$SNP %in% names(snp_df)
snp_in_w <- names(snp_df) %in% weights_df$SNP
# Both should be 100% TRUE !
table(w_in_snp)
table(snp_in_w)
names(snp_df)[!snp_in_w]
weights_df <- weights_df[w_in_snp, ]
snp_df <- snp_df[,snp_in_w]
weights_df <-
  weights_df[match(names(snp_df), weights_df$SNP), ]

table(weights_df$SNP == names(snp_df))
all_weights <- weights_df$homo_tot
names(all_weights) <- weights_df$SNP

all_lin <- sort(unique(metadata$Lineage))
n_lin <- length(all_lin)
covar <- matrix(nrow = nrow(metadata), ncol = n_lin)
colnames(covar) <- all_lin


for (i in 1:n_lin) {
  curr_lin <- all_lin[i]
  covar[,i] <- metadata$Lineage == curr_lin
}

# source('./r_function/generate_GO.R')
# source('./r_function/calc_score.R')
# source('./r_function/cell_division.R')
# source('./r_function/evolve.R')
sapply(X = list.files(path = './r_function/list_version/', 
                      full.names = TRUE), 
       FUN = source)

# General parameters
## Number of generations
n_iter <- 1e2
## Number of individuals
n_ind <- 1e2

# Population parameters
n_eli <- ceiling(n_ind / 20) # Count of elites
n_nov <- ceiling(n_ind / 20) # Count of novel (random) individuals
n_chi <- 4 # Number of children per "top genome"

n_top <- n_ind/n_chi - (n_eli/n_chi) - (n_nov/n_chi)
n_top <- ceiling(n_top) 

n_ind <- n_top * n_chi + (n_eli + n_nov) # Total count of individuals

mutation_rate <- 1e-3
crossing_rate <- 0.5

print_params <- 
  paste0('n_ind = ', n_ind, '\nmu = ', mutation_rate)

all_gen <- vector('list', n_iter)

curr_gen <- 
  generate_G0(n_snps = ncol(snp_df), 
              n_indiv = n_ind, 
              p = runif(n = ncol(snp_df), min = 0, max = mutation_rate))
sapply(curr_gen, sum)

score_list <- vector(mode = 'list', length = n_iter)
model_list <- vector(mode = 'list', length = n_iter)
diversity <- rep(NA_real_, n_iter)

# library(doParallel)
# doParallel::registerDoParallel(cores = 48)
# doParallel::registerDoParallel(cores = 1)

# library("profmem")
# memory_usage <- profmem({

system.time({
# profvis({
  for (i in 1:n_iter) {
    start_gen <- Sys.time()
    cat('Generation ', i, '/', n_iter, '\n')
    # set_snps <- which(curr_gen == 1, arr.ind = TRUE)
    # set_snps <- set_snps[order(set_snps[,1]), ]
    
    # for (j in 1:n_ind) {
    #   snp_list[[j]] <- 
    #     set_snps[set_snps[,1] == j, 2, drop = TRUE]
    # }
    all_gen[[i]] <- 
      sapply(curr_gen, function(x) which(x == 1, arr.ind = TRUE))
    all_gen[[i]] <- sapply(all_gen[[i]], function(x) names(snp_df)[x])
    
    curr_scores_models <- 
      calc_score(genomes = curr_gen, 
                       snps = snp_df, 
                       phenotype = pheno, 
                       fitness = decision_tree_fitness, 
                       covars = covar, 
                       weights = all_weights)
    
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
                    mu = mutation_rate, 
                    cr = crossing_rate)
    # print(length(next_gen))
    curr_gen <- next_gen
    end_gen <- Sys.time()
    gen_duration <- round(end_gen - start_gen, digits = 1)
    print(gen_duration)
  }
})

# Rprofmem(NULL) / 0.99

score_df <-
  data.frame(gen = rep(1:length(score_list), each = n_ind),
             score = unlist(score_list),
             n_snps = unlist(sapply(all_gen, FUN = function(x) {
               sapply(x, length)}
               , simplify = FALSE)),
             diversity = rep(diversity, each= n_ind))



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

if (n_iter > 200) {
  score_df |>
    filter(gen %% (n_iter/100) == 0 | gen < 100) |>
    identity() -> score_df
}

# write.table(x = )
# save(list = c("all_gen", 'model_list', 'score_list', 'score_df'), file = paste0(curr_unique_name,'all_stuff.stuff', collapse = ''))

# load('../output/2023_05_18-14_44_38---genet-Ie1v6Dt2JH4cwLOPWASU8sfFall_gen')

# 
# ggplot(score_df, mapping = aes(x = gen, y = diversity)) +
#   geom_point(aes(x = gen, y = diversity, colour = gen_min)) +
#   geom_smooth() +
#   theme_bw()
# 
# ggplot(score_df) +
#   geom_point(aes(x = 1 - gen_min,  y = diversity)) +
#   geom_smooth(aes(x = 1 - gen_min,  y = diversity)) +
#   theme_bw()
# 
# # ggplot(score_df) +
# #   geom_boxplot(aes(x = gen, y = n_snps, group = gen)) +
# #   theme_bw()
# 
legend_pos <- c(x = 3/4*n_iter,y = max(score_df$conf_high))
ggplot(score_df) +
  geom_ribbon(aes(x = gen, ymin = conf_low, ymax = conf_high), alpha = 0.05) +
  geom_smooth(aes(x = gen, y = score), level = 0.99) +
  geom_point(aes(x = gen, y = gen_min, group = gen, colour = n_snps_min), size = 1, shape = 3) +
  geom_text(aes(x = legend_pos[1], y = legend_pos[2]), label = print_params) +
  scale_colour_viridis_c() +
  theme_bw() +
  geom_blank()

rattle::fancyRpartPlot(model_list[[n_iter]], cex = 0.8)
rpart.plot::rpart.plot(model_list[[n_iter]], cex = 0.8)

tmp_df <- data.frame(pheno,  snp_df, covar)

system.time({
  model <- rpart::rpart(formula = pheno ~ .,
                        data = tmp_df,
                        minbucket = 10,
                        method = 'class', 
                        cost = 1/c(all_weights, rep(mean(all_weights), ncol(covar)))
                        )
})

caret::confusionMatrix(
  data = as.factor(predict(model)[, 1] <=0.5),
  reference = factor(pheno))

caret::confusionMatrix(
  data = as.factor(predict(model_list[[n_iter]])[, 1] <=0.5),
  reference = factor(pheno))


