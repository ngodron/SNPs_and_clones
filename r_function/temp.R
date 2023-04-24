base_dir <- '../../THESE_KLA/'
# base_dir <- '~/Kevin_These/'
# setwd("/home/nicolas/2023/SNPs_and_clones")


snp_df <- 
  read.delim(file = paste0(base_dir,'Current_File/230331_matrix_snp_gene_613.csv'))

metadata <- 
  read.delim(file = paste0(base_dir,'/database/database_613.csv'))

pheno <- 
  merge(y = metadata, x = snp_df[,1:2], by.y = 'Reads', by.x = 'X') 

pheno <- pheno$Disease != 'sputum'
strain_names <- snp_df[,1]
snp_df <- snp_df[,-1]
source('./r_function/generate_GO.R')
source('./r_function/calc_score.R')
source('./r_function/cell_division.R')
source('./r_function/evolve.R')

n_iter <- 1e1
n_ind <- 1e2
n_eli <- 20
n_chi <- 4
n_be <- n_ind/n_chi - n_eli/n_chi
n_be <- ceiling(n_be)
n_ind <- n_be * n_chi + n_eli

all_gen <- vector('list', n_iter)

curr_gen <- 
  generate_G0(n_snps = ncol(snp_df), 
              n_indiv = n_ind, 
              p = runif(n = ncol(snp_df), min = 0, max = 0.0001))
rowSums(curr_gen)

score_list <- vector(mode = 'list', length = n_iter)

system.time({
for (i in 1:n_iter) {
  cat('Generation ', i - 1, '\n')
  print(summary(rowSums(curr_gen)))
  all_gen[[i]] <- curr_gen
  curr_scores <- 
    calc_score(genomes = curr_gen, 
               snps = snp_df, 
               phenotype = pheno, 
               fitness = glm_fitness_mock)
  print(summary(curr_scores))
  score_list[[i]]<- c(curr_scores)
  
  # if (min(curr_scores) < 0.0) {
  #   print('tadaaaa')
  #   all_gen <- all_gen[1:i]
  #   score_list <- score_list[1:i]
  #   break()
  # }
  
  next_gen <- 
    cell_division(genomes = curr_gen, 
                  scores = curr_scores, 
                  n_best = n_be, 
                  n_child = n_chi, 
                  n_elite = n_eli, mu = 0.0001, cr = 0.8)
  curr_gen <- next_gen
}
})

score_df <- 
  data.frame(gen = rep(1:n_iter, each = n_ind), 
             score = unlist(score_list),
             n_snps = unlist(sapply(all_gen, FUN = rowSums, simplify = FALSE)))


require(tidyverse)
conf_int <- 0.05
score_df |> 
  group_by(gen) |> 
  mutate(gen_min = min(score, na.rm = TRUE)) |>
  mutate(is_min = score == gen_min) |>
  mutate(n_snps_min = n_snps[score == gen_min][1]) |>
  mutate(conf_low = quantile(score, probs = conf_int)) |> 
  mutate(conf_high = quantile(score, probs = 1-conf_int)) |> 
  identity() -> score_df

ggplot(score_df) +
  geom_ribbon(aes(x = gen, ymin = conf_low, ymax = conf_high), alpha = 0.05) +
  geom_smooth(aes(x = gen, y = score), level = 0.99) +
  # geom_boxplot(aes(x = gen, y = score, group = gen)) +
  geom_point(aes(x = gen, y = gen_min, group = gen, colour = n_snps_min), size = 1, shape = 3) +
  
  theme_bw() +
  geom_blank()
