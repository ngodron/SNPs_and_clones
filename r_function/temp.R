base_dir <- '../../THESE_KLA/'

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


curr_gen <- 
  generate_G0(n_snps = 1e2, n_indiv = 1e1, p = 0.1)

n_iter <- 1e3

score_list <- vector(mode = 'list', length = n_iter)

for (i in 1:n_iter) {
  print(i)
  curr_scores <- 
    calc_score(genomes = curr_gen, 
               snps = snp_df, 
               phenotype = pheno, 
               fitness = rand_fitness_mock)
  score_list[[i]]<- c(curr_scores)
  
  next_gen <- 
    cell_division(genomes = curr_gen, 
                  scores = curr_scores, 
                  n_best = 4, 
                  n_child = 2, 
                  n_elite = 2, mu = 1, cr = 0)
  
}

score_df <- do.call(rbind, score_list)
plot(y = apply(X = score_df, MARGIN = 1, FUN = mean), x = 1:n_iter)
