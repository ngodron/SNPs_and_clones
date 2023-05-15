sapply(X = list.files(path = './r_function/list_version/', 
      full.names = TRUE), 
      FUN = source)

first_gen <- 
  generate_G0(n_snps = 50, n_indiv = 20, p = 1/10)

curr_gen <- 
  evolve(genomes = first_gen, 
         mut_rate = 0.1, 
         conjug_rate = 0.5, 
         min_mut = 1)

cell_division(genomes = curr_gen, 
              scores = , 
              n_best = , 
              n_child = , 
              n_elite = , 
              n_novel = , 
              mu = , 
              cr = )