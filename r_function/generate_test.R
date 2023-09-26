# Generate test dataset ----
#

library(ape)

set.seed(12102018)
N_indiv <- 2e3
N_feat <- 5e3



# genotypes ----
genotypes <- 
  matrix(data = sample(c(0), N_feat * N_indiv, replace = TRUE), 
         nrow = N_feat, 
         ncol = N_indiv)

pheno_dict <- 
  c('A', 'B')

tree <-
  rtree(n = N_indiv)
int_nodes <- 
  tree$edge[,1]
plot(tree, show.tip.label = FALSE, show.node.label = FALSE)


get_desc <- function(tree, node, curr=NULL){
  if (node > max(tree$edge)) stop('Too far !') 
  if(is.null(curr)) curr <- vector()
  daughters <- tree$edge[which(tree$edge[,1]==node),2]
  curr<-c(curr,daughters)
  # Because nodes start their numbering at length(tip) + 1
  w <- which(daughters >= length(tree$tip.label))
  if(length(w)>0) for(i in 1:length(w))
    curr<-get_desc(tree,daughters[w[i]],curr)
  return(curr)
}

get_dist_root <- function(tree, node, dist = 0) {
  curr_index <- 
    which(tree$edge[,2] == node)
  if(length(curr_index) > 0 ) {
    curr_parent <- 
      tree$edge[curr_index, 1]
    dist <- 
      dist + tree$edge.length[curr_index]
    dist <- 
      get_dist_root(tree = tree, node = curr_parent,  dist = dist)
  }
  return(dist)
}


# Random attribution weighted by dist to root----
## all dist to root ----
all_dtr <- 
  sapply(tree$edge[,2], function(x) get_dist_root(tree, x))
summary(all_dtr)
## presence / absence ----
mut_nodes <- 
  sample(x = tree$edge[,2], 
         size = N_feat, 
         prob = NULL,
         replace = TRUE)

hist(mut_nodes, breaks = 1e2)

mut_children <- 
    sapply(X = mut_nodes, function(x) get_desc(tree = tree, node = x))

## introducing homoplasy ----
n_homop <- rpois(n = length(mut_children), lambda = 2)
for (i in 1:length(mut_children)) {
  if (n_homop[i] > 0){
    to_remove <- 
      unique(sample(x = 1:length(mut_children[[i]]), 
                    size = n_homop[i], replace = TRUE))
  
  if (length(to_remove) < length(mut_children[[i]])) {
    mut_children[[i]] <- mut_children[[i]][-to_remove]
  }
  }
}

mut_tips <- 
   sapply(X = mut_children, function(x) x[x <= N_indiv]) # N_indiv tips, the rest are internal nodes
for (i in 1:N_feat) genotypes[i, mut_tips[[i]] ] <- ! genotypes[i, mut_tips[[i]]]

rownames(genotypes) <- 
  paste0('FEAT_', 1:nrow(genotypes))


## Calc homoplasy 
# source(file = 'r_function/homoplasy_counter.R')
# message('homoplasy !')
# system.time({
#   homoplasy <- 
#     multi_homo(feat_matrix = t(genotypes), 
#            tree = tree)
# })

# Feat effect on phenotype ----
## beta ----
all_beta <-
  rgamma(n = N_feat, shape = 12, scale = 1e1/N_feat)

### centering and folding the distribution ----
all_beta <-
  (all_beta - mean(all_beta))

## homoplasic SNP are more strongly associated ----
all_beta <- 
  all_beta * n_homop

beta_sum <- 
  apply(X = genotypes, MARGIN = 2, FUN = function(x) sum(all_beta[x == 1]))  
  # rnorm(n = N_indiv, mean = 0, sd = 1e-1)
hist(beta_sum)

phenotypes <- ifelse(beta_sum <= 0, 'A', 'B')
table(phenotypes)

# Mutation associated with phenotype
# p_appear <-
#   runif(n = N_assoc, min = 0, max = 1)
# summary(p_appear)
# p_assoc <- 
#   runif(n = N_assoc, min = min_assoc, max = max_assoc)
# summary(p_assoc)
# 
# phenotypes <- c()
# all_p_pheno <- vector(mode = 'integer', length(N_indiv))
# for (i in 1:N_indiv) {
#   curr_geno_asso <- runif(n = N_assoc) <= p_appear
#   genotypes[(1+N_feat):(N_feat + N_assoc), i] <- curr_geno_asso
#   p_pheno <- sum(p_assoc[curr_geno_asso])
#   all_p_pheno[i] <- p_pheno
#   phenotypes[i] <- sample(pheno_dict, size = 1, prob = c(median(p_pheno), p_pheno))
# }
# summary(all_p_pheno)


# Adding the perfect predictor
# genotypes <- rbind(genotypes, 0)
# have_perfect_SNP <-
#   sample(x = 1:ncol(genotypes), size = 50 ,replace = FALSE)
# genotypes[nrow(genotypes), have_perfect_SNP] <- 1
# phenotypes[have_perfect_SNP] <- 'B'
# # phenotypes[!have_perfect_SNP] <- 'A'
# fisher.test(table(phenotypes, genotypes[nrow(genotypes), ]))
# 
# table(phenotypes)
# table(phenotypes, genotypes[nrow(genotypes), ])



plot(tree, show.tip.label = FALSE)
tiplabels(text = 'o', 
          frame = 'none',
          col = as.numeric(as.factor(phenotypes)))

# Population structure covariables ----
pca <- 
  prcomp(genotypes, rank. = 2)

clusters <- 
  dbscan::dbscan(x = pca$rotation[,1:2], minPts = 5, eps = 0.004)$cluster

clusters <- 
  data.frame(data = as.factor(letters[clusters + 1]))
colnames(clusters) <- 'cluster'

plot(pca$rotation[,1], pca$rotation[,2], 
     col = clusters + 1, 
     pch = as.numeric(as.factor(phenotypes)) + 1)

library(ggplot2)
ggplot(data = data.frame(pca$rotation, cluster = as.factor(clusters$cluster), phenotype = phenotypes)) +
  geom_point(aes(x = PC1, y = PC2, colour = cluster, shape = phenotype), 
             size  = 3, alpha = 0.2 ) +
  scale_color_viridis_d()+
  scale_fill_viridis_d()+
  theme_bw()


genotypes <- t(genotypes)

# genotypes <- 
#   genotypes[ , order(all_beta)]


# rm(list = all_things)

# Training and testing ----
data_size <- nrow(genotypes)
half_data <- floor(data_size/2)
to_train <- sample(x = data_size, size = half_data)
to_test <- (1:data_size)[!(1:data_size %in% to_train)]

# Split costs ----
costs <- 
  n_homop * 10
costs[costs == 0] <- 1e-1
costs <- 
  1/costs

costs_homo <-
 c(costs, rep(median(costs), ncol(clusters)))
names(costs_homo) <- 
  c(colnames(genotypes), colnames(clusters)) 

# Genetic algorithm ! ----
source(file = 'r_function/gen_algo.R')
n_gen <- 1e1
system.time({
glou <-
  gen_algo(snp_matrix = genotypes[to_train,], 
           pheno_matrix = phenotypes[to_train], 
           covar_matrix = clusters[to_train, ,drop = FALSE], 
           costs = costs_homo, 
           n_gen = n_gen, 
           n_iter = 1, 
           n_ind = 1e2, 
           n_eli = 2, 
           n_nov = 10, 
           n_chi = 4, 
           mu = 1e-3, 
           cr = 0.5, 
           curr_gen = NULL, 
           fitness_fun = decision_tree_fitness)
})


# Ploting time ! ----
source('r_function/plot_GA.R')
plot_GA(score_df = glou[[5]])
rpart.plot::rpart.plot(glou[[3]][[n_gen]], cex = 0.8)



# Model evaluation: ----
## tautologic ----
predictions_tauto <- 
  predict(glou[[3]][[n_gen]])
predicted_tauto_pheno <- 
  ifelse(predictions_tauto[,1] >= 0.5, 'A', 'B')
caret::confusionMatrix(data = factor(predicted_tauto_pheno), 
                       reference = factor(phenotypes[to_train]))

## Testing dataset ----
predictions <- 
  predict(glou[[3]][[n_gen]], 
          newdata = data.frame(genotypes[to_test, ], clusters[to_test, ,drop = FALSE]))
predicted_pheno <- 
  ifelse(predictions[,1] >= 0.5, 'A', 'B')
caret::confusionMatrix(data = factor(predicted_pheno), 
                       reference = factor(phenotypes[to_test]))
