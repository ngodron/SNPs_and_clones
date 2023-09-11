# Generate test dataset
#
library(ape)

set.seed(12102018)
N_indiv <- 5e2
N_feat <- 2e3
N_assoc <- 50
N_homoplasy <- 1e2
N_covar <- 4

# genotypes
genotypes <- 
  matrix(data = 0, nrow = N_feat + N_homoplasy + N_assoc, ncol = N_indiv)
# Maximum degree of association
min_assoc <- 1e-3
max_assoc <- 0.25

pheno_dict <- 
  c('A', 'B')
p_pheno <- c(1,1)
# pheno <- 
#   sample(pheno_dict, N_indiv, replace = TRUE)

tree <-
  rtree(n = N_indiv)
plot(tree)
nodelabels(frame = 'none')
tiplabels(frame = 'cir')

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

# weight for sampling a node to be the first bearing a mutation
p_nodes <- 
  ((tree$edge[,1]- max(tree$edge[,1]) ) / min(tree$edge[,1]))^20

p_nodes <- 
  sapply(X = tree$edge[,1], function(x) length(get_desc(tree = tree, node = x)))

hist(p_nodes)

node_feat <- 
  sapply(X = 1:N_feat, 
         FUN = function(x) sample(tree$edge[,1], size = 1, prob = p_nodes))
summary(node_feat)

for (i in 1:length(node_feat)) {
  curr_node <- node_feat[i]
  curr_children <- get_desc(tree = tree, node = curr_node)
  curr_tips <- curr_children[curr_children <= N_indiv]
  
  # curr_geno <- sample(x = 0:1, 1)
  genotypes[ i, curr_tips ] <- 1#curr_geno
  # genotypes[ i, !(1:ncol(genotypes) %in% curr_tips) ] <- ! curr_geno
}

for (i in N_feat:(N_feat + N_homoplasy)) {
  curr_feat <- i 
  curr_node <- sample(x = unique(as.numeric(tree$edge)), size = 1)
  curr_children <- get_desc(tree = tree, node = curr_node)
  curr_tips <- curr_children[curr_children <= N_indiv]
  # curr_geno <- sample(x = 0:1, 1)
  genotypes[ i, curr_tips ] <- 1#curr_geno
  # genotypes[ i, !(1:ncol(genotypes) %in% curr_tips) ] <- ! curr_geno
}

p_assoc <- runif(n = N_assoc, min = min_assoc, max = max_assoc)

phenotypes <- c()
for (i in 1:N_indiv) {
  curr_geno_asso <- runif(n = N_assoc) <= p_assoc
  genotypes[(1+N_feat+N_homoplasy):(N_feat+N_homoplasy+N_assoc), i] <- curr_geno_asso
  p_pheno <- sum(p_assoc[curr_geno_asso])
  phenotypes[i] <- sample(pheno_dict, size = 1, prob = c(1,p_pheno))
}
table(phenotypes)



