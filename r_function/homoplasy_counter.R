require(ape)
require(tidytree)


homo_counter <- function(tree, feat) {
  if(sum(feat) == 0 ) {
    message('Uncool feature without 1 !')
    return(list())
  } else if (sum(feat) == length(feat)) {
    message('Uncool feature without 0 !')
    return(list())
  }
  
  feat_ac <- 
    ace(x = feat, phy = tree, type = 'discrete', model = 'ARD')
  feat_full <-
    ifelse(feat_ac$lik.anc[,1] > 0.5, 
         as.numeric(colnames(feat_ac$lik.anc)[1]),
         as.numeric(colnames(feat_ac$lik.anc)[2]))
  feat_full <- c(feat, feat_full)
  tree <- as_tibble(tree)
  tree$feat <- feat_full
  root <- min(tree$parent)
  homo_list <-  list()
  homo_list <- explore(node = root, 
                       exp_tree = tree, 
                       homo_list = homo_list, 
                       i = 1)
  # homo_list <- homo_list[!sapply(homo_list, is.null)]
  return(homo_list)
}

explore <- function(exp_tree, node, homo_list, i) {
  
  curr_child <- exp_tree$node[exp_tree$parent == node]
  curr_child <- curr_child[curr_child != node]
  for (child in curr_child) {
    if (exp_tree$feat[node] != exp_tree$feat[child] ) {
      curr_feat <- exp_tree$feat[child]
      homo_list[[i]] <- c(child, curr_feat)
      i <- i + 1
    }
  }
  for (child in curr_child) {
    homo_list <- explore(exp_tree = exp_tree, node = child, homo_list = homo_list, i = i)
  }
  return(homo_list)
}


multi_homo <- function(feat_matrix, tree) {
  # feat in cols
  i <- 1
  homo_list <- 
  apply(X = feat_matrix,
        MARGIN = 2,
        FUN = function(x) {cat(i, '\n');i <<- i + 1;homo_counter(tree, feat = x)})
  homo_names <- names(homo_list)
  homo_length <- sapply(homo_names, function(x) length(homo_list[[x]]))
  homo_df <- data.frame(Feat = homo_names, N_homo = homo_length)
  return(list(homo_df, homo_list))
}



