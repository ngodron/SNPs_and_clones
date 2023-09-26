library(ggplot2)

loss_to_BW <- function(scores, penalty_power) {
  # This function takes as input the vector of losses from the fitness function and transforms it 
  # into a vector whose sum is 1, each element is the probability for a genome to be selected by 
  # the Biased Wheel at each sampling.
  # scores is expected to be of same length as genomes with values between 0 and 1, minimal scores being best
  
  # A bigger penalty power favours better solutions more heavily
  penalized_scores <- (1+median(scores)-scores) ** penalty_power
  BW_probabilities <- penalized_scores/sum(penalized_scores)
  return(BW_probabilities)
}

scores_vect <- seq(0, 1, .05)
penal_vect <- 1:10
penalized <- list()

for (i in penal_vect) {
  penalized[[i]] <-
    loss_to_BW(scores = scores_vect, penalty_power = i)
}

names(penalized) <- penal_vect
penalized <- data.frame(penalized)

boxplot(penalized,
        main = "seq(0, 1, .05)",
        xlab = "Penalty power",
        ylab = "BW reproduction probabilities")

