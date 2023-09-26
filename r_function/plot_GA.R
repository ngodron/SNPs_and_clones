### Plot things

plot_GA <- function(score_df, conf_int = 0.05) {
library(tidyverse)
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

y_lim <- c(min(score_df$conf_low), max(score_df$conf_high))


if (n_gen > 200) {
  score_df |>
    filter(gen %% (n_gen/100) == 0 | gen < 100) |>
    identity() -> score_df
}
legend_pos <- c(x = 3/4*n_gen,y = max(score_df$conf_high))
out_plot <- 
ggplot(score_df) +
  geom_ribbon(aes(x = gen, ymin = conf_low, ymax = conf_high), alpha = 0.15) +
  geom_smooth(aes(x = gen, y = score), level = 0.99) +
  geom_point(aes(x = gen, y = gen_min, group = gen, colour = n_snps_min), size = 1, shape = 3) +
  # geom_text(aes(x = legend_pos[1], y = legend_pos[2]), label = print_params) +
  scale_colour_viridis_c() +
  ylim(y_lim) +
  theme_bw() +
  xlab("Generations") +
  ylab("Loss") +
  geom_blank()

return(out_plot)
}
