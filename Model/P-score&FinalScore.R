#### VCI MNA + CI-adjusted P-score Analysis ####
#### By Suyi, 20250423 ####

# ----------------- Step 1:  P-score matrix -----------------
n_pairs <- ncol(mu_samples)
p_score_matrix <- matrix(NA, nrow = n_pairs, ncol = n_pairs)

for (i in 1:n_pairs) {
  for (j in 1:n_pairs) {
    if (i != j) {
      p_score_matrix[i, j] <- mean(mu_samples[, i] < mu_samples[, j])  
    }
  }
}

rownames(p_score_matrix) <- colnames(p_score_matrix) <- final_summary_all$param[1:n_pairs]
p_scores <- rowMeans(p_score_matrix, na.rm = TRUE)

# ----------------- Step 2: result_table -----------------
mu_rows <- final_summary_all$type == "mu_or"
result_table <- final_summary_all[mu_rows, ]
result_table$p_score <- p_scores

# ----------------- Step 3: SURCA-----------------
rank_matrix <- t(apply(mu_samples, 1, function(row) rank(-row, ties.method = "average")))  # descending
n_pairs <- ncol(mu_samples)
rank_probs <- matrix(0, nrow = n_pairs, ncol = n_pairs)
for (i in 1:n_pairs) {
  for (r in 1:n_pairs) {
    rank_probs[i, r] <- mean(rank_matrix[, i] == r)
  }
}
sucra_values <- apply(rank_probs, 1, function(p) {
  cumsum_probs <- cumsum(p)
  1 - sum(cumsum_probs[1:(length(p) - 1)]) / (length(p) - 1)
})
result_table$sucra_bayesian <- sucra_values

# ----------------- Step 4:FinalScore -----------------
intercept_val <- mean(intercept_samples) 

result_table$ci_position <- with(result_table, ifelse(
  `97.5%` < intercept_val, "fully_below",
  ifelse(`2.5%` > intercept_val, "fully_above", "overlapping")
))

result_table$directional_boost <- ifelse(
  result_table$ci_position == "fully_below", 1.1,
  ifelse(result_table$ci_position == "fully_above", 0.9, 1) #Sensitivity analysis
)

# FinalScore
result_table$final_score <- result_table$p_score * result_table$directional_boost

final_df <- data.frame(index = seq_along(result_table$final_score), final_score = result_table$final_score)
final_sorted <- final_df[order(final_df$final_score, decreasing = TRUE), ]

final_sorted$sucra <- sapply(seq_len(nrow(final_sorted)), function(i) {
  if (i == 1) return(1)
  sum(final_sorted$final_score[1:(i - 1)]) / (i - 1) 
})

sucra_final <- rep(NA, nrow(result_table))
sucra_final[final_sorted$index] <- final_sorted$sucra
result_table$sucra <- sucra_final

# ----------------- Step 5: Rank of FinalScore -----------------
library(ggplot2)
library(reshape2)

result_table$label <- paste(result_table$from_virus, "→", result_table$to_virus)

ggplot(result_table, aes(x = reorder(label, sucra), y = sucra, fill = sucra)) +
  geom_col(width = 0.7) +

  scale_fill_gradient(low = "#c8d5e5", high = "#3a5d9d", name = "") + #c6dbef 08306b
  
  coord_flip() +
  labs(
    title = "FinalScore Order",
    x = "",
    y = ""
  ) +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_text(size = 14, face = "italic"), 
    axis.text.x = element_text(size = 14)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.15)))

# ----------------- Step 6: P-score Heat map-----------------
labels <- paste(result_table$from_virus, "→", result_table$to_virus)
rownames(p_score_matrix) <- colnames(p_score_matrix) <- labels

p_score_df <- melt(p_score_matrix, varnames = c("Row", "Col"), value.name = "P_score", na.rm = TRUE)

ggplot(p_score_df, aes(x = Col, y = Row, fill = P_score)) +
  geom_tile(color = "white", size = 0.01) +
  scale_fill_gradientn(
    colors = c("#de884b", "#f7f7f7", "#475eac"),  # "#e66101", "#f7f7f7", "#5e3c99" fff8ee 
    limits = c(0, 1), 
    name = ""
  ) 
