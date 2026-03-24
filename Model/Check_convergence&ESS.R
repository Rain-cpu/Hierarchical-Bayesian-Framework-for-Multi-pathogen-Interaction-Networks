####Results Test####
####By Suyi####
####20250415####


# load package ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(ggplot2)
library(bayesplot)

# 1. Convergence and ESS-----------------------------------------------------------
# 1.1 Posterior distribution ---------------------------------------------------------
#print(fit2, pars = c("mu_or", "w_log", "w_vic", "intercept", "mu_raw", "tau"), probs = c(0.025, 0.5, 0.975))

# summary table
summary_table <- summary(fit2, pars = c("mu_or", "w_log", "w_vic", "intercept", "mu_raw", "tau"))$summary


# 1.2 Serial number to virus name --------------------------------------------------------------
# Step 1: add rownames 
summary_table <- as.data.frame(summary_table) %>%
  rownames_to_column("param")

# Step 2: Extract mu_or[i] & mu_raw[i] & index
mu_params <- summary_table %>%
  filter(str_detect(param, "mu_or\\[|mu_raw\\[")) %>%
  mutate(
    type = ifelse(str_detect(param, "mu_or"), "mu_or", "mu_raw"),
    index = as.integer(str_extract(param, "\\d+"))
  )

# Step 3: link pair_lookup
mu_annotated <- mu_params %>%
  left_join(pair_lookup, by = c("index" = "pair_index")) %>%
  separate(pair_key, into = c("from_id", "to_id"), sep = "_", convert = TRUE) %>%
  left_join(pair_id, by = c("from_id" = "id")) %>%
  rename(from_virus = virus) %>%
  left_join(pair_id, by = c("to_id" = "id")) %>%
  rename(to_virus = virus)

# Step 4: scalar parameters
scalar_params <- summary_table %>%
  filter(!str_detect(param, "mu_or\\[|mu_raw\\[")) %>%
  mutate(
    type = "scalar",
    index = NA_integer_,
    from_virus = NA_character_,
    to_virus = NA_character_
  )

# Step 5: Integration
final_summary_all <- bind_rows(mu_annotated, scalar_params) %>%
  select(param, type, index, from_virus, to_virus, everything())

# Optional
#mu_only <- final_summary_all %>% filter(type != "scalar")
#scalar_only <- final_summary_all %>% filter(type == "scalar")


# 1.3 check list -----------------------------------------------------------------
# if Rhat > 1.01 
final_summary_all[final_summary_all[, "Rhat"] > 1.01, ]

# if ESS<1000
final_summary_all[final_summary_all[, "n_eff"] < 1000, ]


#2. Parameter visualization -------------------------------------------------------------------
post <- rstan::extract(fit2)

# scalar parameters
par(mfrow = c(2, 2))  

hist(post$w_log, breaks = 30, main = "Posterior of w_log", xlab = "w_log", col = "skyblue")
hist(post$w_vic, breaks = 30, main = "Posterior of w_vic", xlab = "w_vic", col = "salmon")
hist(post$tau, breaks = 30, main = "Posterior of tau", xlab = "tau", col = "lightgreen")
hist(post$intercept, breaks = 30, main = "Posterior of intercept", xlab = "intercept", col = "orange")

# vectors
par(mfrow = c(3, 3))  

for (i in 1:17) {
  hist(post$mu_or[, i], breaks = 30,
       main = paste0("Posterior of mu_or[", i, "]"),
       xlab = paste0("mu_or[", i, "]"),
       col = "lightgray")
  
  # Change to a new page after every 9 images
  if (i %% 9 == 0 && i != 17) {
    readline(prompt = "Press [Enter] to continue...")
    par(mfrow = c(3, 3))  
  }
}


# mu_or
posterior_samples <- post$mu_or

# set para names
param_names <- paste0("mu_or[", 1:17, "]")
colnames(posterior_samples) <- param_names

# long format
mcmc_data <- bayesplot::mcmc_trace_data(posterior_samples)

# Mean & CrIs
summary_stats <- mcmc_data %>%
  group_by(parameter) %>%
  summarize(
    mean_value = mean(value),
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975)
  )

# 3. Check if there are three types of data to support the data -------------------------------------------------------
# information table
mu_info <- data.frame(
  mu_index = 1:data_list$N,
  OR_support = FALSE,
  LOG_support = FALSE,
  VCI_support = FALSE,
  support_type = NA_character_
)

# Mark whether each mu_or[i] is referenced by the observation layer.
mu_info$OR_support  <- mu_info$mu_index %in% data_list$or_pair
mu_info$LOG_support <- mu_info$mu_index %in% data_list$log_pair
mu_info$VCI_support <- mu_info$mu_index %in% data_list$vic_pair

# Category support types
mu_info$support_type <- with(mu_info, ifelse(
  OR_support & LOG_support & VCI_support, "fully informed",
  ifelse(OR_support & !LOG_support & !VCI_support, "likelihood only",
         ifelse(!OR_support & (LOG_support | VCI_support), "prior only",
                "partially informed")
  )))

#to viruses name
# Step 1: Merge mu_info and pair_lookup
mu_info_mapped <- mu_info %>%
  left_join(pair_lookup, by = c("mu_index" = "pair_index"))  # 获取 pair_key

# Step 2: Split pair_key into two virus_ids
mu_info_mapped <- mu_info_mapped %>%
  separate(pair_key, into = c("from_id", "to_id"), sep = "_", convert = TRUE)

# Step 3: Link to virus name
mu_info_named <- mu_info_mapped %>%
  left_join(pair_id, by = c("from_id" = "id")) %>%
  rename(from_virus = virus) %>%
  left_join(pair_id, by = c("to_id" = "id")) %>%
  rename(to_virus = virus)

# Step 4: keep the cols need.
mu_info_named_2 <- mu_info_named %>%
  select(mu_index, from_virus, to_virus, OR_support, LOG_support, VCI_support, support_type)

# check
#print(mu_info_named_2)


# 4. edge situation ----------------------------------------------------------------

# A. Sampling positive/negative ratio -------------------------------------------------------------
# Extracting posterior sampling results
mu_samples <- rstan::extract(fit2)$mu_or
intercept_samples <- rstan::extract(fit2)$intercept

# Step 1: set param index
N <- ncol(mu_samples)
mu_list <- paste0("mu_or[", 1:N, "]")

# Step 2: Calculate  posterior difference index
calc_delta <- function(i) {
  delta <- mu_samples[, i] - intercept_samples
  prob_pos <- mean(delta > 0)
  prob_neg <- mean(delta < 0)
  diff <- prob_neg - prob_pos
  return(c(prob_positive = prob_pos, prob_negative = prob_neg, diff = diff))
}

delta_mat <- t(sapply(1:N, calc_delta))  # matrix: [N x 3]
delta_df <- as.data.frame(delta_mat) %>%
  mutate(index = 1:N, param = paste0("mu_or[", index, "]"))

# Step 3: add from_virus, to_virus
delta_df <- delta_df %>%
  left_join(pair_lookup, by = c("index" = "pair_index")) %>%
  separate(pair_key, into = c("from_id", "to_id"), sep = "_", convert = TRUE) %>%
  left_join(pair_id, by = c("from_id" = "id")) %>%
  rename(from_virus = virus) %>%
  left_join(pair_id, by = c("to_id" = "id")) %>%
  rename(to_virus = virus) %>%
  select(param, index, from_virus, to_virus, prob_positive, prob_negative, diff)

#print(delta_df)



# B. intercept ----------------------------------------------------------------
# Step 1:  mu_or
mu_plot_data <- final_summary_all %>%
  filter(type == "mu_or") %>%
  mutate(label = paste0(from_virus, " → ", to_virus))

# Step 2: add intercept
intercept_val <- final_summary_all %>%
  filter(param == "intercept") %>%
  pull(mean)

mu_plot_data <- mu_plot_data %>%
  mutate(
    covers_intercept = intercept_val >= `2.5%` & intercept_val <= `97.5%`,
    relation = case_when(
      `97.5%` < intercept_val ~ "Higher than CrI",
      `2.5%` > intercept_val ~ "Lower than CrI",
      TRUE ~ "Inside CrI"
    )
  )

#write.csv(mu_plot_data, "test.csv")

# Step 3: plot

library(ggplot2)

ggplot(mu_plot_data, aes(x = mean, y = reorder(label, mean))) +
  geom_errorbarh(aes(xmin = `2.5%`, xmax = `97.5%`, color = relation), 
                 height = 0.25, size = 0.8) +
  geom_point(aes(color = relation), size = 3) +
  geom_vline(xintercept = intercept_val, linetype = "dashed", 
             color = "black", size = 0.7) +
  scale_color_manual(values = c(
    "Inside CrI" = "#999999",
    "Higher than CrI" = "#3161a1", 
    "Lower than CrI" = "#fa8d55"    
  )) +
  labs(
    x = "Posterior mean",
    y = "Virus pairs",
    title = "Posterior means of VII with 95% CrI",
    color = "Intercept relation"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0, face = "bold"),
    axis.text.y  = element_text(size = 15, face = "italic"), 
    axis.text.x = element_text(size = 15),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )
# 9. output -----------------------------------------------------------------
