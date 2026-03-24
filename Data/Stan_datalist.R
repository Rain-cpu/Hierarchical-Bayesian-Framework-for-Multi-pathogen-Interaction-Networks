####Stan Data Structures####
####By Suyi####
####20250413####


# 1. Number of virus pairs ---------------------------------------------------
N <- max(stan_input$pair_index)


# 2. OR layer input -------------------------------------------------------
or_layer <- stan_input %>% filter(layer == "OR")
K_or <- nrow(or_layer)
or_data <- list(
  K_or = K_or,
  theta_or = or_layer$log_value_scaled,
  sigma_or = or_layer$sigma,
  or_pair = or_layer$pair_index
)


# 3. LOG layer input ------------------------------------------------------
log_layer <- stan_input %>% filter(layer == "LOG")
K_log <- nrow(log_layer)
log_data <- list(
  K_log = K_log,
  theta_log = log_layer$log_value_scaled,
  sigma_log = log_layer$sigma,
  log_pair = log_layer$pair_index
)


# 4. VIC layer input ------------------------------------------------------
vic_layer <- stan_input %>% filter(layer == "VIC")
K_vic <- nrow(vic_layer)
vic_data <- list(
  K_vic = K_vic,
  theta_vic = vic_layer$log_value_scaled,
  sigma_vic = vic_layer$sigma,
  vic_pair = vic_layer$pair_index
)


# 5. Flu-all OR constraint input ------------------------------------------
flu_all <- list(
  has_flu_all = ifelse(nrow(flu_all_or) > 0, 1, 0),
  theta_flu_all = flu_all_or$log_value_scaled,
  sigma_flu_all = flu_all_or$sigma,
  flu_A_idx = flu_all_or$flu_A_idx,
  flu_B_idx = flu_all_or$flu_B_idx,
  flu_is_directed = flu_all_or$flu_is_directed
)


# 6. Combine into full Stan input list ------------------------------------
data_list <- list(
  N = N,
  K_or = or_data$K_or,
  theta_or = or_data$theta_or,
  sigma_or = or_data$sigma_or,
  or_pair = or_data$or_pair,
  
  K_log = log_data$K_log,
  theta_log = log_data$theta_log,
  sigma_log = log_data$sigma_log,
  log_pair = log_data$log_pair,
  
  K_vic = vic_data$K_vic,
  theta_vic = vic_data$theta_vic,
  sigma_vic = vic_data$sigma_vic,
  vic_pair = vic_data$vic_pair,
  
  K_flu = length(flu_all$theta_flu_all),
  has_flu_all = flu_all$has_flu_all,
  theta_flu_all = flu_all$theta_flu_all,
  sigma_flu_all = flu_all$sigma_flu_all,
  flu_A_idx = flu_all$flu_A_idx,   # maps to either IAV or IBV in mu_or
  flu_B_idx = flu_all$flu_B_idx,
  flu_is_directed = flu_all$flu_is_directed    # treated as same for now, can extend later
)


# Optional: Save input list for reuse -------------------------------------
save(data_list, file = "data_list.RData")
