####Model with non-centered parameterization####
####By Suyi####
####20250414####


# 1. load data & package ------------------------------------------------------------
library(rstan)
load("data_list.RData")
#data_list$K_flu <- length(data_list$theta_flu_all)


# 2. model ----------------------------------------------------------------
model_code <- "
data {
  int<lower=1> N;

  int<lower=1> K_or;
  int<lower=1,upper=N> or_pair[K_or];
  vector[K_or] theta_or;
  vector<lower=0>[K_or] sigma_or;

  int<lower=1> K_log;
  int<lower=1,upper=N> log_pair[K_log];
  vector[K_log] theta_log;
  vector<lower=0>[K_log] sigma_log;

  int<lower=1> K_vic;
  int<lower=1,upper=N> vic_pair[K_vic];
  vector[K_vic] theta_vic;
  vector<lower=0>[K_vic] sigma_vic;

  int<lower=1> K_flu;
  vector[K_flu] theta_flu_all;
  vector<lower=0>[K_flu] sigma_flu_all;
  int<lower=1,upper=N> flu_A_idx[K_flu];
  int<lower=1,upper=N> flu_B_idx[K_flu];
}

parameters {
  real w_log;
  real w_vic;
  real intercept;
  real<lower=0> tau;
  vector[N] mu_raw;
}

transformed parameters {
  vector[N] mu_or;

  for (i in 1:N) {
    real mean_log = 0;
    real weight_log = 0;
    real mean_vic = 0;
    real weight_vic = 0;

    for (k in 1:K_log)
      if (log_pair[k] == i) {
        mean_log += theta_log[k] / square(sigma_log[k]);
        weight_log += 1 / square(sigma_log[k]);
      }
    for (k in 1:K_vic)
      if (vic_pair[k] == i) {
        mean_vic += theta_vic[k] / square(sigma_vic[k]);
        weight_vic += 1 / square(sigma_vic[k]);
      }

    real prior_mean = 0;
    if (weight_log > 0 && weight_vic > 0)
      prior_mean = intercept + w_log * (mean_log / weight_log) + w_vic * (mean_vic / weight_vic);
    else if (weight_log > 0)
      prior_mean = intercept + w_log * (mean_log / weight_log);
    else if (weight_vic > 0)
      prior_mean = intercept + w_vic * (mean_vic / weight_vic);
    else
      prior_mean = intercept;

    mu_or[i] = prior_mean + mu_raw[i] * tau;
  }
}

model {
  // Priors
  w_log ~ normal(1, 1);
  w_vic ~ normal(1, 1);
  intercept ~ normal(0, 2);
  tau ~ normal(0, 1);
  mu_raw ~ normal(0, 1);

  // OR likelihood
  for (k in 1:K_or)
    theta_or[k] ~ normal(mu_or[or_pair[k]], sigma_or[k]);

  // Flu average constraint likelihood
  for (k in 1:K_flu)
     theta_flu_all[k] ~ normal(1/2 * mu_or[flu_A_idx[k]] + 1/2 * mu_or[flu_B_idx[k]], sigma_flu_all[k]);
   // IBV占比越小，index越小
    //theta_flu_all[k] ~ normal(2/3 * mu_or[flu_A_idx[k]] + 1/3 * mu_or[flu_B_idx[k]], sigma_flu_all[k]);
}
"
