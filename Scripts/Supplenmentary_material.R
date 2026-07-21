#Quantitative diagnostics#####
library(posterior)
library(bayesplot)

# Extract posterior draws
#draws <- stan_out$draws()

# Extract sampler diagnostics
#diag <- stan_out$sampler_diagnostics()

# Save both
#saveRDS(draws, "./Output/draws_Mhet.rds")
#saveRDS(diag,  "./Output/diag_Mhet.rds")

# Load draws and diagnostics
#files too large but can be found in the link
#https://drive.google.com/file/d/1PU8dxpQIef8DJ0ioWQ0uH9MFdCkc1JKs/view?usp=sharing
draws <- readRDS("./Output/draws_Mhet.rds")
diag  <- readRDS("./Output/diag_Mhet.rds")
# Rhat and ESS
summ <- summarise_draws(draws)

max_rhat     <- max(summ$rhat,      na.rm = TRUE)
min_ess_bulk <- min(summ$ess_bulk,  na.rm = TRUE)
min_ess_tail <- min(summ$ess_tail,  na.rm = TRUE)

cat("Maximum Rhat:        ", max_rhat, "\n")
cat("Minimum bulk ESS:    ", min_ess_bulk, "\n")
cat("Minimum tail ESS:    ", min_ess_tail, "\n")
# Divergent transitions (3D array)

var_div <- which(dimnames(diag)$variable == "divergent__")
divergent_vals <- diag[, , var_div]
n_divergent <- sum(divergent_vals)
cat("Divergent transitions:", n_divergent, "\n")
# E-BFMI (manual computation)
var_energy <- which(dimnames(diag)$variable == "energy__")
energy_vals <- diag[, , var_energy]

# Flatten all chains
energy_vec <- as.vector(energy_vals)

# BFMI formula: Var(E) / mean(diff(E)^2)
bfmi_val <- var(energy_vec) / mean(diff(energy_vec)^2)

cat("E-BFMI:               ", bfmi_val, "\n")

#end diagnostic

#Compute WAIC and LOO difference for Mhet  and Mhom model#####
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

# Save them
#saveRDS(loo_result, file = "./Output/Mhet_LOO.rds")
#saveRDS(waic_result, file = "./Output/Mhet_WAIC.rds")


#readin models
loo1  <- readRDS("./Output/Mhet_LOO.rds")
waic1 <- readRDS("./Output/Mhet_WAIC.rds")

loo2  <- readRDS("./Output/Mhom_LOO.rds")
waic2 <- readRDS("./Output/Mhom_WAIC.rds")

# LOO comparison
loo_comp <- loo_compare(loo1, loo2)
print(loo_comp)

# Difference ModelA - ModelB
delta_elpd <-
  loo1$estimates["elpd_loo","Estimate"] -
  loo2$estimates["elpd_loo","Estimate"]

# Correct paired SE from loo_compare()
se_delta_elpd <- abs(loo_comp[2,"se_diff"])

cat("LOO difference =", delta_elpd, "\n")
cat("SE of LOO difference =", se_delta_elpd, "\n")


# WAIC comparison

delta_waic <-
  waic1$estimates["waic","Estimate"] -
  waic2$estimates["waic","Estimate"]

# Pointwise WAIC differences
waic_point_diff <-
  waic1$pointwise[, "elpd_waic"] -
  waic2$pointwise[, "elpd_waic"]

# paired SE
se_delta_waic <-
  2 * sqrt(length(waic_point_diff) * var(waic_point_diff))

cat("WAIC difference =", delta_waic, "\n")
cat("SE of WAIC difference =", se_delta_waic, "\n")



# comparison table

comparison_table <- data.frame(
  
  Metric = c("WAIC", "elpd_loo"),
  
  ModelA_Estimate = c(
    waic1$estimates["waic","Estimate"],
    loo1$estimates["elpd_loo","Estimate"]
  ),
  
  ModelA_SE = c(
    waic1$estimates["waic","SE"],
    loo1$estimates["elpd_loo","SE"]
  ),
  
  ModelB_Estimate = c(
    waic2$estimates["waic","Estimate"],
    loo2$estimates["elpd_loo","Estimate"]
  ),
  
  ModelB_SE = c(
    waic2$estimates["waic","SE"],
    loo2$estimates["elpd_loo","SE"]
  ),
  
  Difference = c(
    delta_waic,
    delta_elpd
  ),
  
  SE_Difference = c(
    se_delta_waic,
    se_delta_elpd
  )
  
)

print(comparison_table)
#end

#Compute WAIC and LOO difference for Mhet  and Mhom-het model#####
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

# Save them
#saveRDS(loo_result, file = "./Output/Mhet_LOO.rds")
#saveRDS(waic_result, file = "./Output/Mhet_WAIC.rds")

# Save them
#saveRDS(loo_result, file = "./Output/Mhom-het_LOO.rds")
#saveRDS(waic_result, file = "./Output/Mhom-het_WAIC.rds")

#readin models
loo1  <- readRDS("./Output/Mhet_LOO.rds")
waic1 <- readRDS("./Output/Mhet_WAIC.rds")

loo2  <- readRDS("./Output/Mhom-het_LOO.rds")
waic2 <- readRDS("./Output/Mhom-het_WAIC.rds")

# LOO comparison
loo_comp <- loo_compare(loo1, loo2)
print(loo_comp)

# Difference ModelA - ModelB
delta_elpd <-
  loo1$estimates["elpd_loo","Estimate"] -
  loo2$estimates["elpd_loo","Estimate"]

# Correct paired SE from loo_compare()
se_delta_elpd <- abs(loo_comp[2,"se_diff"])

cat("LOO difference =", delta_elpd, "\n")
cat("SE of LOO difference =", se_delta_elpd, "\n")


# WAIC comparison

delta_waic <-
  waic1$estimates["waic","Estimate"] -
  waic2$estimates["waic","Estimate"]

# Pointwise WAIC differences
waic_point_diff <-
  waic1$pointwise[, "elpd_waic"] -
  waic2$pointwise[, "elpd_waic"]

# paired SE
se_delta_waic <-
  2 * sqrt(length(waic_point_diff) * var(waic_point_diff))

cat("WAIC difference =", delta_waic, "\n")
cat("SE of WAIC difference =", se_delta_waic, "\n")



# comparison table

comparison_table <- data.frame(
  
  Metric = c("WAIC", "elpd_loo"),
  
  ModelA_Estimate = c(
    waic1$estimates["waic","Estimate"],
    loo1$estimates["elpd_loo","Estimate"]
  ),
  
  ModelA_SE = c(
    waic1$estimates["waic","SE"],
    loo1$estimates["elpd_loo","SE"]
  ),
  
  ModelB_Estimate = c(
    waic2$estimates["waic","Estimate"],
    loo2$estimates["elpd_loo","Estimate"]
  ),
  
  ModelB_SE = c(
    waic2$estimates["waic","SE"],
    loo2$estimates["elpd_loo","SE"]
  ),
  
  Difference = c(
    delta_waic,
    delta_elpd
  ),
  
  SE_Difference = c(
    se_delta_waic,
    se_delta_elpd
  )
  
)

print(comparison_table)
#end



#MCMC diagnostics for two states#####

library(bayesplot)
library(ggplot2)
library(posterior)
library(bayesplot)
library(cowplot)

#saving stan code
# Save the full stan_out object
#saveRDS(stan_out, file = "./Output/stan_out_Mhet_McMC.rds")

#read stanout
# LOAD SAVED OBJECTS
#files too large but can be found in the link
#https://drive.google.com/file/d/1xv_df7_xDidF6UrdlO1OYv38sHwh1kZh/view?usp=sharing
#stan_out <- readRDS("./Output/stan_out_Mhet_McMc.rds")
stan_out <- readRDS("./Output/stan_out_Mhet_McMC.rds")

#convert to draws
draws_array <- as_draws_array(stan_out$draws())

# Identify parameters with ANY NA across all iterations & chains
param_names <- dimnames(draws_array)$variable

na_params <- sapply(param_names, function(p) {
  any(is.na(draws_array[, , p]))
})

# Remove NA parameters
clean_draws <- draws_array[, , !na_params]



# Parameters for the two states you want
params_to_check <- c("beta[3]", "v[3]", "I0[3]", "p_reported",  # Edo
                     "beta[10]", "v[10]", "I0[10]", "p_reported") # Rivers


for (param in params_to_check) {
  
  trace_plot <- mcmc_trace(clean_draws, pars = param)
  
  # Clean filename: replace "[" and "]" with "_"
  safe_name <- gsub("\\[|\\]", "_", param)
  
  ggsave(
    filename = paste0("Output/Figure/Traceplot_", safe_name, ".png"),
    plot = trace_plot,
    width = 8, height = 5, dpi = 1000
  )
}

for (param in params_to_check) {
  
  trace_plot <- mcmc_trace(clean_draws, pars = param)
  
  # Clean filename: replace "[" and "]" with "_"
  safe_name <- gsub("\\[|\\]", "_", param)
  
  ggsave(
    filename = paste0("Output/Figure/Traceplot_", safe_name, ".pdf"),
    plot = trace_plot,
    width = 8, height = 5
  )
}


##chain mixxing plots
for (param in params_to_check) {
  
  dens_plot <- mcmc_dens_overlay(clean_draws, pars = param)
  
  safe_name <- gsub("\\[|\\]", "_", param)
  
  ggsave(
    filename = paste0("Density_", safe_name, ".png"),
    plot = dens_plot,
    width = 8, height = 5, dpi = 300
  )
}


##save all density plot
for (param in params_to_check) {
  
  dens_plot <- mcmc_dens_overlay(clean_draws, pars = param)
  
  safe_name <- gsub("\\[|\\]", "_", param)
  
  ggsave(
    filename = paste0("Output/Figure/Density_plot_", safe_name, ".png"),
    plot = dens_plot,
    width = 8,
    height = 5,
    dpi = 1000
  )
}

##save all density plot
for (param in params_to_check) {
  
  dens_plot <- mcmc_dens_overlay(clean_draws, pars = param)
  
  safe_name <- gsub("\\[|\\]", "_", param)
  
  ggsave(
    filename = paste0("Output/Figure/Density_plot_", safe_name, ".pdf"),
    plot = dens_plot,
    width = 8,
    height = 5
  )
}

#Edo pair plot
edo_pairs <- mcmc_pairs(
  clean_draws,
  pars = c("beta[3]", "v[3]", "I0[3]", "p_reported")
)


ggsave(
  filename = "Output/Figure/Figure_S3_notrelabeled.pdf",
  plot = edo_pairs,
  width = 7,
  height = 5,
  bg = "white"
)
ggsave(
  filename = "Output/Figure/Figure_S3_notrelabel.png",
  plot = edo_pairs,
  width = 7,
  height = 5,
  bg = "white",
  dpi=1000
)

#Rivers pair plot
Rivers_pairs <- mcmc_pairs(
  clean_draws,
  pars = c("beta[10]", "v[10]", "I0[10]", "p_reported")
)


#save Rivers
ggsave(
  filename = "Output/Figure/Figure_S4_notrelabel.pdf",
  plot = Rivers_pairs,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "Output/Figure/Figure_S4_notrelabel.png",
  plot = Rivers_pairs,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 1000
)
#

#Relabeling for r
clean_draws_pairs <- clean_draws

vars <- variables(clean_draws_pairs)

vars[vars == "p_reported"] <- "r"

variables(clean_draws_pairs) <- vars

edo_pairs <- mcmc_pairs(
  clean_draws_pairs,
  pars = c("beta[3]", "v[3]", "I0[3]", "r")
)
ggsave(
  filename = "Output/Figure/Figure_S3.pdf",
  plot = edo_pairs,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "Output/Figure/Figure_S3.png",
  plot = edo_pairs,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 1000
)


###
Rivers_pairs <- mcmc_pairs(
  clean_draws_pairs,
  pars = c("beta[10]", "v[10]", "I0[10]", "r")
)
ggsave(
  filename = "Output/Figure/Figure_S4.pdf",
  plot = Rivers_pairs,
  width = 7,
  height = 5,
  bg = "white"
)

ggsave(
  filename = "Output/Figure/Figure_S4.png",
  plot = Rivers_pairs,
  width = 7,
  height = 5,
  bg = "white",
  dpi = 1000
)#end relabeling

# Figure S1: Combined Trace Plots
# Edo
trace_beta3 <- mcmc_trace(clean_draws, pars = "beta[3]") +
  ggtitle(expression(beta~"(Edo)"))

trace_I03 <- mcmc_trace(clean_draws, pars = "I0[3]") +
  ggtitle(expression(I[0]~"(Edo)"))

trace_v3 <- mcmc_trace(clean_draws, pars = "v[3]") +
  ggtitle(expression(nu~"(Edo)"))

edo_trace <- plot_grid(
  trace_beta3,
  trace_I03,
  trace_v3,
  ncol = 1,
  labels = c("A","B","C")
)

# Rivers
trace_beta10 <- mcmc_trace(clean_draws, pars = "beta[10]") +
  ggtitle(expression(beta~"(Rivers)"))

trace_I010 <- mcmc_trace(clean_draws, pars = "I0[10]") +
  ggtitle(expression(I[0]~"(Rivers)"))

trace_v10 <- mcmc_trace(clean_draws, pars = "v[10]") +
  ggtitle(expression(nu~"(Rivers)"))

rivers_trace <- plot_grid(
  trace_beta10,
  trace_I010,
  trace_v10,
  ncol = 1,
  labels = c("D","E","F")
)

# Combine
FigureS1 <- plot_grid(
  edo_trace,
  rivers_trace,
  ncol = 2
)

ggsave(
  "Output/Figure/Figure_S1_TracePlots.pdf",
  FigureS1,
  width = 12,
  height = 10
)

ggsave(
  "Output/Figure/Figure_S1_TracePlots.png",
  FigureS1,
  width = 12,
  height = 10,
  dpi = 1000
)


#Figure S2: Combined Density Plots###


# Edo
dens_beta3 <- mcmc_dens_overlay(clean_draws, pars = "beta[3]") +
  ggtitle(expression(beta~"(Edo)"))

dens_I03 <- mcmc_dens_overlay(clean_draws, pars = "I0[3]") +
  ggtitle(expression(I[0]~"(Edo)"))

dens_v3 <- mcmc_dens_overlay(clean_draws, pars = "v[3]") +
  ggtitle(expression(nu~"(Edo)"))

edo_density <- plot_grid(
  dens_beta3,
  dens_I03,
  dens_v3,
  ncol = 1,
  labels = c("A","B","C")
)

# Rivers
dens_beta10 <- mcmc_dens_overlay(clean_draws, pars = "beta[10]") +
  ggtitle(expression(beta~"(Rivers)"))

dens_I010 <- mcmc_dens_overlay(clean_draws, pars = "I0[10]") +
  ggtitle(expression(I[0]~"(Rivers)"))

dens_v10 <- mcmc_dens_overlay(clean_draws, pars = "v[10]") +
  ggtitle(expression(nu~"(Rivers)"))

rivers_density <- plot_grid(
  dens_beta10,
  dens_I010,
  dens_v10,
  ncol = 1,
  labels = c("D","E","F")
)

# Combine
FigureS2 <- plot_grid(
  edo_density,
  rivers_density,
  ncol = 2
)

ggsave(
  "Output/Figure/Figure_S2_DensityPlots.pdf",
  FigureS2,
  width = 12,
  height = 10
)

ggsave(
  "Output/Figure/Figure_S2_DensityPlots.png",
  FigureS2,
  width = 12,
  height = 10,
  dpi = 1000
)

##end mcmc plots


#posterior predictive check and posterior mean curve######
library(posterior)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)

# Load all required objects
stan_out        <- readRDS("./Output/stan_out_Mhet_McMc.rds")
cases_matrix    <- readRDS("./Output/cases_matrix.rds")
t_last          <- readRDS("./Output/t_last.rds")
selected_states <- readRDS("./Output/selected_states.rds")

n_states <- length(t_last)
max_days <- nrow(cases_matrix)

# 1. Extract posterior draws for pred_cases[s,i]

draws_pred <- as_draws_df(stan_out$draws("pred_cases"))

t_last_df <- tibble(state = 1:n_states, t_last = t_last)

#  Convert pred_cases[s,i] into long format


long_pred <- draws_pred %>%
  pivot_longer(
    cols = starts_with("pred_cases["),
    names_to = "var",
    values_to = "value"
  ) %>%
  mutate(
    var   = str_remove_all(var, "pred_cases\\[|\\]"),
    state = as.integer(str_extract(var, "^[0-9]+")),
    day   = as.integer(str_extract(var, "(?<=,)[0-9]+"))
  ) %>%
  left_join(t_last_df, by = "state") %>%
  filter(day <= t_last)

# Posterior predictive mean + 95% intervals


summary_pred <- long_pred %>%
  group_by(state, day) %>%
  summarise(
    mean   = mean(value),
    median = median(value),
    lower  = quantile(value, 0.025),
    upper  = quantile(value, 0.975),
    .groups = "drop"
  )

#  Observed data

cases_subset <- cases_matrix[1:max(t_last), ]

obs_df <- map_dfr(1:n_states, function(s) {
  tibble(
    state    = s,
    day      = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

#  Join predictions + observations


plot_df <- summary_pred %>%
  left_join(obs_df, by = c("state", "day"))

state_labels <- setNames(selected_states, 1:n_states)

#. Build PPC plot WITH LEGEND


p_ppc <- ggplot(plot_df, aes(x = day)) +
  
  geom_ribbon(
    aes(ymin = lower, ymax = upper, fill = "95% Interval"),
    alpha = 0.35
  ) +
  
  geom_line(
    aes(y = mean, color = "Posterior Mean"),
    size = 1
  ) +
  
  geom_point(
    aes(y = Observed, color = "Observed Data"),
    size = 0.8
  ) +
  
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +
  
  scale_fill_manual(
    name = "Posterior Predictive",
    values = c("95% Interval" = "lightblue")
  ) +
  
  scale_color_manual(
    name = "Posterior Predictive",
    values = c(
      "Posterior Mean" = "blue",
      "Observed Data"  = "black"
    )
  ) +
  
  labs(
    x = "Day",
    y = "Cases",
    title = "Posterior Predictive Intervals and Mean Curve by State"
  ) +
  
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

print(p_ppc)

ggsave(
  filename = "Output/Figure/Figure_S5_plot.pdf",
  plot = p_ppc,
  width = 10,
  height = 8
)

ggsave(
  filename = "Output/Figure/Figure_S5_plot.png",
  plot = p_ppc,
  width = 10,
  height = 8,
  dpi = 1000
)
#end



#prior vs posterior same y-axis######

#saveRDS(stan_out, file = "./Output/stan_out_Mhet_prior.rds")
#saveRDS(data_sir, file = "./Output/data_sir_Mhet_prior.rds")   # optional but recommended


#load in
stan_out <- readRDS("./Output/stan_out_Mhet_prior.rds")
data_sir <- readRDS("./Output/data_sir_Mhet_prior.rds")

#rebuild block
draws_df <- as_draws_df(stan_out$draws())

v3_post  <- draws_df[["v[3]"]]
v10_post <- draws_df[["v[10]"]]

beta3_post  <- draws_df[["beta[3]"]]
beta10_post <- draws_df[["beta[10]"]]

p_reported_post <- draws_df[["p_reported"]]

I0_edo_post    <- draws_df[["I0[3]"]]
I0_rivers_post <- draws_df[["I0[10]"]]

R0_post <- draws_df[["R0"]]


#Prior vesrsus posterior
#Number of prior samples
n_prior <- 5000

# Priors for mixture parameters
theta_ss_prior <- rbeta(n_prior, 1, 1)
slab_sd_prior  <- rgamma(n_prior, 1, 10)

# spike_sd comes from your Stan data list
spike_sd <- data_sir$spike_sd

# Sample v from spike-and-slab prior
v_prior <- numeric(n_prior)
for (i in 1:n_prior) {
  if (runif(1) < theta_ss_prior[i]) {
    v_prior[i] <- rnorm(1, 0, spike_sd)      # spike
  } else {
    v_prior[i] <- rnorm(1, 0, slab_sd_prior[i])  # slab
  }
}

# Plot with rescaled prior density

library(ggplot2)

# Kernel densities
d_prior <- density(v_prior)
d_v3    <- density(v3_post)
d_v10   <- density(v10_post)

# Scale prior down so all curves are visible
scale_factor <- max(c(max(d_v3$y), max(d_v10$y))) /
  max(d_prior$y)

prior_df <- data.frame(
  x = d_prior$x,
  y = d_prior$y * scale_factor,
  type = "Prior (scaled)"
  #type = "Prior"
)

v3_df <- data.frame(
  x = d_v3$x,
  y = d_v3$y,
  type = "Posterior v[3] (Edo)"
)

v10_df <- data.frame(
  x = d_v10$x,
  y = d_v10$y,
  type = "Posterior v[10] (Rivers)"
)

plot_df <- rbind(prior_df, v3_df, v10_df)

p <- ggplot(plot_df,
            aes(x = x,
                y = y,
                colour = type,
                fill = type)) +
  geom_line(size = 1.1) +
  geom_area(alpha = 0.25,
            position = "identity") +
  labs(
    title = "Prior vs Posterior for v",
    x = "v",
    y = "Density"
  ) +
  annotate(
    "text",
    x = min(plot_df$x),
    y = max(plot_df$y),
    hjust = 0,
    #label = paste("Prior density scaled by",
    # round(scale_factor,4))
    label = paste("Prior density rescaled for visual comparison ")
  ) +
  #theme_minimal(base_size = 14)#hide y-axis
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

print(p)

ggsave(
  "Output/Figure/Figure_S7.png",
  p,
  width = 8,
  height = 5,
  dpi = 1000
)

ggsave(
  "Output/Figure/Figure_S7.pdf",
  p,
  width = 8,
  height = 5
)


#for beta
library(posterior)
library(ggplot2)

# Extract posterior draws

beta3_post  <- as_draws_df(stan_out$draws("beta[3]"))$`beta[3]`
beta10_post <- as_draws_df(stan_out$draws("beta[10]"))$`beta[10]`

# Simulate prior samples
# Stan prior: beta ~ normal(0.5, 0.1)

n_prior <- 5000
beta_prior <- rnorm(n_prior, mean = 0.5, sd = 0.1)


# Beta plot with scaled prior

d_prior <- density(beta_prior)
d_b3    <- density(beta3_post)
d_b10   <- density(beta10_post)

scale_factor <- max(c(max(d_b3$y), max(d_b10$y))) /
  max(d_prior$y)

prior_df <- data.frame(
  x = d_prior$x,
  y = d_prior$y * scale_factor,
  type = "Prior (scaled)"
  # type = "Prior"
)

b3_df <- data.frame(
  x = d_b3$x,
  y = d_b3$y,
  type = "Posterior beta[3] (Edo)"
)

b10_df <- data.frame(
  x = d_b10$x,
  y = d_b10$y,
  type = "Posterior beta[10] (Rivers)"
)

plot_df <- rbind(prior_df, b3_df, b10_df)

p_beta <- ggplot(plot_df,
                 aes(x = x,
                     y = y,
                     colour = type,
                     fill = type)) +
  geom_line(size = 1.1) +
  geom_area(alpha = 0.25,
            position = "identity") +
  labs(
    title = "Prior vs Posterior for beta",
    x = expression(beta),
    y = "Density"
  ) +
  annotate(
    "text",
    x = min(plot_df$x),
    y = max(plot_df$y),
    hjust = 0,
    #label = paste("Prior density scaled by",
    #  round(scale_factor,4))
    
    label = paste("Prior density  rescaled for visual comparison ")
  ) +
  #theme_minimal(base_size = 14)#hide y-axis
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

print(p_beta)

ggsave(
  "Output/Figure/Figure_S6.png",
  p_beta,
  width = 8,
  height = 5,
  dpi = 1000
)

ggsave(
  "Output/Figure/Figure_S6.pdf",
  p_beta,
  width = 8,
  height = 5
)
#end prior vs posterior


# POSTERIOR DISTRIBUTION PLOTS#####

library(posterior)
library(ggplot2)

# Save posterior draws
#saveRDS(stan_out$draws(), "Output/draws_model.rds")

# Save Stan data list (needed for priors)
#saveRDS(data_sir, "Output/data_sir.rds")

# Save metadata
#saveRDS(n_states, "Output/n_states.rds")
#saveRDS(t_last, "Output/t_last.rds")
#saveRDS(cases_matrix, "Output/cases_matrix.rds")

# LOAD SAVED OBJECTS
#files too large but can be found in the link
#https://drive.google.com/file/d/1HkCv3cgly33f7g8jrk6h5OhzKhNDA6wP/view?usp=sharing
draws <- readRDS("Output/draws_model.rds")
posterior_df <- as_draws_df(draws)

# Function to create posterior distribution plots
#

make_post_plot <- function(par){
  
  df <- data.frame(value = posterior_df[[par]])
  
  ggplot(df, aes(x = value)) +
    geom_density(
      fill = "steelblue",
      alpha = 0.5,
      colour = "steelblue",
      linewidth = 1
    ) +
    labs(
      title = paste("Posterior Distribution of", par),
      x = par,
      y = "Density"
    ) +
    theme_minimal(base_size = 15)
  
}

# Figure S8a: beta for Edo and Rivers


Figure_S8a <- plot_grid(
  make_post_plot("beta[3]"),
  make_post_plot("beta[10]"),
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

ggsave(
  "Output/Figure/Figure_S8a.pdf",
  Figure_S8a,
  width = 12,
  height = 5,
  bg = "white"
)

ggsave(
  "Output/Figure/Figure_S8a.png",
  Figure_S8a,
  width = 12,
  height = 5,
  dpi = 1000,
  bg = "white"
)

# Figure S8b: mean I0 and mean R0


Figure_S8b <- plot_grid(
  make_post_plot("meanI0"),
  make_post_plot("meanR0"),
  labels = c("A", "B"),
  ncol = 2,
  align = "hv"
)

ggsave(
  "Output/Figure/Figure_S8b.pdf",
  Figure_S8b,
  width = 12,
  height = 5,
  bg = "white"
)

ggsave(
  "Output/Figure/Figure_S8b.png",
  Figure_S8b,
  width = 12,
  height = 5,
  dpi = 1000,
  bg = "white"
)

cat("Combined posterior distribution plots saved in Output/Figure/\n")
