library(bayesplot)
library(cmdstanr) 
library(tidybayes)
library(tidyverse)
library(gridExtra)
library(outbreaks)
library(readr)
library(ggplot2)
library(dplyr)
library(magrittr)
library(loo)
library(posterior)
library(writexl)
rm(list = ls()) #clears the memory

########FOR REAL AND S  FIXING N AND ESTIMATING v Model 1 ######
                 
set.seed(1275) #real
#set.seed(1275) simu
total_number_of_time_series=19#for simulation 15NL and 3L must bge highr than noof outbvreak
#MAX_TIME<-100 #To test experimets
MAX_TIME<-189 ## To test real data
N=matrix(data=0, nrow = total_number_of_time_series)
tru_incidence_matrix <- obs_incidence_matrix<- matrix(data=0, nrow = MAX_TIME , ncol = total_number_of_time_series)#all in a matrix with that row n column
Incidence=matrix(data=0, nrow = MAX_TIME,ncol=total_number_of_time_series)#so we have day row n outbreak col
obs_incidence_matrix2 <- matrix(0, nrow = MAX_TIME, ncol = total_number_of_time_series)


SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma  # recompute beta from R0 and gamma
    dS <- -beta * I * (S / N[i])^(1 + v^2)
    dI <-  beta * I * (S / N[i])^(1 + v^2) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

#----model2---
#SIRode <- function(time, state, parameters) {
# with(as.list(c(state, parameters)), {
#beta <- R0 * gamma  # recompute beta from R0 and gamma
#dS <- -beta * I * (S / (pN*N[i]))^(1 + 0^2)
#dI <-  beta * I * (S / (pN*N[i]))^(1 + 0^2) - gamma * I
# dR <-  gamma * I
#  return(list(c(dS, dI, dR)))
# })
#}


#Code that sets rp seperately for the different states, has a decay process on rp
getwd()
cmdstanr_max_rows=15
max_rows = 15
getwd()
setwd("C:/Users/2847125J/OneDrive - University of Glasgow/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025/Attempt2/Attempt2")
#setwd("C:/Users/Christabel/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025/Attempt2/Attempt2")
#setwd("~/Desktop/Dans Folder 2023/People/Christabel James/Stan code/Attempt3")

The_mode="R"#change only here from  S to R N chge prior for rp
if (The_mode=="R"){
  df_N <- read.csv("Nigeria17.csv") # Columns: Day, Edo,Lagos, Ondo, Ogun, Osun, Oyo, 
  populations <- c(Ondo = 4671695, Edo = 4235595, Osun = 4705589, Ogun = 5217716, Oyo = 7840864, Lagos = 12550598, Kano = 13076892, Kwara = 3192893, Delta = 5663362, Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366, Bayelsa = 2277961, Benue = 5741815, Fct = 3564126, Gombe = 3256962, Rivers = 7303924, Bauchi = 6537314, `Cross River` = 3866269)
  #end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200)#real data 15 states
  end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200, Bayelsa =89, Benue =86, Fct= 225, Gombe = 143, Rivers= 157, Bauchi = 129, `Cross River` =89)#real data 15 states
  
  # Select States
  #selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kano", "Kwara", "Ondo", "Osun", "Oyo")  # alphabetical order
  #selected_states <- c("Bayelsa", "Benue", "Delta", "Ebonyi", "Edo", "Enugu", 
  #"Gombe", "Kaduna", "Kano", "Kwara", "Ondo", "Osun", "Oyo", "Rivers")#arrange alphabetical order
  
  # Select States
  #selected_states <- c("Bayelsa", "Delta", "Ebonyi", "Edo", "Enugu", 
  # "Gombe", "Kaduna", "Kano", "Kwara", "Ondo", "Osun", "Oyo", "Rivers")#arrge alphabe 13 no Benue
  
  # Select States (ten)
  selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara","Ondo", "Osun", "Oyo", "Rivers")#arrge alphabe 13 no Benue
  
  # Prepare Data 
  n_states <- length(selected_states)
  #n_days <- length(df_N$Day)
  t_last <- end_days[selected_states]#end days for selcted states alone
  
  cases_matrix <- df_N %>%#matrix of case data
    select(all_of(selected_states)) %>%#selects the column for specified state
    as.matrix(nrow=t_last,ncol=n_states)#converts to matrix
  
  
  N_vector <- populations[selected_states]
  #end real data
  
} else{#begin simulated data
  
  # beta <- 0.2731#mean real data 7 state
  # beta <- 0.2679#mean real data 10  state
  beta <- 0.4#
  gamma <- 0.2
  #pN=0.4
  # v <- 2.6619#mean real data 7 statw
  #v <-2.8694# real 10
  v <-2 #Fixing homogenous
  R0 <- beta / gamma
  #I0 <- 1124.576#meanI0 7 sta
  #I0 <- 1169.240#meanI0 10 sta
  #reporting_prob <-  0.0185#mean pr testing good here as rp is large
  #reporting_prob <-  0.00185#mean pr testing under est v here as rp is small
  reporting_prob <-  0.015#mean pr testing 1/pr-
  I0 <- rpois(1,1/ reporting_prob)#DIFFERENT Initial conditions
  #I0 <- 540.5405#1/pr same initial conditions
  #reporting_prob <- 0.015#good esti but not monotone
  #reporting_prob <- 0.0066#mean pr real 7
  #reporting_prob <- 0.0071#mean pr real 10
  sim_r <- 2.1
  noof_outbreaks <- 10 # using 10 simulated outbreaks as u
  n_states <- noof_outbreaks
  selected_states <- paste0("SimState_", 1:noof_outbreaks)
  
  # Delta, Ebonyi, Edo, Enugu, Kaduna, Kwara, Ondo, Osun, Oyo, Rivers
  pop_order <- c(
    5663362, # Delta
    2880383, # Ebonyi
    4235595, # Edo
    4411119, # Enugu
    8252366, # Kaduna
    3192893, # Kwara
    4671695, # Ondo
    4705589, # Osun
    7840864, # Oyo
    7303924  # Rivers
  )
  
  # pre-allocate vectors / matrices
  N_vector <- integer(noof_outbreaks)        # hold integer populations passed to Stan
  t_last <- integer(noof_outbreaks)
  
  # Make sure the global N that you used earlier is numeric/integer and sized at least noof_outbreaks:
  # (old code had N = matrix(...). We'll keep using N[] as before but ensure assignment)
  # If N was a matrix, you can keep it; here we coerce to numeric vector for S0 computation:
  if (!is.numeric(N)) N <- as.numeric(N)
  
  if (length(N) < noof_outbreaks) N <- numeric(noof_outbreaks)  # ensure N has at least noof_outbreaks values
  
  for (i in 1:noof_outbreaks) {
    N[i] <- pop_order[i]
    S0 <- N[i] - I0
    R0_init <- 0
    y_init <- c(S = S0, I = I0, R = R0_init)
    
    times <- seq(0, MAX_TIME, by = 1)
    
    parameters <- c(beta = beta, gamma = gamma, v = v)
    
    out <- ode(y = y_init, times = times, func = SIRode, parms = parameters)
    out_df <- as.data.frame(out)
    
    incidence <- -diff(out_df$S) # new infections
    incidence[incidence < 0] <- 0
    
    # True process
    tru_incidence_matrix[, i] <- incidence
    
    # Observed process (Neg-Bin)
    #obs_incidence_matrix[, i] <- rnbinom(MAX_TIME, mu = reporting_prob * incidence + 1e-6, size = sim_r)# nBINOM
    obs_incidence_matrix[, i] <- rbinom(MAX_TIME, size = (trunc(incidence)+1), prob = reporting_prob)#binom
    obs_incidence_matrix2[, i] <- obs_incidence_matrix[, i]
    
    N_vector[i] <- as.integer(N[i])    # store population and t_last
    t_last[i] <- MAX_TIME
  }
  
  cases_matrix <- obs_incidence_matrix[, 1:noof_outbreaks] # matrix of cases
  
}#end simulated

#I0 <- 1
#y0_array <- t(sapply(N_vector, function(N) c(N - I0, I0, 0)))  # [n_states, 3]#here R=0, I=1 and S=N-1


#N = as.integer(N_vector)# ensure integers#SImulation added converts N vectors to integers after tlast

data_sir <- list(max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,N = N_vector,cases = cases_matrix[1:max(t_last), ], spike_sd = 1e-3)#pas as input stan no phi and REAL estim IO
#data_sir <- list(phi=5,max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,y0 = y0_array,N = N_vector,cases = cases_matrix[1:max(t_last), ], spike_sd = 1e-3)#pas as input to stan
#data_sir <- list(phi=5,max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,y0 = y0_array,N = N_vector,cases = cases_matrix[1:max(t_last), ])#pas as input to stan

#NCovid <- cmdstan_model("NCovid_slabnspike_v_mean3.stan")# compile model suscep no random effc
#NCovid <- cmdstan_model("NCovid_Multi_Normal_FixvPr_Real_Estimate Pr.stan")#model 1 normal estimate pr
#NCovid <-cmdstan_model("NCovid_Multi_Normal_FixvPr_Real_Estim Pr_RemoveIO.stan")#norm esti RemoveIO
#NCovid <-cmdstan_model("NCovid_Multi_Normal_dummy_Real_Estimate Pr.stan")#estimat dummy
#NCovid <- cmdstan_model("NCovid_Multi_Normal_FixvPr_Realuse.stan")#model 1 normal fix pr
#NCovid <- cmdstan_model("NCovid_slabnspike_v_mean3_Fix pr.stan")#model 1
#NCovid <- cmdstan_model("fixing_params_for_spike_mean3.stan")
#NCovid <- cmdstan_model("Stan_code_ with_final_epidemic_size_Hetero1.stan")#epidemic size
#NCovid <- cmdstan_model("NCovid_Multi_Normal_FixvPr_Real_Estim Pr_RemoveIO_Global.stan")#state speicif _GOODJANuary2026

NCovid <- cmdstan_model("NCovid_Multi_Normal_FixvPr_Real_Estim Pr_RemoveIO_Global1.stan")#state speicif _GOODFEBRUARY2026#NCovid <- cmdstan_model("NCovid_Multi_nb_rpt_v_w.stan")# compile model suscep no random effc
#NCovid<- cmdstan_model("NCovid_Multi_nb_rpt_v_Rdeff_NEW.stan")#RDF suscep _NEW
#NCovid<- cmdstan_model("NCovid_Multi_nb_rpt_v_w_conec.stan")#conncetivity no randaom eff

iters=2000 #4000 for simu
stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4, seed = 0)
#stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4, seed = 0, adapt_delta = 0.99, max_treedepth = 15)# helps for more convergen


#print(stan_out, variables = c("MeanV", "meanbeta", "meanp_reported","meanR0", "meanI0", "theta_ss","v", "beta","R0", "I0"), digits = 8,max_rows = 70)
print(stan_out, variables = c("MeanV", "meanbeta", "meanp_reported","meanR0", "meanI0", "theta_ss", "I0", "v", "beta"), digits = 8,max_rows = 70)
#summary_df <- as.data.frame(stan_out$summary())#summary to see if much
#write_xlsx(summary_df, "stan_summary_001.xlsx")
#write_xlsx(summary_df, path = "posterior_summary.xlsx")

####saving results out
#all apram
summary_df <- as.data.frame(stan_out$summary())

#summary for MeanV, meanbeta, meanR0 
summary_selected <- summary_df %>%
  filter(variable %in% c("MeanV", "meanbeta", "meanR0"))

#all full draws
draws <- as_draws_df(stan_out)

# only MeanV, meanbeta, meanR0 columns
draws_selected <- draws %>% select(MeanV, meanbeta, meanR0)

#multi sheet all
write_xlsx(
  list(
    full_summary = summary_df,
    selected_summary = summary_selected,
    selected_draws = draws_selected
  ),
  "stan_results_alldummyrea_removeI0_Nor_10goodst.xlsx"
)#END

###------save not to reload
# Save full object (so you can reload later without rerunning MCMC)
saveRDS(stan_out, file = "stan_outremoveI0_Nor_10godst.rds")

# Reload later
stan_out <- readRDS("stan_outremoveI0_Nor_10godst.rds")#end

# Convert draws to a dataframe
posterior_df <- as_draws_df(stan_out)

# Save all draws
write.csv(posterior_df, "posterior_samplesremoveI0_Nor_10godst.csv", row.names = FALSE)

#save summaries only
summ <- stan_out$summary()
write.csv(summ, "posterior_summaryremoveI0_Nor_10godst.csv", row.names = FALSE)#end
##---end saving


##.... to check for R0 CI##
draws <- stan_out$draws()
R0_samples <- as_draws_df(draws)[["R0"]]
quantile(R0_samples, probs = c(0.025, 0.5, 0.975))#2.5% and 7.5% and mean 50%

#----Model selction----
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

print(loo_result)
print(waic_result)

stan_out$diagnostic_summary()#Check divergence

#posterior_samples <- as.data.frame(stan_out$draws(format = "draws_df"))
posterior_samples <- stan_out$draws(format = "draws_df")#trace plots
#######################. code for traces ######################
all_vars <- variables(stan_out$draws())

# Exclude generated quantities with many entries or possible NAs
exclude_patterns <- c("pred_cases", "log_lik", "y", "incidence")
valid_vars <- all_vars[!sapply(all_vars, function(v) any(sapply(exclude_patterns, grepl, v)))]

# Split into chunks of size 10
chunk_size <- 2
var_chunks <- split(valid_vars, ceiling(seq_along(valid_vars) / chunk_size))

# Plot each chunk
for (i in seq_along(var_chunks)) {
  message(paste("Plotting trace for chunk", i, "/", length(var_chunks)))
  
  current_chunk <- var_chunks[[i]]
  
  # Try plotting and catch any errors due to NAs or missing variables
  tryCatch({
    p <- mcmc_trace(stan_out$draws(current_chunk), 
                    facet_args = list(ncol = 2)) +
      ggtitle(paste("Trace Plots: Chunk", i))
    print(p)
  }, error = function(e) {
    message("Error in chunk ", i, ": ", conditionMessage(e))
  })
}


#  Posterior vs Prior (not fully working)
pars <- c("beta", "v", "p_reported")
#pars <- c("beta", "v", "phi_inv", "p_reported")
n <- 10000
prior <- tibble(
  beta = abs(rnorm(n * n_states, 0.5, 0.1)),
  #v = rgamma(n * n_states, 1, 1),#old
  v = rgamma(n * n_states, 1, 10),#spike new used
  #phi_inv = rexp(n * n_states, 50),
  p_reported = rbeta(n * n_states, 1, 10000)#Real
  # p_reported = rbeta(n * n_states, 1, 1000)#simulated old
  #p_reported = rbeta(n * n_states, 0.001, 1)#simulated current try
)



prior_df <- prior %>%
  mutate(state = rep(selected_states, each = n)) %>%
  pivot_longer(cols = -state, names_to = "parameter", values_to = "value") %>%
  mutate(type = "Prior")

posterior_df <- posterior_samples %>%
  select(matches(paste0("^(", paste(pars, collapse = "|"), ")\\["))) %>%
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  separate(param, into = c("parameter", "index"), sep = "\\[", remove = TRUE) %>%
  mutate(index = str_remove(index, "\\]")) %>%
  mutate(state = selected_states[as.numeric(index)], type = "Posterior")

combined_df <- bind_rows(prior_df, posterior_df) %>%
  mutate(parameter = factor(parameter, levels = pars))

##pLot of posterior vs prior
ggplot(combined_df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_grid(parameter ~ state, scales = "free") +
  #facet_wrap(~parameter, scales = "free_x")+
  scale_fill_manual(values = c("Prior" = "gray70", "Posterior" = "blue")) +
  theme_minimal() +
  labs(title = "Prior vs Posterior Distributions", x = "Value", y = "Density", fill = NULL)



#######################. code for predictions ######################
# Extract posterior draws
draws_df <- as_draws_df(stan_out$draws("pred_cases"))

# Dimensions
n_states <- length(t_last)
max_days <- dim(cases_matrix)[1]

# Name t_last for joining
t_last_df <- tibble(state = 1:n_states, t_last = t_last)

# Pivot draws to long format and parse indices
long_draws <- draws_df %>%
  pivot_longer(cols = starts_with("pred_cases["), names_to = "var", values_to = "value") %>%
  mutate(
    var = str_remove_all(var, "pred_cases\\[|\\]"),
    state = as.integer(str_extract(var, "^[0-9]+")),
    day = as.integer(str_extract(var, "(?<=,)[0-9]+"))
  ) %>%
  left_join(t_last_df, by = "state") %>%
  filter(day <= t_last[state])  # Keep only valid days per state

# Summarise draws: median and 95% credible interval
summary_df <- long_draws %>%
  group_by(state, day) %>%
  summarise(
    median = median(value),
    #lower = quantile(value, 0.25),   # 25% → for 50% CI
    #upper = quantile(value, 0.75),   # 75% → for 50% CI
    lower = quantile(value, 0.025),
    upper = quantile(value, 0.975),
    .groups = "drop"
  )

# Observed data prep
cases_subset <- as.matrix(cases_matrix[1:max(t_last), ])

obs_df <- purrr::map_dfr(1:n_states, function(s) {
  tibble(
    state = s,
    day = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

# Join predictions and observations
plot_df <- summary_df %>%
  left_join(obs_df, by = c("state", "day"))
state_labels <- setNames(selected_states, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))#replac 1234 with abc for four states
#state_labels <- setNames(selected_states, c(1, 2, 3, 4, 5,6,7,8,9,10,11,12,13))#replace 1234 with names if 10 states

# Plot using facets
ggplot(plot_df, aes(x = day)) +
  geom_point(aes(y = Observed), color = "black", size = 0.8) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.15)+#to lighten the colour for bette visbility
  #geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) +##CI
  #facet_wrap(~state, scales = "free_y") +
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +#FORlabeling states
  labs(x = "Day", y = "Cases", title = "Observed vs Predicted Cases by State") +
  theme_minimal(base_size = 12)#end
 

###TO SAVE MODEL 1

library(posterior)

draws_df_m1 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m1, file = "draws_model1_heterogeneity.rds")

###new used
library(posterior)

draws_df_m11 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m11, file = "draws_model11_heterogeneity.rds")

###save t_last
t_last <- c(152, 189, 127, 197, 174, 228, 144, 154, 140, 119)

saveRDS(t_last, "t_last.rds")

##save case matrix only once
library(tidyverse)

df_N <- read.csv("Nigeria17.csv")

selected_states <- c(
  "Delta", "Ebonyi", "Edo", "Enugu", "Kaduna",
  "Kwara", "Ondo", "Osun", "Oyo", "Rivers"
)

cases_matrix <- df_N %>%
  select(all_of(selected_states)) %>%
  as.matrix()

saveRDS(cases_matrix, "cases_matrix.rds")
saveRDS(selected_states, "selected_states.rds")


###TO SAVE MODEL Conncetivity

library(posterior)

draws_df_m3 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m3, file = "draws_model1_heterogeneity.rds")



##### FINAL EPIDEMIC SIZE
final_size_draws <- stan_out$draws("final_size")
library(posterior)
final_size_df <- as_draws_df(final_size_draws)
apply(final_size_df, 2, mean)
apply(final_size_df, 2, quantile, probs = c(0.025, 0.5, 0.975))

##mean
state_names <- selected_states
mean_final_size <- setNames(
  apply(final_size_df, 2, mean),
  state_names
)

mean_final_size

####mean posterior
overall_final_size <- rowMeans(final_size_df)

mean(overall_final_size)#315344
quantile(overall_final_size, probs = c(0.025, 0.5, 0.975))

#
#Prior shape before posterio
prior %>%
  pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_density(fill = "gray70", alpha = 0.5) +
  facet_wrap(~parameter, scales = "free") +
  theme_minimal() +
  labs(title = "Prior Distributions Only", x = "Value", y = "Density")


###### For Global beta, v and P_repoeted
# Extract posterior draws as a data frame
posterior_samples <- as_draws_df(stan_out$draws())

# Parameters that actually exist in your model
pars <- c("beta", "v", "p_reported", "I0")

# ---- Build posterior df ----
posterior_df <- posterior_samples %>%
  select(any_of(pars) | starts_with("I0[")) %>%  # Include vector parameter I0
  pivot_longer(cols = everything(), names_to = "param", values_to = "value") %>%
  mutate(
    # Separate vector index (for I0) if present
    parameter = str_remove(param, "\\[.*\\]"),
    index = if_else(str_detect(param, "\\["),
                    as.numeric(str_extract(param, "(?<=\\[)[0-9]+(?=\\])")),
                    NA_real_),
    state = if_else(!is.na(index), selected_states[index], "Global"),
    type = "Posterior"
  )

# ---- Build prior df ----
n <- 10000
prior_list <- list(
  beta = tibble(value = abs(rnorm(n, 0.5, 0.1)), parameter = "beta", state = "Global"),
  v = tibble(value = rgamma(n, 1, 1), parameter = "v", state = "Global"),
  p_reported = tibble(value = rbeta(n, 1, 1000), parameter = "p_reported", state = "Global"),
  I0 = tibble(value = rnorm(n * length(selected_states), 1001, 500),
              parameter = "I0",
              state = rep(selected_states, each = n))
)

prior_df <- bind_rows(prior_list) %>% mutate(type = "Prior")

# ---- Combine ----
combined_df <- bind_rows(prior_df, posterior_df) %>%
  mutate(parameter = factor(parameter, levels = pars))

# ---- Plot ----
ggplot(combined_df, aes(x = value, fill = type)) +
  geom_density(alpha = 0.5) +
  facet_grid(parameter ~ state, scales = "free") +
  scale_fill_manual(values = c("Prior" = "gray70", "Posterior" = "blue")) +
  theme_minimal() +
  labs(title = "Prior vs Posterior Distributions", x = "Value", y = "Density", fill = NULL)





#mean(c(1.572460 ,1.723179 ,1.898559 ,1.635993 ,1.697177 ,1.636589 ,1.622174 ,1.546629 ,1.666607 ,1.743275 ,1.667721 ,1.876117 ,1.645747 ,1.768622)) 
#logit(0.966142)
exp( -0.878637)#0.4153 beta
exp(0.966142)#v02.627787
inv_logit(-3.765857)
plogis(p0) 
plogis(-3.765857) #0.022624
1/(1+exp(-p0))
1/(1+exp(-(-3.765857)))


#CI for beta0
exp(-0.884670)#lower=0.4128504
exp(-0.872294)#upper 0.4179916
#Ci for vo
exp(0.883079)# lower  2.418334
exp(1.049764)#upper  2.856977
#CI for p0
plogis(-3.885124)#lower 0.02013167
plogis(-3.640278)#upper 0.02557386


x <- seq(0, 0.1, length.out = 200)
plot(x, dbeta(x, 0.02, 4), type="l")

# Parameters
shape <- 1001
rate <- 1.0

# Create a sequence of x values around the mean (which is shape/rate)
x <- seq(900, 1100, by = 0.1)

# Compute the density
y <- dgamma(x, shape = shape, rate = rate)

# Plot it
plot(x, y, type = "l", main = "Gamma(1001, 1.0) Distribution",
     xlab = "x", ylab = "Density", col = "blue", lwd = 2)


#If fixing paramters
# FIXED VALUES
true_beta <- rep(0.4, n_states)
true_p_reported <- rep(0.001, n_states)
true_phi_inv <- 5


data_sir <- list(phi=5, phi_inv = true_phi_inv, max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,y0 = y0_array,N = N_vector,cases = cases_matrix[1:max(t_last), ], spike_sd = 1e-3, beta = true_beta,
                 p_reported = true_p_reported)#pas as input to stan



##correlation btwn v and I0
#to save just once
library(posterior)

#posterior_draws <- as_draws_df(
#stan_out$draws(c("v", "I0"))## helps rename stan_out
#)

#saveRDS(
# posterior_draws,
# file = "posterior_v_I0_states.rds"
#)

###start here
rm(list = ls())
gc()

library(dplyr)
library(tidyr)
library(ggplot2)
#loading up data
draws <- readRDS("posterior_v_I0_states.rds")

draws_long <- draws %>%
  pivot_longer(
    cols = matches("^(v|I0)\\["),
    names_to = c("param", "state"),
    names_pattern = "(v|I0)\\[(\\d+)\\]",
    values_to = "value"
  ) %>%
  pivot_wider(
    names_from = param,
    values_from = value
  ) %>%
  mutate(state = as.integer(state))

###corre by state
state_correlations <- draws_long %>%
  group_by(state) %>%
  summarise(cor_v_I0 = cor(v, I0))

print(state_correlations)


#plot
ggplot(draws_long, aes(x = I0, y = v)) +
  geom_point(alpha = 0.15, size = 0.4) +   # smaller points old was 0.8
  #geom_smooth(method = "lm", color = "red", se = FALSE) + #remove red lines
  facet_wrap(~ state, scales = "free") +
  theme_bw() +
  theme(
    # axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm")
  ) +
  labs(
    #title = "Posterior relationship between I0 and v across states",
    x = expression(I[0]),
    y = "v"
  )

#end
 