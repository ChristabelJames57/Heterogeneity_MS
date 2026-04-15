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
library(deSolve)
rm(list = ls()) #clears the memory

set.seed(1275) 
total_number_of_time_series=19
MAX_TIME<-189
N=matrix(data=0, nrow = total_number_of_time_series)
tru_incidence_matrix <- obs_incidence_matrix<- matrix(data=0, nrow = MAX_TIME , ncol = total_number_of_time_series)#all in a matrix with that row n column
Incidence=matrix(data=0, nrow = MAX_TIME,ncol=total_number_of_time_series)
obs_incidence_matrix2 <- matrix(0, nrow = MAX_TIME, ncol = total_number_of_time_series)



SIRode <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    beta <- R0 * gamma  
    dS <- -beta * I * (S /(pN* N[i]))
    dI <-  beta * I * (S / (pN*N[i])) - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

getwd()
cmdstanr_max_rows=15
max_rows = 15

The_mode="R"#change only here from  S to R N 
if (The_mode=="R"){
  df_N <- read.csv("./Data/Nigeria17.csv") # Columns: Day, Edo,Lagos, Ondo, Ogun, Osun, Oyo, 
  populations <- c(Ondo = 4671695, Edo = 4235595, Osun = 4705589, Ogun = 5217716, Oyo = 7840864, Lagos = 12550598, Kano = 13076892, Kwara = 3192893, Delta = 5663362, Ebonyi = 2880383, Enugu = 4411119, Kaduna = 8252366, Bayelsa = 2277961, Benue = 5741815, Fct = 3564126, Gombe = 3256962, Rivers = 7303924, Bauchi = 6537314, `Cross River` = 3866269)
  end_days<-c(Ondo = 152, Edo = 189, Osun = 127, Ogun = 197, Oyo = 174, Lagos = 228,Kano =144, Kwara = 154, Delta = 140, Ebonyi = 119, Enugu = 128, Kaduna = 200, Bayelsa =89, Benue =86, Fct= 225, Gombe = 143, Rivers= 157, Bauchi = 129, `Cross River` =89)#real data 15 states
  
 
  # Select States (ten)
  selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara","Ondo", "Osun", "Oyo", "Rivers")#arrge alphabe 13 no Benue
  
  # Prepare Data 
  n_states <- length(selected_states)
  #n_days <- length(df_N$Day)
  t_last <- end_days[selected_states]#end days for selcted states alone
  
  cases_matrix <- df_N %>%#matrix of case data
    select(all_of(selected_states)) %>%#selects the column for specified states
    as.matrix(nrow=t_last,ncol=n_states)#converts to matrix
  
  
  N_vector <- populations[selected_states]
  #end real data
  
} else{#begin simulated data

  beta <- 0.4#
  gamma <- 0.2
  
  R0 <- beta / gamma
  reporting_prob <- 0.0075#using now
  I0 <- rpois(1,1/ reporting_prob)#DIFFERENT Initial conditions

  sim_r <- 2.1
  noof_outbreaks<-5#as you want but must be less than population 12
  n_states <- noof_outbreaks
  selected_states <- paste0("SimState_", 1:noof_outbreaks)
  N_vector<-rep(0, each = noof_outbreaks)
  t_last<-rep(0,each = noof_outbreaks)
  
  
  for (i in 1:noof_outbreaks) {
    N[i] <- 4235595#Edo
    S0 <- N[i] - I0
    R0_init <- 0
    y_init <- c(S = S0, I = I0, R = R0_init)
    
    times <- seq(0, MAX_TIME, by = 1)
    
    
    parameters <- c(beta = beta, gamma = gamma, pN=pN)#remove pN if not estimating it that is for simualtion. Use her if Real
 
    
    out <- ode(y = y_init, times = times, func = SIRode, parms = parameters)
    out_df <- as.data.frame(out)
    
    
    incidence <- -diff(out_df$S)# Calculate incidence as new infections (delta S)
    incidence[incidence < 0] <- 0 
    
    #True process for incidence (deltaS)
    tru_incidence_matrix[,i ] <- incidence
    
    #Observed process
    obs_incidence_matrix[,i ] <- rnbinom(MAX_TIME, mu = reporting_prob * incidence + 1e-6, size = sim_r)
    obs_incidence_matrix2[, i ] <- obs_incidence_matrix[,i ]#so we have day row n outbreak col
    N_vector[i] <- N[i]
    t_last[i]<-MAX_TIME
    
  }
  
  cases_matrix <- obs_incidence_matrix[,1:noof_outbreaks] #matrix of 4 cols
  
}#end simulated


data_sir <- list(max_days = max(t_last),n_states = n_states,t0 = 0,t_last = t_last,N = N_vector,cases = cases_matrix[1:max(t_last), ])

NCovid <-cmdstan_model("./Scripts/Mhom_model.stan")

iters=2000

##Regularises where the initial value begins
init_fun <- function() {
  list(
    pN = rep(0.6, 10),               # vector of same length
    beta = rep(0.4, 10),       # now a vector of length n_states
    p_reported = 0.02,                   # scalar is fine (global parameter)
    I0 = rep(1/0.02, 10)  # vector of same length
  )
}

stan_out <- NCovid$sample(data = data_sir, iter_warmup = iters, iter_sampling = iters, parallel_chains = 4,init = init_fun)

print(stan_out, variables = c("MeanpN","meanbeta", "meanR0", "meanp_reported", "meanI0", "beta", "R0","pN", "I0"), digits = 8,max_rows = 70)



#----Model selction----
log_lik_array <- stan_out$draws("log_lik", format = "matrix")
loo_result <- loo(log_lik_array)
waic_result <- waic(log_lik_array)

print(loo_result)
print(waic_result)

stan_out$diagnostic_summary()#Check divergence



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
state_labels <- setNames(selected_states, c(1, 2, 3, 4, 5,6,7,8,9,10))#replace 1234 with names if 10 states

# Plot using facets
ggplot(plot_df, aes(x = day)) +
  geom_point(aes(y = Observed), color = "black", size = 0.8) +
  geom_line(aes(y = median), color = "blue", size = 1) +
  facet_wrap(~state, labeller = as_labeller(state_labels), scales = "free_y") +#FORlabeling states
  labs(x = "Day", y = "Cases", title = "Observed vs Predicted Cases by State") +
  theme_minimal(base_size = 12)


###TO SAVE MODEL CODE
draws_df_m22 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m22, file = "draws_model22_noheterogeneity.rds")
#end

###TO new SAVE MODEL CODE ESTIMATING PN SAME prior
draws_df_m33 <- as_draws_df(stan_out$draws("pred_cases"))
saveRDS(draws_df_m33, file = "draws_model33_noheterogeneity.rds")
#end

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



