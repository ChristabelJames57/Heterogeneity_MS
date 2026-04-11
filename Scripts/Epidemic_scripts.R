##Figure 1: Incidence data for the first wave 
####manuscript edit  After adding MHET-hom model and third talk 11/04/2026
library(tidyverse)
library(plotly)
library(scales)
library(ggforce)
library(zoo)  
rm(list = ls()) # Clear workspace

getwd()
#setwd("C:/Users/Christabel/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025")
setwd("C:/Users/2847125J/OneDrive - University of Glasgow/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025")
Data<-read_csv("covid19_data.csv")
Data
head(Data)
Data1 <- Data %>% #Ton use the new data
  rename(Date_of_reg=`Date of re`,Age_u=`Age unit`) %>% #To rename the varaiables, new name=old name
  separate(Date_of_reg, into = c("month","day","year"), "/") %>% # Does in layers
  mutate(across(c(month,day,year),as.numeric)) %>% #create a data column for day month and name numeric 
  filter(year>=2020)# 

Data_wave1 <- Data1 %>%
  mutate(Date = make_date(year, month, day)) 

##To filter for firSt wave with selected states
selected_states <- c("Bayelsa", "Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Oyo", "Rivers")  # Replace with the states you're interested in
Data_wave1 <- Data_wave1 %>%
  filter(Date >= as.Date("2020-02-27") & Date <= as.Date("2020-10-31"),
         `State of Residence` %in% selected_states )# means to filter both dates and states
Data_wave1
summary(Data_wave1$Date)#shows the min and max too
Data_wave1 <-Data_wave1 %>% 
  unite("date",c(year, month, day),sep = "-",remove = FALSE )
Data_wave1$date<-ymd(Data_wave1$date)#to convert to date format in R
str(Data_wave1)
colnames(Data_wave1)
names(Data_wave1)
unique(Data_wave1$year)# shows all the unique values in year
head(Data_wave1)#
sum(is.na(Data_wave1$Date))#check for missing dates
Data_wave1 <- Data_wave1 %>% filter(!is.na(Date))#To filter missing date


#write.csv(Data_wave1, "../Heterogeneity_MS/wave1_data.csv")# to sedn to manuscript place

##Figure 1 starts

getwd()
setwd("C:/Users/2847125J/OneDrive - University of Glasgow/Downloads/BEll/GLASGOW 2023/REsearch/Heterogeneity_MS/Data")
Data_wave1<-read.csv("wave1_data.csv")# to now read in for manuscript
###Data has day, month and year can create a proper date column
Data_wave1 <- Data_wave1 %>%
  mutate(Date = make_date(year, month, day)) 
names(Data_wave1)

##Plot of number of cases over time
# Aggregate data by date
D_State <- Data_wave1 %>%
  group_by(Date,`State of Residence`, `Case classification`, `Vaccination status`, `Outcome of case`, `Quarantine`) %>% #can do multiple group by of variables
  summarize(Case_Count = n()) %>% 
  na.omit(`State of Residence`) #To remove the NA in state

#To check if any NA
any(is.na(D_State))# in the entire data frame of D_state. if NA will return true if no will return False
any(is.na(D_State$date))# to check for NA in date column
any(is.na(D_State$Case_Count))#check for NA in the case_count column
any(is.na(D_State$`Case classification`))#check for NA in the case classification col
any(is.na(D_State$`State of Residence`))#check for NA in states

#Aggregate Cases by date and state
cases_by_state <- Data_wave1 %>% ##Shows the data frame to use for plotting
  group_by(State = `State of Residence`, Date) %>%
  summarise(Case_Count = n(), .groups = "drop")


#Title: Time series of cases by state
##On a single plot of all 10
library(ggplot2)
library(dplyr)
###clean
cases_by_state_clean <- na.omit(cases_by_state)#remove NA before ploting


# Define a vector with the specific states to plot
selected_states <- c("Bayelsa", "Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Oyo", "Rivers")  # Replace with the states you're interested in

# Filter the data to include only the selected states
cases_by_state_filtered <- cases_by_state_clean %>%
  filter(State %in% selected_states)

plot <- ggplot(cases_by_state_filtered, aes(x = Date, y = Case_Count)) +
  geom_line(color = "blue", linewidth = 0.4) +
  
  # Facet layout similar to the image
  facet_wrap(~ State, scales = "free_y", ncol = 4) +
  
  labs(
    #title = "Time Series of Cases by State",
    x = "Date",
    y = "Number of Cases"
  ) +
  
  theme_minimal(base_size = 11) +
  
  theme(
    # Title styling
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    ),
    
    # Facet labels
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    
    # Axis text formatting
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_text(size = 8),
    
    # Clean panel look
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines")
  )

print(plot)
#END Figure 1




### Fig. 2 MAP OF SELECTED STATES 

#########map for all 10 states 
library(sf)
library(ggplot2)
library(dplyr)

# Folder to save the shapefile
output_dir <- "C:/Users/2847125J/OneDrive - University of Glasgow/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025"

# Create a full path for the zip file
zip_path <- file.path(output_dir, "gadm41_NGA_shp.zip")


#  Download the Nigeria shapefile (Level 1)
download.file("https://geodata.ucdavis.edu/gadm/gadm4.1/shp/gadm41_NGA_shp.zip", destfile = zip_path)

# Unzip it into the same directory
unzip(zip_path, exdir = output_dir)

# Load the shapefile using sf

shapefile_path <- file.path(output_dir, "gadm41_NGA_1.shp")
nigeria_states <- st_read(shapefile_path)


#prepare data
# Check the column names
print(names(nigeria_states))

# Clean and format the state names
nigeria_states$name <- tools::toTitleCase(tolower(nigeria_states$NAME_1))

selected_states <- c("Ondo", "Edo", "Oyo", "Delta", "Enugu", "Kaduna", "Ebonyi",  
                     "Kwara", "Rivers", "Bayelsa")#10 states

# Mark states as "Selected" or "Other"
nigeria_states <- nigeria_states %>%
  mutate(selected = ifelse(name %in% selected_states, "Selected", "Other"))

# Plot the map

ggplot(nigeria_states) +
  geom_sf(aes(fill = selected), color = "black", size = 0.2) +
  scale_fill_manual(values = c("Selected" = "skyblue", "Other" = "grey90")) +
  theme_minimal() +
  #labs(title = "Map of Nigeria Highlighting 17 Selected States",
  #fill = "State Category")
  labs(title = "",
       fill = "State Category")

ggsave("nigeria_selected_states_map.png", width = 10, height = 8)# look for this title

###END Fig.2 




###Fig. 3---Recovery of R0 parameter from simulation----### 
library(ggplot2)
library(dplyr)
library(readr)

# Load data
df_R0 <- read_csv("MeanR0_summaryNEW.csv")

# Treat v as factor
df_R0 <- df_R0 %>%
  mutate(v = factor(v))

# Position dodge for consistent separation
pd <- position_dodge(width = 0.0008)

ggplot(df_R0, aes(x = p_r, y = MeanR0, colour = v, group = v)) +
  
  # Mean estimate
  geom_point(
    size = 3,
    position = pd
  ) +
  
  # 95% credible intervals
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0004,
    linewidth = 0.7,
    position = pd
  ) +
  
  # True R0 reference line
  geom_hline(
    yintercept = 2,
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.8
  ) +
  
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red"
    )
  ) +
  
  labs(
    x = expression(p[reported]),
    y = expression(hat(R)[0])
  ) +
  
  scale_y_continuous(
    limits = c(1.95, 2.05),
    breaks = c(1.95, 2.00, 2.05)
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Keep clean look
    panel.grid = element_blank(), # remove all grid lines
    
    plot.margin = margin(10, 10, 10, 10)
  )

# Save
ggsave(
  "R0_recovery_by_v_clean.pdf",
  width = 7,
  height = 5,
  dpi = 300
)
#end Fig. 3


### Fig. 4 (Mean p_reported) along p_reported — edited for manuscript
library(ggplot2)
library(dplyr)
library(readr)

# Read data
df <- read_csv("MeanPreported_summaryNEW.csv")

# Treat v as factor
df <- df %>%
  mutate(v = factor(v))

# Plot
ggplot(df, aes(x = p_r, y = Meanp_r, colour = v)) +
  
  # Mean estimates
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  
  # 95% credible interval (horizontal dashes)
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = LowerCI, yend = LowerCI
    ),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(
      x = p_r - 0.00015, xend = p_r + 0.00015,
      y = UpperCI, yend = UpperCI
    ),
    linewidth = 0.9
  ) +
  
  # True r reference lines
  geom_hline(
    aes(yintercept = p_r_true),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Facet by ν (labels removed below)
  facet_wrap(~ v, ncol = 1) +
  
  # Manual colours
  scale_colour_manual(
    values = c(
      "0" = "blue",
      "1" = "green",
      "2" = "red"
    )
  ) +
  
  labs(
    x = expression(p[r]),
    y = expression(hat(p[r]))
  ) +
  
  # ---- UNIFORM THEME ----
theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    panel.grid = element_blank(),
    
    # REMOVE facet labels completely
    strip.text = element_blank(),
    strip.background = element_blank(),
    
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )

##end Figur4 


### Figure. 5 (Mean v) along p_reported — edited for manuscript 
library(ggplot2)
library(dplyr)
library(readr)

# Read data
df <- read_csv("MeanV_summary3.csv")

# Treat v_true as a factor
df <- df %>%
  mutate(v_true = as.factor(v_true))

# Plot
ggplot(df, aes(x = p_r, y = MeanV, color = v_true, group = v_true)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0005,
    linewidth = 0.8
  ) +
  
  # True ν reference lines
  geom_hline(
    yintercept = c(0, 1, 2),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  scale_color_manual(
    values = c("#1f77b4", "#2ca02c", "#d62728")
  ) +
  
  labs(
    x = expression(p[reported]),
    y = expression(hat(nu))
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    #axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    
    panel.grid = element_blank(),  # removes all grid lines
    # remove extra space previously reserved for caption
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
#end Fig 5.


##Figure 6 updated with three models including Mhom-het
library(tidyverse)
library(posterior)
getwd()
setwd("C:/Users/2847125J/OneDrive - University of Glasgow/Downloads/BEll/GLASGOW 2023/REsearch/Hetero_2025/Attempt2/Attempt2")
# Load saved objects
draws_df_m11 <- readRDS("draws_model11_heterogeneity.rds")
draws_df_m33 <- readRDS("draws_model33_noheterogeneity.rds")
draws_df_m44 <- readRDS("draws_model44_PNheterogeneity.rds")  # Mhom-het


t_last          <- readRDS("t_last.rds")
cases_matrix    <- readRDS("cases_matrix.rds")
selected_states <- readRDS("selected_states.rds")


#TO SUMMARIZE REUSEABLE
summarise_predictions <- function(draws_df, t_last, model_name) {
  
  n_states <- length(t_last)
  
  t_last_df <- tibble(
    state = 1:n_states,
    t_last = t_last
  )
  
  draws_df %>%
    pivot_longer(
      cols = starts_with("pred_cases["),
      names_to = "var",
      values_to = "value"
    ) %>%
    mutate(
      var = str_remove_all(var, "pred_cases\\[|\\]"),
      state = as.integer(str_extract(var, "^[0-9]+")),
      day   = as.integer(str_extract(var, "(?<=,)[0-9]+"))
    ) %>%
    left_join(t_last_df, by = "state") %>%
    filter(day <= t_last[state]) %>%
    group_by(state, day) %>%
    summarise(
      median = median(value),
      lower  = quantile(value, 0.025),
      upper  = quantile(value, 0.975),
      .groups = "drop"
    ) %>%
    mutate(model = model_name)
}

##SUMMARIZE ALL THREE MODELS

summary_m11 <- summarise_predictions(
  draws_df_m11, 
  t_last, 
  model_name = "Heterogeneity"
)

summary_m33 <- summarise_predictions(
  draws_df_m33, 
  t_last, 
  model_name = "No heterogeneity"
)

summary_m44 <- summarise_predictions(
  draws_df_m44, 
  t_last, 
  model_name = "Mhom-het"
)

summary_all <- bind_rows(summary_m11, summary_m33, summary_m44)

##PREPARE OBSERVED DATA SAHRED BY BOTH MODEL
cases_subset <- as.matrix(cases_matrix[1:max(t_last), ])

obs_df <- purrr::map_dfr(1:length(t_last), function(s) {
  tibble(
    state = s,
    day = 1:(t_last[s] - 1),
    Observed = cases_subset[1:(t_last[s] - 1), s]
  )
})

##MERGE PLOT PREDICTIONS
plot_df <- summary_all %>%
  left_join(obs_df, by = c("state", "day"))

#SINGLE PANEL FOR ALLTHREE
state_labels <- setNames(
  selected_states,
  as.character(1:length(selected_states))
)

ggplot(plot_df, aes(x = day)) +
  geom_point(
    aes(y = Observed),
    color = "black",
    size = 0.8
  ) +
  geom_line(
    #aes(y = median, color = model),
    aes(y = median, color = model,linetype = model),#add line type to show model since overlapping
    linewidth = 1
  ) +
  
  scale_color_manual(
    values = c(
      "Heterogeneity"     = "red",
      "No heterogeneity"  = "blue",
      "Mhom-het"          = "darkgreen"
    )
  ) +
  
  scale_x_continuous(
    limits = c(1, max(t_last)),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(
    ~ state,
    labeller = as_labeller(state_labels),
    scales = "free_y"
  ) +
  labs(
    x = "Day",
    #y = "Observed cases",
    y = "Reported cases",
    color = "Model"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",#IF U WANT CHANGE NON TO RIGHT
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    #axis
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

#end Figure 6


#####NEW Deterministic outbreak comparison Fig. 7 MANSUSCRIPT####
#mean states and four models

library(deSolve)
library(ggplot2)
library(dplyr)
library(grid)

### MODEL: SIR with heterogeneity
SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    incidence <- beta * (S / N)^p * I^q
    
    dS.dt <- -incidence
    dI.dt <-  incidence - gamma * I
    dR.dt <-  gamma * I
    
    return(list(
      c(dS.dt, dI.dt, dR.dt),
      reported = incidence
    ))
  })
}

### FIXED PARAMETERS
gamma <- 0.2
q <- 1
times <- seq(1, 330, by = 1)

### POPULATION
N  <- 5315779
PN <- 0.508

##### 1. Mhom (p = 1)

beta1 <- 0.272
I01   <- 1052.8
S01   <- N - I01

params0 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1,
  q = q,
  N = N
)

out0 <- as.data.frame(
  ode(c(S = S01, I = I01, R = 0),
      times, SIRode, params0)
)

out0$scenario <- "Mhom (p = 1)"


##### 2. Mhet

v1 <- 2.561

params1 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1 + v1^2,
  q = q,
  N = N
)

out1 <- as.data.frame(
  ode(c(S = S01, I = I01, R = 0),
      times, SIRode, params1)
)

out1$scenario <- "Mhet"


##### 3. Mhom-het

beta2 <- 0.272
v2    <- 2.001
I02   <- 728.3

N_eff <- N * PN
S02   <- N_eff - I02

params2 <- list(
  beta = beta2,
  gamma = gamma,
  p = 1 + v2^2,
  q = q,
  N = N_eff
)

out2 <- as.data.frame(
  ode(c(S = S02, I = I02, R = 0),
      times, SIRode, params2)
)

out2$scenario <- "Mhom-het"


##### 4. Mhom (p = 0.508)

S03 <- N_eff - I01

params3 <- list(
  beta = beta1,
  gamma = gamma,
  p = 1,
  q = q,
  N = N_eff
)

out3 <- as.data.frame(
  ode(c(S = S03, I = I01, R = 0),
      times, SIRode, params3)
)

out3$scenario <- "Mhom (p = 0.508)"


##### COMBINE ALL

df_plot <- bind_rows(out0, out1, out2, out3)

#nforce legend order
df_plot$scenario <- factor(
  df_plot$scenario,
  levels = c(
    "Mhom-het",
    "Mhet",
    "Mhom (p = 0.508)",
    "Mhom (p = 1)"
  )
)

##### COLOURS

scenario_colours <- c(
  "Mhom-het"         = "darkgreen",
  "Mhet"             = "red",
  "Mhom (p = 0.508)" = "purple",
  "Mhom (p = 1)"     = "blue"
)

###### PLOT

ggplot(df_plot,
       aes(x = time, y = reported,
           colour = scenario)) +
  
  geom_line(linewidth = 1.2) +
  
  scale_colour_manual(values = scenario_colours) +
  
  labs(
    x = "Time",
    y = "Reported Infected",
    colour = "Scenario"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "right"
  )

# END Figure 7

#### Figure 8: INSEt figure: #Gamma dsn with a mean of 1 and cv=2.561 and Herd Immunity threhsold 

CV <- 2.561
shape <- 1 / (CV^2)
scale <- 1 / shape

x <- seq(0, 20, length.out = 1000)
df_gamma <- data.frame(
  x = x,
  y = dgamma(x, shape = shape, scale = scale)
)

p_inset <- ggplot(df_gamma, aes(x = x, y = y)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.8) +
  labs(
    x = "Susceptibility",
    y = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

# Convert inset to grob
g_inset <- ggplotGrob(p_inset)

# == FINAL PLOT WITH INSET


p_final <- ggplot(hit_df, aes(
  x = v,
  y = HIT,
  group = state
)) +
  
  geom_line(colour = "black", linewidth = 0.9) +
  
  geom_point(
    data = hit_points,
    aes(x = v_hat, y = HIT),
    size = 3,
    colour = "red"
  ) +
  
  labs(
    x = expression("Coefficient of variation (" * v * ")"),
    y = "Herd immunity threshold"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    panel.grid = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  
  # Inset ## to merge both figures
  annotation_custom(
    grob = g_inset,
    xmin = 2.5, xmax = 3.9,
    ymin = 0.08, ymax = 0.5
  ) +
  
  # labeLS 
  #annotate("text", x = 0.2, y = 0.34, label = "(a)", size = 6) +
  annotate("text", x = 0.35, y = 0.36, label = "(a)", size = 6)+ #moves label a slightly up
  annotate("text", x = 2.6, y = 0.48, label = "(b)", size = 5)


p_final

#save
ggsave(
  "Figure9_with_inset.png",
  plot = p_final,
  width = 8,
  height = 6,
  dpi = 300,
  bg = "white"
)#end figure 8

###--  Figure 9: DECLINE IN EFFECTIVE REPRODUCTION NUMBER (mean states) four models (panels), Fig 9 with lines


library(deSolve)
library(dplyr)
library(ggplot2)
library(tidyr)
library(grid)
library(scales)
library(patchwork)

###### SIR MODEL ######

SIRode <- function(t, state, params) {
  with(as.list(c(state, params)), {
    
    incidence <- beta * (S / N)^p * I^q
    
    dS <- -incidence
    dI <-  incidence - gamma * I
    dR <-  gamma * I
    
    list(c(dS, dI, dR))
  })
}

#### PARAMETERS ####

gamma <- 0.2
q <- 1
times <- seq(1, 330, by = 1)

###### MEAN POPULATION ######

N  <- 5315779
PN <- 0.508

##### 1. Mhom (p = 1) #####

beta  <- 0.272  
v     <- 2.561 
I0    <- 1052.8 

S0 <- N - I0
init0 <- c(S = S0, I = I0, R = 0)

R0_basic0 <- beta / gamma
p0 <- 1

out0 <- as.data.frame(
  ode(init0, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = p0,
                   q = q,
                   N = N))
)

out0$model <- "Mhom (p = 1)"
out0$Rt <- R0_basic0 * (out0$S / N)^p0


##### 2. Mhet #####

p1 <- 1 + v^2

out1 <- as.data.frame(
  ode(init0, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = p1,
                   q = q,
                   N = N))
)

out1$model <- "Mhet"
out1$Rt <- R0_basic0 * (out1$S / N)^p1


##### 3. Mhom-het #####

beta2 <- 0.272 
v2    <- 2.001 
I02   <- 728.3 

N_eff2 <- N * PN
S02    <- N_eff2 - I02

init2 <- c(S = S02, I = I02, R = 0)

R0_basic2 <- beta2 / gamma
p2 <- 1 + v2^2

out2 <- as.data.frame(
  ode(init2, times, SIRode,
      parms = list(beta = beta2,
                   gamma = gamma,
                   p = p2,
                   q = q,
                   N = N_eff2))
)

out2$model <- "Mhom-het"
out2$Rt <- R0_basic2 * (out2$S / N_eff2)^p2


##### 4. Mhom (p = 0.508) #####

N_eff0 <- N * PN
S03    <- N_eff0 - I0

init3 <- c(S = S03, I = I0, R = 0)

out3 <- as.data.frame(
  ode(init3, times, SIRode,
      parms = list(beta = beta,
                   gamma = gamma,
                   p = 1,
                   q = q,
                   N = N_eff0))
)

out3$model <- "Mhom (p = 0.508)"
out3$Rt <- R0_basic0 * (out3$S / N_eff0)


##### COMBINE #####

df_all <- bind_rows(out0, out1, out2, out3)

df_all$model <- factor(
  df_all$model,
  levels = c(
    "Mhom-het",
    "Mhet",
    "Mhom (p = 0.508)",
    "Mhom (p = 1)"
  )
)

##### COMPUTE TIME TO Rt = 1 #####

time_to_R1 <- df_all %>%
  arrange(model, time) %>% ##like sorting
  group_by(model) %>% ##each model stands seperately
  summarise(
    t_R1 = {     #creating one number per model
      idx <- which(Rt <= 1)[1]   ##time to which epidemic is declining first outomce is 1
      
      if (is.na(idx) || idx == 1) { #Rt never drops below 1 and already  less than  or equal to 1 at start
        NA   ##both returns NA
      } else {
        t1 <- time[idx - 1]  ##Gets the two points arounnd
        t2 <- time[idx]
        r1 <- Rt[idx - 1]
        r2 <- Rt[idx]
        
        t1 + (1 - r1) * (t2 - t1) / (r2 - r1)
      }
    }
  )

print(time_to_R1)

##### COLOURS #####

scenario_colours <- c(
  "Mhom-het"         = "darkgreen",
  "Mhet"             = "red",
  "Mhom (p = 0.508)" = "purple",
  "Mhom (p = 1)"     = "blue"
)

######## PANEL A ########

p1 <- ggplot(df_all, aes(x = R, y = Rt, colour = model)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_x_continuous(labels = comma) +
  scale_colour_manual(values = scenario_colours) +
  labs(
    x = "Recovered",
    y = expression(R[t]),
    colour = "Scenario",
    title = "A"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

######## PANEL B (UPDATED) ########

p2 <- ggplot(df_all, aes(x = time, y = Rt, colour = model)) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # vertical lines at Rt = 1 crossing
  geom_vline(
    data = time_to_R1,
    aes(xintercept = t_R1, colour = model),
    linetype = "dotted",
    linewidth = 0.8,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(values = scenario_colours) +
  labs(
    x = "Time (days)",
    y = expression(R[t]),
    colour = "Scenario",
    title = "B"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.line = element_line(color = "black"),
    panel.grid = element_blank()
  )

##### SAVE #####

ggsave(
  "Rt_decline_panels1.png",
  p1 / p2,
  width = 8,
  height = 6,
  dpi = 300
)

p1 / p2
#end Figure 9

