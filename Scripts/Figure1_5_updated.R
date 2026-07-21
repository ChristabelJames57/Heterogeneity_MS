#Figure 1: Incidence data for the first wave #####
library(tidyverse)
library(plotly)
library(ggforce)
library(zoo)
library(sf)
library(posterior)
library(grid)
library(deSolve)
library(scales)
library(patchwork)

rm(list = ls()) # Clear workspace

##Figure_1 
Data_wave1 <- read_csv("./Data/wave12_data.csv")#new

###Data has day, month and year can create a proper date column
Data_wave1 <- Data_wave1 %>%
  mutate(Date = make_date(year, month, day)) 
names(Data_wave1)


# Aggregate data by date

D_State <- Data_wave1 %>%
  group_by(Date,`State of Residence`) %>% #can do multiple group by of variables
  summarize(Case_Count = n()) %>% 
  na.omit(`State of Residence`) #To remove the NA in state

#To check if any NA
any(is.na(D_State))# in the entire data frame of D_state. if NA will return true if no will return False
any(is.na(D_State$date))# to check for NA in date column
any(is.na(D_State$Case_Count))#check for NA in the case_count column


#Aggregate Cases by date and state
cases_by_state <- Data_wave1 %>% ##Shows the data frame to use for plotting
  group_by(State = `State of Residence`, Date) %>%
  summarise(Case_Count = n(), .groups = "drop")

###clean
cases_by_state_clean <- na.omit(cases_by_state)#remove NA before ploting


# Define a vector with the specific states to plot
selected_states <- c("Delta", "Ebonyi", "Edo", "Enugu", "Kaduna", "Kwara", "Ondo", "Osun", "Oyo", "Rivers")  # Replace with the states you're interested in


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
ggsave("Output/Figure/Figure_1.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_1.png", width = 7, height = 5, bg = "white", dpi=1000)
#END Figure 1

# Figure_2 MAP OF SELECTED STATES ######
# Load the shapefile using sf

nigeria_states <- st_read("Data/Shapefile/gadm41_NGA_1.shp")


#prepare data
# Check the column names
print(names(nigeria_states))

# Clean and format the state names
nigeria_states$name <- tools::toTitleCase(tolower(nigeria_states$NAME_1))

selected_states <- c("Ondo", "Edo", "Oyo", "Delta", "Enugu", "Kaduna", "Ebonyi",  
                     "Kwara", "Rivers", "Osun")#10 states updated

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

ggsave("Output/Figure/Figure_2.pdf", width = 7, height = 5, bg = "white")
###END Fig.2 



#Figure_3 Recovery of R0 parameter#####

# Load data
#df_R0 <- read_csv("./Data/MeanR0_summaryNEW.csv")
df_R0 <- read_csv("./Data/MeanR0_summaryNEW_Added.csv")#updated

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
      "2" = "red",
      "3" = "black"#added
    )
  ) +
  
  labs(
    #x = expression(p[reported]),
    # y = expression(hat(R)[0])
    x = expression(italic(r)),
    y = expression(hat(R)[0])
  ) +
  
  scale_y_continuous(
    limits = c(1.93, 2.05),
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
ggsave("Output/Figure/Figure_3.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_3.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_3



#Figure_4 Recovery of r parameter##### 

# Read data
#df <- read_csv("./Data/MeanPreported_summaryNEW.csv")
df <- read_csv("./Data/MeanPreported_summaryNEW_Added.csv")#updated
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
      "2" = "red",
      "3" = "black"
    )
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(r))
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
# Save
ggsave("Output/Figure/Figure_4.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_4.png", width = 7, height = 5, bg = "white",dpi=1000)

##end Figue_4 



#Figure_5 Recovery of v parameter#####

# Read data
#df <- read_csv("./Data/MeanV_summary3.csv")
df <- read_csv("./Data/MeanV_summary3_Added.csv")  # updated

# Treat v_true as a factor with fixed levels
df <- df %>%
  mutate(v_true = factor(v_true, levels = c(0, 1, 2, 3)))

# Plot
p <- ggplot(df, aes(x = p_r, y = MeanV, color = v_true, group = v_true)) +
  
  geom_line(linewidth = 1) +
  
  geom_point(size = 3) +
  
  geom_errorbar(
    aes(ymin = LowerCI, ymax = UpperCI),
    width = 0.0005,
    linewidth = 0.8
  ) +
  
  # True ν reference lines
  geom_hline(
    yintercept = c(0, 1, 2, 3),
    linetype = "dashed",
    color = "gray50"
  ) +
  
  # Colours:
  # v = 0 → Blue
  # v = 1 → Green
  # v = 2 → Red
  # v = 3 → Black
  scale_color_manual(
    values = c(
      "0" = "#1f77b4",
      "1" = "#2ca02c",
      "2" = "#d62728",
      "3" = "black"
    )
  ) +
  
  labs(
    x = expression(italic(r)),
    y = expression(hat(nu))
  ) +
  
  theme_minimal(base_size = 14) +
  
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 13),
    
    # Axes
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.25, "cm"),
    
    # Remove grid
    panel.grid = element_blank(),
    
    # Margins
    plot.margin = margin(
      t = 10,
      r = 10,
      b = 10,
      l = 10
    )
  )
print(p)
# Save
ggsave("Output/Figure/Figure_5.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_5.png", width = 7, height = 5, bg = "white", dpi=1000)
#end Figure_5.


#R0 — Bias, Coverage, RMSE, Uncertainty#####

true_R0 <- 2   # known simulation truth

df_R0_stats <- df_R0 %>%
  mutate(
    Bias        = MeanR0 - true_R0,
    Coverage    = ifelse(LowerCI <= true_R0 & UpperCI >= true_R0, 1, 0),
    RMSE        = sqrt((MeanR0 - true_R0)^2),
    Uncertainty = UpperCI - LowerCI
  ) %>%
  group_by(v, p_r) %>%
  summarise(
    MeanR0      = mean(MeanR0),
    Bias        = mean(Bias),
    Coverage    = mean(Coverage),
    RMSE        = mean(RMSE),
    Uncertainty = mean(Uncertainty),
    .groups = "drop"
  )

write.csv(df_R0_stats, "Output/Figure/R0_summary_stats.csv", row.names = FALSE)
#end R0


#r — Bias, Coverage, RMSE, Uncertainty#####

df_r_stats <- df %>%
  mutate(
    Bias        = Meanp_r - p_r_true,
    Coverage    = ifelse(LowerCI <= p_r_true & UpperCI >= p_r_true, 1, 0),
    RMSE        = sqrt((Meanp_r - p_r_true)^2),
    Uncertainty = UpperCI - LowerCI
  ) %>%
  group_by(v, p_r) %>%
  summarise(
    Meanp_r     = mean(Meanp_r),
    Bias        = mean(Bias),
    Coverage    = mean(Coverage),
    RMSE        = mean(RMSE),
    Uncertainty = mean(Uncertainty),
    .groups = "drop"
  )

write.csv(df_r_stats, "Output/Figure/r_summary_stats.csv", row.names = FALSE)
#end r


#v — Bias, Coverage, RMSE, Uncertainty#####
df_v <- read_csv("./Data/MeanV_summary3_Added.csv")

df_v_stats <- df_v %>%
  mutate(
    Bias        = MeanV - v_true,
    Coverage    = ifelse(LowerCI <= v_true & UpperCI >= v_true, 1, 0),
    RMSE        = sqrt((MeanV - v_true)^2),
    Uncertainty = UpperCI - LowerCI
  ) %>%
  group_by(v_true, p_r) %>%
  summarise(
    MeanV       = mean(MeanV),
    Bias        = mean(Bias),
    Coverage    = mean(Coverage),
    RMSE        = mean(RMSE),
    Uncertainty = mean(Uncertainty),
    .groups = "drop"
  )

write.csv(df_v_stats, "Output/Figure/v_summary_stats.csv", row.names = FALSE)
#end v
#Output/Figure/R0_summary_stats.csv
#Output/Figure/r_summary_stats.csv
#Output/Figure/v_summary_stats.csv


#trial
R0_out <- df_R0_stats %>%
  rename(Estimate = MeanR0) %>%
  mutate(
    Parameter = "R0",
    TrueValue = 2
  )

r_out <- df_r_stats %>%
  rename(
    Estimate = Meanp_r,
    TrueValue = p_r
  ) %>%
  mutate(Parameter = "r")

v_out <- df_v_stats %>%
  rename(
    Estimate = MeanV,
    TrueValue = v_true
  ) %>%
  mutate(Parameter = "v")

summary_stats <- bind_rows(R0_out, r_out, v_out)

write.csv(
  summary_stats,
  "Output/Figure/Combined_summary_stats.csv",
  row.names = FALSE
)


#all in one#####
# Add parameter labels
df_R0_stats2 <- df_R0_stats %>%
  mutate(Parameter = "R0")

df_r_stats2 <- df_r_stats %>%
  mutate(Parameter = "r")

df_v_stats2 <- df_v_stats %>%
  mutate(Parameter = "v")

# Combine all three
combined_stats <- bind_rows(df_R0_stats2, df_r_stats2, df_v_stats2)

# Save
write.csv(combined_stats, "Output/Figure/combined_summary_stats.csv", row.names = FALSE)

#Output/Figure/combined_summary_stats.csv

