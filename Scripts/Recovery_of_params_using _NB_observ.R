
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

#Figure_S9 Recovery of R0 parameter from Nrgative Binom observation#####

# Load data
#df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_use2.csv")#gam(1,10)for Normal and NB
#df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_Gamma.csv")#using gama (1,1)
df_R0 <- read_csv("./Data/MeanR0_summaryNEW_NB_Gam01.csv")#using gama (0.1,0.1)

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
    #x = expression(p[reported]),
    # y = expression(hat(R)[0])
    x = expression(italic(r)),
    y = expression(hat(R)[0])
  ) +
  
  scale_y_continuous(
    #limits = c(1.80, 2.05),
    limits = c(1.80, 2.06),
    breaks = c(1.95, 2.00, 2.06)
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
ggsave("Output/Figure/Figure_S10a.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_S10a.png", width = 7, height = 5, bg = "white", dpi=1000)

#ggsave("Output/Figure/Figure_S9a.pdf", width = 7, height = 5, bg = "white")
#ggsave("Output/Figure/Figure_S9a.png", width = 7, height = 5, bg = "white", dpi=1000)

#end Figure_10



#Figure_11 Recovery of r  parameter from NEgative Binom observation#####

# Read data
df <- read_csv("./Data/MeanPreported_summaryNEW_NB_use2.csv")#0 to 3 Gamma(1,10)
#df <- read_csv("./Data/MeanPreported_summaryNEW_NB_Gamma.csv")# 0 to 2 Gamma(1,1)
#df <- read_csv("./Data/MeanPreported_summaryNEW_NB_Gam01.csv")# 0 to 2 Gamma(0.1,0.1)

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
#ggsave("Output/Figure/Figure_11_NB_01_01.pdf", width = 7, height = 5, bg = "white")
#ggsave("Output/Figure/Figure_11_NB_01_01.png", width = 7, height = 5, bg = "white", dpi=1000)

#ggsave("Output/Figure/Figure_S9b.pdf", width = 7, height = 5, bg = "white")
#ggsave("Output/Figure/Figure_S9b.png", width = 7, height = 5, bg = "white", dpi=1000)


#Figure_12 Recovery of V  parameter from NEgative Binom observation#####

# Read data
#df <- read_csv("./Data/MeanV_summary3_NB_use2.csv")#stops at 2 Nb or norm Gammma(1,10)
#df <- read_csv("./Data/MeanV_summary3_NB_Gamma.csv")#stops at 2 Gamma(1,1)
df <- read_csv("./Data/MeanV_summary3_NB_Gam01.csv")#stops at 2 Gamma(0.1,0.1)


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
    
    values = c("#1f77b4", "#2ca02c", "#d62728")#0 to 2
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
# Save
ggsave("Output/Figure/Figure_S10b.pdf", width = 7, height = 5, bg = "white")
ggsave("Output/Figure/Figure_S10b.png", width = 7, height = 5, bg = "white", dpi=1000)

#ggsave("Output/Figure/Figure_S9b.pdf", width = 7, height = 5, bg = "white")
#ggsave("Output/Figure/Figure_S9b.png", width = 7, height = 5, bg = "white", dpi=1000)

#end Figure_12
