##Figure 1: Incidence data for the first wave 
####manuscript edit  After adding MHET-hom model and third talk 11/04/2026
library(tidyverse)
library(plotly)
library(scales)
library(ggforce)
library(zoo)  
library(sf)

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
ggsave("Output/Figure/Figure1.png", width = 10, height = 8)#check for dimension
#END Figure 1


### Fig. 2 MAP OF SELECTED STATES 
# Load the shapefile using sf

nigeria_states <- st_read("Data/Shapefile/gadm41_NGA_1.shp")


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

ggsave("Output/Figure/Figure2_map.png", width = 10, height = 8)# look for this title

###END Fig.2 


