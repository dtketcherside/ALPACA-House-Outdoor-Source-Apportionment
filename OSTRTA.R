# Load Libraries ----------
library(dplyr)
library(data.table)
library(lubridate)
library(tidyr)
library(ggplot2)
library(patchwork)
library(openair)
library(scales)
library(magick)
library(grid)

# Load Data ---------
## Cold start vehicle exhaust experiments not removed for some species.
cold_start_exp_times <- as.POSIXct(c('2022-02-23 20:00:00','2022-02-24 19:00:00'),tz = 'America/Sitka')

CTC_data <- fread('data/ALPACA_Outdoor_1min_20240607.txt',tz='') %>%
  mutate_at(vars(Time_AKST),~ymd_hms(.,tz = 'America/Sitka',truncated = 3)) %>%
  filter(!between(Time_AKST,cold_start_exp_times[1],cold_start_exp_times[1]+minutes(59)+seconds(59))) %>%
  filter(!between(Time_AKST,cold_start_exp_times[2],cold_start_exp_times[2]+minutes(59)+seconds(59)))

NO3_rad_data <- fread('data/PACT1D_ALPACA_wSnow3m.txt',tz='',check.names = T) %>%
  mutate_at(vars(datetime..AKST.),~dmy_hms(.,tz = 'America/Sitka',truncated = 3))
colnames(NO3_rad_data) <- gsub('..molec.cm..3.','',colnames(NO3_rad_data))
colnames(NO3_rad_data)[1] <- 'Time_AKST'

NO3_prod_data <- fread('data/NO3 Production Rates.csv',tz='') %>%
  mutate_at(vars(Time_AKST),~mdy_hms(.,tz = 'America/Sitka',truncated = 3))

hr_data <- fread('data/time_series.csv',tz='') %>%
  mutate_at(vars(Time_AKST),~mdy_hm(.,tz = 'America/Sitka',truncated = 3)) %>%
  filter(!between(Time_AKST,cold_start_exp_times[1],cold_start_exp_times[1]+minutes(59)+seconds(59))) %>%
  filter(!between(Time_AKST,cold_start_exp_times[2],cold_start_exp_times[2]+minutes(59)+seconds(59)))

met_data <- fread('data/PAFA.csv',tz='') %>%
  mutate_at(vars(Date_Time),~mdy_hm(.,tz = 'America/Sitka',truncated = 3)) %>% 
  filter(between(Date_Time,min(hr_data$Time_AKST),max(hr_data$Time_AKST)))
colnames(met_data)[1] <- 'Time_AKST'

PMF_TS_data <- fread('data/PMF_TS_contributions.csv',tz='') %>%
  mutate_at(vars(Time_AKST),~mdy_hms(.,tz = 'America/Sitka',truncated = 3))

PMF_spec_cont <- fread('data/PMF_species_perc_dist.csv')

tabs4 <- fread('data/Table S4 - Factor Contributions to Species Fit.csv') %>%
  pivot_longer(!Species_Name,names_to = 'PMF_factor', values_to = 'perc_cont')

SI_plot_data <- fread('data/Figure S3 Data.csv') %>%
  pivot_wider(id_cols = Species,names_from = data_ID,values_from = composition)

# Load Custom Functions ----------
source('Custom Functions.R')
# Run Project ------
source('Make and Save Figures.R')