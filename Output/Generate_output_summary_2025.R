library(dplyr)
library(ggplot2)
library(readxl)
library(stringr)
library(spdep)
library(sp)
library(MASS)
library(surveillance)
library(tidyr)
library(tidyverse)
library(sf)
library(scales)
library(dplyr)
library(lubridate)
library(readr)
library(dplyr)
library(writexl)
library(fs)
library(scales)
###Read the quantile summary
## Please note:  
#- Horizon = 1 represents a 1-month ahead prediction.  
#- Horizon = 2 indicates a 2-month ahead forecast.  
#- Horizon = 3 corresponds to a 3-month ahead forecast.  

## The "date" in the quantile summary refers to the target date for which the prediction is made.

quantile_summary <- read.csv('./Output/Results_summary/final_summary_quantiles_fcode_128_RHC_corrected.csv')

quantile_summary <- quantile_summary %>%
  mutate(date = as.Date(date, origin = "1970-01-01"))

###Join the quantile summary with actual data to see dengue cases at that date along with the populaiton
obs_case <- readRDS('./Model/Data/Full_data_set_with_covariates_and_lags.rds') %>%
  dplyr::select(fcode,date, obs_dengue_cases, pop_total)


obs_case <- obs_case %>%
  mutate(date = as.Date(date, origin = "1970-01-01"))

quantiles_and_obs_cases <- left_join(quantile_summary,obs_case,by=c('date','fcode'))

quantiles_and_obs_cases$date <- as.Date(quantiles_and_obs_cases$date,origin = "1970-01-01")

####add also the mean and sd of the historical for each row
all.baselines <- readRDS('./Model/Data/all_baselines.rds')  

set.seed(10000)

all.baselines$baseline_mean <- NA
all.baselines$baseline_sd <- NA

for (i in 1:nrow(all.baselines)) {
  set.historic_log_mean <- all.baselines$mean_log_baseline[i]
  set.historic_log_sd <- all.baselines$sd_log_baseline[i]
  
  historic_samp_mu <- rnorm(10000, mean = set.historic_log_mean, sd = set.historic_log_sd)
  historic_samp <- rpois(10000, lambda = exp(historic_samp_mu))
  
  all.baselines$baseline_mean[i] <- mean(historic_samp)  
  all.baselines$baseline_sd[i] <- sd(historic_samp)
}

quantiles_and_obs_cases<- left_join(quantiles_and_obs_cases,all.baselines,by=c('date','fcode'))

quantiles_and_obs_cases$issue_date <- with(quantiles_and_obs_cases, date %m-% months(horizon))


#### in horizon use 1=1month, 2=2 months, 3= 3months 
quantiles_and_obs_cases <- quantiles_and_obs_cases %>%
  mutate(horizon = case_when(
    horizon == 1 ~ "1 month",
    horizon == 2 ~ "2 months",
    horizon == 3 ~ "3 months",
    TRUE ~ as.character(horizon)
  ))


# Rename columns
quantiles_and_obs_cases <- quantiles_and_obs_cases %>%
  rename(
    pred_date = date,
    Forecast_horizon = horizon,
    pred_cases = mean,
    pred_std=sd,
    pred_cases_95lb = lower_95CI,
    pred_cases_95ub = upper_95CI,
    pred_cases_75lb = lower_75CI,
    pred_cases_75ub = upper_75CI,
    hmean_cases=baseline_mean,
    hsd_cases=baseline_sd
  )

# Check the updated column names
names(quantiles_and_obs_cases)

quantiles_and_obs_cases <- quantiles_and_obs_cases %>%
  dplyr::select(-obs_dengue_cases.x,-obs_dengue_cases.y, -pop_total.y) %>%
  rename(
    pop_total = pop_total.x
  )


# Reorder columns
quantiles_and_obs_cases_subset <- quantiles_and_obs_cases %>%
  dplyr::select(fcode, issue_date, pred_date,Forecast_horizon, hmean_cases, hsd_cases,pred_std,
                pred_cases, pred_cases_95lb, pred_cases_95ub,pred_cases_75lb,pred_cases_75ub, obs_dengue_cases, pop_total)




# Check the updated column order
names(quantiles_and_obs_cases_subset)

########################add the risk score, z-score and probabilities for the year 2023
##selected date here is the forecasted date
#high_risk_threshold: the threshold used for intervention i.e. mean + 2sd
##low_risk_threshold :the low risk threshold (historical mean)
##medium_risk_threshold :the medium risk threshold ( historical mean + 1 std)
##probs_low_risk	P(forecast_dist ≤ mean) — probability that cases fall below or at the historical mean
##probs_med_risk	P(mean < forecast_dist ≤ mean + 1 SD) — probability that cases fall in the medium-risk range
## probs_high_risk	P(mean + 1 SD < forecast_dist ≤ mean + 2 SD) — probability of moderately elevated outbreak risk
## probs_outbreak_risk	P(forecast_dist > mean + 2 SD) — probability that cases are significantly above mean+2SD, indicating outbreak risk

risklevel<- read.csv('./Output/Results_summary/risk_score_and_z_scores_with_prob_corrected.csv')


##this file is generated in fun_threshold_risk_scores.R
risklevel <- risklevel %>%
  mutate(selected_horizon = case_when(
    selected_horizon == 1 ~ "1 month",
    selected_horizon == 2 ~ "2 months",
    selected_horizon== 3 ~ "3 months",
    TRUE ~ as.character(selected_horizon)
  ))




risklevel <- risklevel %>%
  dplyr::select(
    selected_fcode,
    selected_date,
    selected_horizon,
    risk_value_absoulte = risk.threshold1,  # Risk score absolute value
    risk_prob_absoulte = probability_matching_max_threshold,  # Probability P(x >= risk_value)
    risk_value_max_absolute = risk.threshold1_max_point,  # the number of cases/incidence rate corresponding to the maximum risk score
    risk_value_z = risk.threshold1z,  # Z-score risk value
    risk_prob_z = probability_matching_max_threshold_z , # Probability P(x >= risk_value) for z-score
    risk_value_max_z = risk.threshold1z_max_point,  # the number of cases/incidence rate corresponding to the maximum z-risk score
    high_risk_threshold = ucl,
    medium_risk_threshold = med,
    low_risk_threshold = historic_mean,
    outbreak_risk_prob= probs_outbreak_risk,
    high_risk_prob= probs_high_risk,
    medium_risk_prob= probs_med_risk,
    low_risk_prob= probs_low_risk
  )

risklevel$selected_date<- as.Date(risklevel$selected_date, origin = "1970-01-01")

data_summary <- left_join(
  quantiles_and_obs_cases_subset, risklevel , by = c("pred_date" = "selected_date", 
                                                     "Forecast_horizon" = "selected_horizon", 
                                                     "fcode" = "selected_fcode"))
dim(data_summary)


data_summary <- data_summary %>%
  dplyr::select(
    fcode=fcode,
    issue_date=issue_date,
    pred_date,
    pred_std,
    pred_horizon = Forecast_horizon,
    pred_mean = pred_cases,
    obs_dengue_cases=obs_dengue_cases,
    pop_total=pop_total,
    
    pred_75CI_low  = pred_cases_75lb,
    pred_75CI_high = pred_cases_75ub,
    
    pred_95CI_low  = pred_cases_95lb,
    pred_95CI_high = pred_cases_95ub,
    
    pred_risk_abs = risk_value_max_absolute,
    pred_risk_zscore = risk_value_max_z,
    
    epi_threshold = high_risk_threshold,
    prob_outbreak=  outbreak_risk_prob
  )


dim(data_summary)


data_summary <- data_summary %>%
  mutate(year = lubridate::year(issue_date)) %>%
  arrange(issue_date) %>%
  distinct()

##save the results for each file separately 

#write.csv(data_summary,'data_summary.csv')

output_dir <- "./Output/Output_Forecasted_by_Issue_Date"
dir_create(output_dir)

for (yr in unique(data_summary$year)) {
  
  year_folder <- file.path(output_dir, as.character(yr))
  dir_create(year_folder)
  
  year_data <- data_summary %>% 
    filter(year == yr)
  
  unique_dates <- unique(year_data$issue_date)
  
  for (i in seq_along(unique_dates)) {
    
    iss_date <- unique_dates[i]
    subset_data <- year_data %>% filter(issue_date == iss_date)
    
    file_name <- paste0("ed_", format(iss_date, "%Y_%m"), "_pred.csv")
    file_path <- file.path(year_folder, file_name)
    
    write_csv(subset_data, file_path)
    print(paste("Saved:", file_path))
  }
}


##plot one fcode

df_plot <- data_summary %>%
  filter(fcode == "TTYTKV_BA_TRI_BT",
         pred_horizon == "3 months") 

ggplot(df_plot, aes(x = pred_date)) +
  geom_line(aes(y = obs_dengue_cases), linewidth = 0.9) +
  geom_line(aes(y = pred_mean), linewidth = 0.9, color = "red") +
  labs(
    title = "Observed vs Predicted Mean (Horizon = 3 months) - BA_TRI",
    x = "Prediction target date",
    y = "Dengue cases"
  ) +
  theme_bw()

