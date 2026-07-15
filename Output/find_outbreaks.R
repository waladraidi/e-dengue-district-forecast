# =============================================================================
# Alert detection across incidence thresholds
# Adds randomization group: Intervention / Control / Excluded
# =============================================================================

library(dplyr)
library(readr)
library(lubridate)

# ---- Load data --------------------------------------------------------------

obs_case <- readRDS('./Model/Data/Full_data_set_with_covariates_and_lags.rds') %>%
  dplyr::select(fcode, date, obs_dengue_cases, pop_total)

risk_scores <- read.csv('./Output/Results_summary/risk_score_and_z_scores_with_prob_corrected.csv') %>%
  mutate(
    selected_date    = as.Date(selected_date),
    selected_horizon = as.numeric(selected_horizon),
    issue_date       = selected_date %m-% months(selected_horizon),
    selected_horizon = case_when(
      selected_horizon == 1 ~ "1 month",
      selected_horizon == 2 ~ "2 months",
      selected_horizon == 3 ~ "3 months",
      TRUE ~ as.character(selected_horizon)
    )
  )

quantile_summary <- read.csv('./Output/Results_summary/final_summary_quantiles_fcode_128_RHC_corrected.csv') %>%
  mutate(
    date    = as.Date(date),
    horizon = case_when(
      horizon == 1 ~ "1 month",
      horizon == 2 ~ "2 months",
      horizon == 3 ~ "3 months",
      TRUE ~ as.character(horizon)
    )
  )

randomization <- read.csv('./Output/Results_summary/ED Randomization result_26022026_FINAL.csv') %>%
  dplyr::select(fcode, Randomization) %>%
  mutate(
    group = case_when(
      Randomization == "I" ~ "Intervention",
      Randomization == "C" ~ "Control",
      TRUE ~ "Excluded"
    )
  ) %>%
  dplyr::select(fcode, group)

# ---- Join all data ----------------------------------------------------------

d <- risk_scores %>%
  inner_join(quantile_summary,
             by = c("selected_date" = "date",
                    "selected_fcode" = "fcode",
                    "selected_horizon" = "horizon")) %>%
  inner_join(obs_case,
             by = c("selected_date" = "date",
                    "selected_fcode" = "fcode")) %>%
  # Remove duplicate columns from join
  dplyr::select(-obs_dengue_cases.y, -pop_total.y) %>%
  rename(obs_dengue_cases = obs_dengue_cases.x,
         pop_total        = pop_total.x) %>%
  # Add randomization group
  left_join(randomization, by = c("selected_fcode" = "fcode")) %>%
  mutate(group = if_else(is.na(group), "Excluded", group)) %>%
  mutate(
    pred_mean_inc        = mean / pop_total * 100000,
    ensemble_dist_wgt    = mean,
    ensemble_dist_wgt_inc = mean / pop_total * 100000,
    year  = lubridate::year(selected_date),
    month = lubridate::month(selected_date)
  )

# ---- Generate alerts for each incidence threshold ---------------------------
# ---- Generate alerts for each incidence threshold ---------------------------

dir.create("./Output/Results_summary/alerts_by_incidence",
           showWarnings = FALSE)

thresholds <- c(0,10, 15, 20, 25, 30, 35, 40, 45, 50)

#thresholds <- c(20)
all_alerts <- list()

for (th in thresholds) {
  
  a_flagged <- d %>%
    group_by(selected_date, issue_date, selected_fcode, selected_horizon) %>%
    mutate(
      threshold_incidence = th,
      
      epidemic_flag_risk      = as.numeric(risk.threshold1_max_point > ucl),
      epidemic_flag_z_score   = as.numeric(risk.threshold1z_max_point > ucl),
      epidemic_flag_pred_mean = as.numeric(ensemble_dist_wgt > ucl),
      
      epidemic_methods_count = rowSums(
        across(c(epidemic_flag_risk,
                 epidemic_flag_z_score,
                 epidemic_flag_pred_mean)),
        na.rm = TRUE
      ),
      
      epidemic_flag_from_two_methods = as.numeric(
        epidemic_methods_count >= 2 & pred_mean_inc > th),
      
      epidemic_flag_from_three_methods = as.numeric(
        epidemic_methods_count >= 3 & pred_mean_inc > th)
    ) %>%
    ungroup()
  
  alerts_only <- a_flagged %>%
    dplyr::select(
      threshold_incidence,
      issue_date,
      selected_date,
      selected_fcode,
      selected_horizon,
      group,
      year,
      epidemic_flag_risk,
      epidemic_flag_z_score,
      epidemic_flag_pred_mean,
      epidemic_flag_from_two_methods,
      epidemic_flag_from_three_methods,
      mean,
      median,
      obs_dengue_cases,
      pop_total
    ) %>%
    distinct()
  
  # Save each threshold separately
  write_csv(
    alerts_only,
    paste0("./Output/Results_summary/alerts_by_incidence/alerts_by_incidence", th, ".csv")
  )
  
  # Store for combined file
  all_alerts[[as.character(th)]] <- alerts_only
  
  message(sprintf("Threshold %d: %d rows, %d alerts (3-method)",
                  th, nrow(alerts_only),
                  sum(alerts_only$epidemic_flag_from_three_methods, na.rm = TRUE)))
}

# ---- Bind all thresholds together ------------------------------------------

alerts_all_thresholds <- bind_rows(all_alerts)

# ---- Save one combined file ------------------------------------------------

write_csv(
  alerts_all_thresholds,
  "./Output/Results_summary/alerts_by_incidence/alerts_by_incidence_all_thresholds.csv"
)

message("Done. Separate files and combined file saved.")