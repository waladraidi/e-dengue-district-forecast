# =============================================================================
# Outbreak risk scores from ensemble forecast draws + INLA baselines
# =============================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(ggplot2)

source('./Model/R/99_load.R')


# ---- Config -----------------------------------------------------------------

config <- list(
  baselines_path = "./Model/Data/all_baselines.rds",
  draws_path     = "./Output/Results_summary/draws_from_the_ensemble_fcode_128_RHC.rds",
  output_path    = "./Output/Results_summary/risk_score_and_z_scores_with_prob_corrected.csv",
  n_draws_cap    = 9998,
  horizons       = c(1, 2, 3),
  seed           = 100000
)

set.seed(config$seed)


# ---- Load -------------------------------------------------------------------

d2 <- readRDS('./Model/Data/Full_data_set_with_covariates_and_lags.rds')

obs_cases <- d2 %>%
  mutate(date = as.Date(paste(year, month, "01", sep = "-"))) %>%
  select(fcode, date, obs_dengue_cases, pop_total)

all_baselines <- readRDS(config$baselines_path) %>%
  mutate(date = as.Date(date))

draws <- readRDS(config$draws_path) %>%
  bind_rows() %>%
  mutate(
    date  = as.Date(date),
    cases = exp(value) * pop_total / 100000
  )


# ---- QC: posterior predictive check vs observed cases -----------------------

summarise_posterior_row <- function(log_mean, sd_log, n_draws) {
  # log_mean = mean_log_baseline (log-CASES scale, includes log(offset))
  lambda   <- exp(rnorm(n_draws, log_mean, sd_log))
  historic <- rpois(n_draws, lambda)
  hist_mean <- mean(historic)
  hist_sd   <- sd(historic)
  tibble(
    hist_mean = hist_mean,
    hist_sd   = hist_sd,
    hist_q025 = unname(quantile(historic, 0.025)),
    hist_q500 = unname(quantile(historic, 0.500)),
    hist_q975 = unname(quantile(historic, 0.975)),
    med_thr   = hist_mean + 1 * hist_sd,
    ucl_thr   = hist_mean + 2 * hist_sd
  )
}

summarise_baseline_posterior <- function(baselines, n_draws = 10000) {
  baselines %>%
    rowwise() %>%
    mutate(.summary = list(summarise_posterior_row(
      mean_log_baseline_cases,   # matches column name from baseline_inla_refactored.R
      sd_log_baseline,
      n_draws
    ))) %>%
    ungroup() %>%
    tidyr::unnest(.summary)
}

qc_posterior_plot <- function(fcode_select, posterior_summary) {
  d <- posterior_summary %>% filter(fcode == fcode_select)
  ggplot(d, aes(x = date)) +
    geom_ribbon(aes(ymin = hist_q025, ymax = hist_q975),
                fill = "steelblue", alpha = 0.25) +
    geom_line(aes(y = hist_q500), color = "steelblue") +
    geom_line(aes(y = ucl_thr),   color = "steelblue", lty = 2) +
    geom_line(aes(y = obs_dengue_cases), color = "red") +
    theme_classic() +
    ylab("Cases") +
    ggtitle(sprintf(
      "%s - posterior baseline (blue, 95%% band; dashed = +2 SD) vs observed (red)",
      fcode_select
    ))
}

qc_posterior_sample <- function(baselines, n = 20, n_draws = 10000,
                                plot_dir = "./Output/QC/posterior_plots", seed = 1) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  fcodes  <- unique(baselines$fcode)
  set.seed(seed)
  pick    <- sample(fcodes, min(n, length(fcodes)))
  picked  <- baselines %>% filter(fcode %in% pick)
  summary <- summarise_baseline_posterior(picked, n_draws = n_draws)
  for (fc in pick) {
    p <- qc_posterior_plot(fc, summary)
    ggsave(file.path(plot_dir, sprintf("posterior_qc_%s.png", fc)),
           plot = p, width = 8, height = 4, dpi = 120)
  }
  invisible(pick)
}

# Uncomment to run QC plots:
# qc_posterior_sample(all_baselines, n = 20, n_draws = 1000)


# ---- Risk-score function ----------------------------------------------------

calculate_risk_scores <- function(selected_fcode, selected_date, selected_horizon,
                                  draws, baselines, n_draws_cap = 10000) {

  draws_filt <- draws %>%
    filter(date == selected_date,
           fcode == selected_fcode,
           horizon == selected_horizon) %>%
    slice_head(n = n_draws_cap)

  baseline_filt <- baselines %>%
    filter(date == selected_date, fcode == selected_fcode)

  if (nrow(draws_filt) == 0 || nrow(baseline_filt) == 0) {
    warning(sprintf(
      "Skipping fcode='%s' date=%s horizon=%d (draws=%d baseline=%d)",
      selected_fcode, selected_date, selected_horizon,
      nrow(draws_filt), nrow(baseline_filt)
    ))
    return(NULL)
  }

  n <- nrow(draws_filt)

  # mean_log_baseline is on the log-CASES scale (includes log(offset))
  log_mean_cases <- baseline_filt$mean_log_baseline   # matches refactored baseline
  sd_hist        <- baseline_filt$sd_log_baseline

  draw_df <- tibble(
    forecast = draws_filt$cases,
    lambda   = exp(rnorm(n, log_mean_cases, sd_hist)),
    historic = rpois(n, lambda),
    RR       = (forecast + 1) / (historic + 1)
  )

  hist_mean <- mean(draw_df$historic, na.rm = TRUE)
  hist_sd   <- sd(draw_df$historic,   na.rm = TRUE)
  med_thr   <- hist_mean + 1 * hist_sd
  ucl_thr   <- hist_mean + 2 * hist_sd

  forecast_pmf <- draw_df %>%
    count(forecast, name = "n_obs") %>%
    mutate(probability = n_obs / sum(n_obs))

  prob_in_band <- function(lo, hi) {
    sum(forecast_pmf$probability[forecast_pmf$forecast >  lo &
                                   forecast_pmf$forecast <= hi],
        na.rm = TRUE)
  }

  probs_outbreak <- prob_in_band(ucl_thr,   Inf)
  probs_high     <- prob_in_band(med_thr,   ucl_thr)
  probs_med      <- prob_in_band(hist_mean, med_thr)
  probs_low      <- sum(forecast_pmf$probability[
    forecast_pmf$forecast <= hist_mean], na.rm = TRUE)

  thresholds  <- seq_len(max(draw_df$forecast, na.rm = TRUE))
  probability <- vapply(thresholds, function(t)
    sum(forecast_pmf$probability[forecast_pmf$forecast > t]),
    numeric(1))

  risk_curve <- tibble(
    threshold    = thresholds,
    probability  = probability,
    risk_score   = probability * thresholds,
    risk_score_z = probability * (thresholds - hist_mean) / hist_sd
  )

  i_max   <- which.max(risk_curve$risk_score)
  i_max_z <- which.max(risk_curve$risk_score_z)

  data.frame(
    selected_fcode   = selected_fcode,
    selected_date    = selected_date,
    selected_horizon = selected_horizon,
    prob_RR1_gt_1    = mean(draw_df$RR > 1),

    risk.threshold1                      = risk_curve$risk_score[i_max],
    risk.threshold1_max_point            = risk_curve$threshold[i_max],
    probability_matching_max_threshold   = risk_curve$probability[i_max],

    risk.threshold1z                     = risk_curve$risk_score_z[i_max_z],
    risk.threshold1z_max_point           = risk_curve$threshold[i_max_z],
    probability_matching_max_threshold_z = risk_curve$probability[i_max_z],

    ucl                 = ucl_thr,
    med                 = med_thr,
    historic_mean       = hist_mean,
    probs_outbreak_risk = probs_outbreak,
    probs_high_risk     = probs_high,
    probs_med_risk      = probs_med,
    probs_low_risk      = probs_low
  )
}


# ---- Driver -----------------------------------------------------------------

selected_fcodes   <- unique(draws$fcode)
selected_vintages <- unique(draws$vintage_date)

grid <- expand.grid(
  fcode   = selected_fcodes,
  vintage = selected_vintages,
  horizon = config$horizons,
  stringsAsFactors = FALSE
) %>%
  mutate(
    target_date = as.Date(vintage, origin = "1970-01-01") %m+% months(horizon)
  )

all_results <- vector("list", nrow(grid))
for (k in seq_len(nrow(grid))) {
  all_results[[k]] <- calculate_risk_scores(
    selected_fcode   = grid$fcode[k],
    selected_date    = grid$target_date[k],
    selected_horizon = grid$horizon[k],
    draws            = draws,
    baselines        = all_baselines,
    n_draws_cap      = config$n_draws_cap
  )
}

final_results <- bind_rows(all_results)

dir.create(dirname(config$output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(final_results, config$output_path, row.names = FALSE)

message("Done! Results saved to: ", config$output_path)


# ---- QC: baseline vs observed cases -----------------------------------------

qc_baseline_vs_obs <- function(baselines, obs,
                               fcodes = unique(baselines$fcode)) {
  b <- baselines %>%
    filter(fcode %in% fcodes) %>%
    mutate(
      base_cases     = exp(mean_log_baseline),                        # matches refactored baseline
      base_cases_2sd = exp(mean_log_baseline + 2 * sd_log_baseline)   # matches refactored baseline
    )
  o <- obs %>% filter(fcode %in% fcodes)

  ggplot() +
    geom_line(data = b, aes(x = date, y = base_cases),
              color = "steelblue") +
    geom_line(data = b, aes(x = date, y = base_cases_2sd),
              color = "steelblue", lty = 2) +
    geom_line(data = o, aes(x = date, y = obs_dengue_cases),
              color = "red") +
    facet_wrap(~ fcode, scales = "free_y") +
    theme_classic() +
    ylab("Cases") +
    ggtitle("Baseline (blue, dashed = +2 SD) vs observed (red)")
}

# Uncomment to inspect the first few fcodes:
# print(qc_baseline_vs_obs(all_baselines, obs_cases,
#                          head(unique(all_baselines$fcode), 6)))