# =============================================================================
# Baseline estimation via INLA seasonal decomposition  (refactored)
# -----------------------------------------------------------------------------
# For each (fcode, forecast_year), fits a Poisson INLA model with an AR(1)
# short-term trend, an RW(2) long-term trend, and a cyclic RW(1) seasonal
# component. The "baseline" is the linear predictor MINUS the AR(1) part,
# i.e. long-run + seasonal expectation only. Saved as one RDS per fit:
#     ./Output/Results/baselines/baseline_<year>_<fcode>.rds
# =============================================================================

library(tidyverse)
library(INLA)

source('./Model/R/99_load.R')


# ---- Config -----------------------------------------------------------------

config <- list(
  output_dir           = "./Output/Results/baselines",
  combined_output_path = "./Model/Data/all_baselines.rds",
  fit_start_date       = as.Date("2004-09-01"),    # earliest date used for fitting
  forecast_extend_from = as.Date("2026-05-01"),    # extend skeleton to cover
  forecast_extend_to   = as.Date("2026-08-01"),    #   forecast months with NA cases
  hyper_rw1            = list(prec = list(prior = "pc.prec", param = c(0.3, 0.01))),
  inla_threads         = 8,
  inla_seed            = 8123
)

dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)


# ---- Load data --------------------------------------------------------------

d2 <- readRDS('./Model/Data/Full_data_set_with_covariates_and_lags.rds')

ds <- d2 %>%
  mutate(
    date    = as.Date(paste(year, month, "01", sep = "-")),
    offset1 = pop_total / 100000
  )


# ---- Extend skeleton for future forecast months -----------------------------

extend_skeleton <- function(ds, from_date, to_date) {
  new_months <- seq.Date(from_date, to_date, by = "month")
  
  latest_offset <- ds %>%
    group_by(fcode) %>%
    filter(date == max(date)) %>%
    dplyr::select(fcode, offset1) %>%
    ungroup()
  
  new_rows <- ds %>%
    distinct(fcode) %>%
    crossing(date = new_months) %>%
    left_join(latest_offset, by = "fcode")
  
  bind_rows(ds, new_rows) %>%
    arrange(fcode, date) %>%
    dplyr::select(fcode, date, obs_dengue_cases, offset1)
}

ds_ext <- extend_skeleton(
  ds,
  from_date = config$forecast_extend_from,
  to_date   = config$forecast_extend_to
)


# ---- Baseline INLA model ----------------------------------------------------

prepare_panel <- function(data, fit_start_date) {
  data %>%
    filter(date >= fit_start_date) %>%
    arrange(fcode, date) %>%
    mutate(
      t        = lubridate::interval(min(date), date) %/% months(1) + 1,
      time_id1 = t - min(t, na.rm = TRUE) + 1,
      monthN   = month(date),
      year     = lubridate::year(date),
      yearN    = as.numeric(as.factor(year))
    )
}

fit_baseline_inla <- function(panel, fcode_select, forecast_year, hyper_rw1,
                              num_threads = 8) {
  
  c2 <- panel %>%
    filter(fcode == fcode_select, year <= forecast_year) %>%
    mutate(obs_dengue_cases_fit = if_else(year >= forecast_year,
                                          NA_real_, obs_dengue_cases))
  
  if (nrow(c2) == 0) {
    warning(sprintf("No data for fcode='%s' year<=%d - skipping.",
                    fcode_select, forecast_year))
    return(NULL)
  }
  
  form <- obs_dengue_cases_fit ~ 1 +
    f(time_id1, model = "ar1", constr = TRUE) +
    f(yearN,    model = "rw2", constr = TRUE) +
    f(monthN,   model = "rw1", hyper = hyper_rw1, cyclic = TRUE,
      scale.model = TRUE, constr = TRUE)
  
  year_mat  <- model.matrix(~ -1 + as.factor(yearN),  data = c2)
  month_mat <- model.matrix(~ -1 + as.factor(monthN), data = c2)
  
  lc_no_ar1 <- INLA::inla.make.lincombs(
    "(Intercept)" = rep(1, nrow(c2)),
    "time_id1"    = rep(0, nrow(c2)),
    "yearN"       = year_mat,
    "monthN"      = month_mat
  )
  
  mod <- INLA::inla(
    form,
    data    = c2,
    family  = "poisson",
    E       = c2$offset1,
    lincomb = lc_no_ar1,
    control.compute   = list(dic = FALSE, waic = FALSE, config = TRUE,
                             return.marginals = FALSE),
    control.predictor = list(compute = TRUE, link = 1),
    control.fixed     = list(mean.intercept = 0, prec.intercept = 1e-4,
                             mean = 0, prec = 1),
    inla.mode   = "experimental",
    num.threads = num_threads
  )
  
  mod$summary.lincomb.derived %>%
    transmute(
      mean_log_baseline_cases = mean + log(c2$offset1),
      sd_log_baseline         = sd
    ) %>%
    bind_cols(tibble(
      date             = c2$date,
      fcode            = fcode_select,
      obs_dengue_cases = c2$obs_dengue_cases,
      year             = c2$year
    )) %>%
    filter(year == forecast_year) %>%
    dplyr::select(date, fcode, obs_dengue_cases, mean_log_baseline_cases, sd_log_baseline)
}

baseline_path <- function(fcode, forecast_year, dir) {
  file.path(dir, sprintf("baseline_%d_%s.rds", forecast_year, fcode))
}


# ---- Driver -----------------------------------------------------------------

run_baselines <- function(data, fcodes, forecast_years, config) {
  set.seed(config$inla_seed)
  panel <- prepare_panel(data, config$fit_start_date)
  
  for (yr in forecast_years) {
    for (fc in fcodes) {
      message(sprintf("Fitting baseline: year=%d  fcode=%s", yr, fc))
      out <- fit_baseline_inla(panel, fc, yr,
                               config$hyper_rw1, config$inla_threads)
      if (!is.null(out)) {
        saveRDS(out, baseline_path(fc, yr, config$output_dir))
      }
    }
  }
}

forecast_years <- 2026:max(ds$year)   # matches original script loop start
fcodes         <- unique(ds$fcode)

run_baselines(ds_ext, fcodes, forecast_years, config)


# ---- Combine per-fcode baselines into a single file -------------------------

combine_baselines <- function(forecast_years, dir) {
  pattern <- sprintf("^baseline_(%s)_.+\\.rds$",
                     paste(forecast_years, collapse = "|"))
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  bind_rows(lapply(files, readRDS))
}

all_baselines <- combine_baselines(forecast_years, config$output_dir)
saveRDS(all_baselines, config$combined_output_path)

message("Done! Combined baselines saved to: ", config$combined_output_path)


# ---- QC plots ---------------------------------------------------------------

qc_plot <- function(forecast_year, fcode_select, config) {
  out <- readRDS(baseline_path(fcode_select, forecast_year, config$output_dir)) %>%
    mutate(
      base_cases     = exp(mean_log_baseline_cases),
      base_cases_2sd = exp(mean_log_baseline_cases + 2 * sd_log_baseline)
    )
  
  ggplot(out, aes(x = date)) +
    geom_line(aes(y = base_cases)) +
    geom_line(aes(y = base_cases_2sd), lty = 2, color = "gray") +
    geom_line(aes(y = obs_dengue_cases), color = "red") +
    theme_classic() +
    ggtitle(fcode_select)
}

qc_plot_sample <- function(forecast_year, config, n = 20,
                           plot_dir = "./Output/QC/baseline_plots", seed = 1) {
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  pattern <- sprintf("^baseline_%d_(.+)\\.rds$", forecast_year)
  files   <- list.files(config$output_dir, pattern = pattern)
  fcodes  <- sub(pattern, "\\1", files)
  
  set.seed(seed)
  pick <- sample(fcodes, min(n, length(fcodes)))
  
  for (fc in pick) {
    p <- qc_plot(forecast_year, fc, config)
    ggsave(file.path(plot_dir, sprintf("qc_%d_%s.png", forecast_year, fc)),
           plot = p, width = 8, height = 4, dpi = 120)
  }
  invisible(pick)
}

# Uncomment to run QC:
# qc_plot_sample(2026, config)