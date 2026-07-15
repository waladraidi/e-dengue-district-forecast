##New models set 


# "flexible harmonics" family 
# ---------------------------------------------------------

# mod1: MOST FLEXIBLE - Full hierarchical seasonality
# Captures seasonality at 3 levels (district, year, year×district) with both annual + semi-annual cycles
mod1 <- '
  obs_dengue_cases_hold ~ 1 +
  lag3_y +
  sin12 + cos12 +
  f(fcodeID0, model="iid") +
  # District seasonality (annual + semi-annual)
  f(fcodeID1, sin12, model="iid") +
  f(fcodeID2, cos12, model="iid") +
  f(fcodeID3, sin6,  model="iid") +
  f(fcodeID4, cos6,  model="iid") +
  f(fcodeID5, lag3_y, model="iid") +
  # Year seasonality (annual + semi-annual)
  f(yearID1, sin12, model="iid") +
  f(yearID2, cos12, model="iid") +
  f(yearID3, sin6,  model="iid") +
  f(yearID4, cos6,  model="iid") +
  # Year×district
  f(year_fcode_ID1, sin12, model="iid") +
  f(year_fcode_ID2, cos12, model="iid") +
  f(year_fcode_ID3, sin6,  model="iid") +
  f(year_fcode_ID4, cos6,  model="iid")
'

# mod2: SIMPLEST - Basic district seasonality only
# Only annual cycles at district level, no year effects or semi-annual patterns
mod2 <-  '
     obs_dengue_cases_hold  ~ 1 +
     lag3_y +
    sin12 + cos12 +
    f(fcodeID0, model = "iid") +
    
    # District seasonality - annual
    f(fcodeID1, sin12, model = "iid") +
    f(fcodeID2, cos12, model = "iid") 
    '

# mod3: CLIMATE-FOCUSED - Adds climate (Thermal) + cumulative incidence
mod3<- '
  obs_dengue_cases_hold ~ 1 +
  climate_scale_lag3 +
  log_cum_inc_12m + log_cum_inc_24m + log_cum_inc_36m +
  lag3_y +
  sin12 + cos12 +
  f(fcodeID0, model="iid") +
  f(fcodeID1, sin12, model="iid") +
  f(fcodeID2, cos12, model="iid") +
  f(fcodeID3, climate_scale_lag3, model="iid") +
  f(monthN, model="rw1", hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE)
'

# mod4: AR2 TEMPORAL - Second-order autoregression
# Uses AR(2) process to capture more complex temporal momentum (2 lag periods)
# Includes spatial (BESAG) + climate

mod4 <- 'obs_dengue_cases_hold ~ lag3_y + log_cum_inc_12m + log_cum_inc_24m + log_cum_inc_36m +
  f(fcodeID, model="besag", constr=TRUE, graph=MDR.adj, hyper=hyper.besag, scale.model=TRUE) +
  lag3_avg_min_daily_temp + lag3_monthly_cum_ppt +
  f(t, replicate=fcodeID3, model="ar", order=2, hyper=hyper.ar2, constr=TRUE) +  # AR(2) instead of AR(1)
  f(monthN, model="rw1", hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE)'


# mod5: DISTRICT-SPECIFIC AR1 - Dual temporal structure
# Combines global AR(1) with district-specific AR(1) effects
mod5 <- 'obs_dengue_cases_hold ~ lag3_y + log_cum_inc_12m + log_cum_inc_24m + log_cum_inc_36m +
  f(fcodeID, model="besag", constr=TRUE, graph=MDR.adj, hyper=hyper.besag, scale.model=TRUE) +
  lag3_avg_min_daily_temp + lag3_monthly_cum_ppt +
  f(t, model="ar1", hyper=hyper.ar1, constr=TRUE) +  # Global AR1
  f(t2, replicate=fcodeID, model="ar1", hyper=hyper.ar1, constr=TRUE) +  # District-specific AR1
  f(monthN, model="rw1", hyper=hyper2.rw, cyclic=TRUE, scale.model=TRUE, constr=TRUE)'


