library(parallel)
library(stats)
library(dplyr)
library(readr)
library(readxl)
library(tidyr)
library(tidyverse)  
library(zoo)
library(pbapply)
library(MASS)
library(scoringutils)
library(sf)
library(spdep)
library(lattice)
library(stringr)
library(janitor)
library(surveillance)
library(lubridate)
library(INLA)


source('./Model/R/99_helper_funcs.R')
source('./Model/R/99_define_inla_spacetime_mods.R')
source('./Model/R/01_fun_inla_spacetime.R')
source('./Model/R/02_fun_hhh4.R')
#source('./Model/R/03_fun_lag_district_pca.R')

##read the file from "Input" folder
d2<-  readRDS('./Model/Data/Full_data_set_with_covariates_and_lags.rds')


MDR_NEW <- readRDS( './Model/Data/MDR_NEW.rds')


spat_IDS <- readRDS( "./Model/Data/spatial_IDS.rds")

# Create initial neighbor structure
neighb <- poly2nb(st_make_valid(MDR_NEW), queen = TRUE, snap = sqrt(0.001))

# Define the indices based on the spatial_IDS data
# From the file analysis:
# Index 54: TTYTKV_KIEN_HAI_KG
# Index 76: TTYTKV_PHU_QUOC_KG  
# Index 81: TTYTKV_RACH_GIA_KG

kien_hai_idx <- 54
phu_quoc_idx <- 76
rach_gia_idx <- 81

# Print current neighbor status
#cat("\nBefore adding neighbors:\n")
#cat("KIEN_HAI (", kien_hai_idx, ") neighbors:", neighb[[kien_hai_idx]], "\n")
#cat("PHU_QUOC (", phu_quoc_idx, ") neighbors:", neighb[[phu_quoc_idx]], "\n")
#cat("RACH_GIA (", rach_gia_idx, ") neighbors:", neighb[[rach_gia_idx]], "\n")

# Manually add neighbors
# For regions with no neighbors (represented as 0 or integer(0)), we need to replace them

# Add Rach Gia and Phu Quoc as neighbors to Kien Hai
if(length(neighb[[kien_hai_idx]]) == 0 || neighb[[kien_hai_idx]][1] == 0) {
  neighb[[kien_hai_idx]] <- as.integer(c(rach_gia_idx, phu_quoc_idx))
} else {
  neighb[[kien_hai_idx]] <- sort(unique(as.integer(c(neighb[[kien_hai_idx]], rach_gia_idx, phu_quoc_idx))))
}

# Add Rach Gia and Kien Hai as neighbors to Phu Quoc
if(length(neighb[[phu_quoc_idx]]) == 0 || neighb[[phu_quoc_idx]][1] == 0) {
  neighb[[phu_quoc_idx]] <- as.integer(c(rach_gia_idx, kien_hai_idx))
} else {
  neighb[[phu_quoc_idx]] <- sort(unique(as.integer(c(neighb[[phu_quoc_idx]], rach_gia_idx, kien_hai_idx))))
}

# Add Kien Hai and Phu Quoc as neighbors to Rach Gia (making it bidirectional)
neighb[[rach_gia_idx]] <- sort(unique(as.integer(c(neighb[[rach_gia_idx]], kien_hai_idx, phu_quoc_idx))))

# Clean the entire neighbor structure - remove any 0 values and ensure all are integers
neighb <- lapply(neighb, function(x) {
  x <- x[x != 0]
  as.integer(x)
})

# Update the class attributes
class(neighb) <- c("nb", "list")
attr(neighb, "region.id") <- as.character(1:length(neighb))
attr(neighb, "call") <- NULL

# Print updated neighbor status
#cat("\nAfter adding neighbors:\n")
#cat("KIEN_HAI (", kien_hai_idx, ") neighbors:", neighb[[kien_hai_idx]], "\n")
#cat("PHU_QUOC (", phu_quoc_idx, ") neighbors:", neighb[[phu_quoc_idx]], "\n")
#cat("RACH_GIA (", rach_gia_idx, ") neighbors:", neighb[[rach_gia_idx]], "\n")


# Create INLA graph
cat("\nCreating INLA graph...\n")
nb2INLA("MDR.graph", neighb)
MDR.adj <- paste(getwd(), "/MDR.graph", sep = "")



#date.test2 <- seq.Date(from=as.Date(max(d2$date)) %m-% months(5) ,to=as.Date(max(d2$date)) %m-% months(3) , by='month')
#date.test2 <- as.Date(max(d2$date%m-% months(3)))

date.test2 <- seq.Date(from=as.Date('2026-10-01')  ,to=as.Date('2026-05-01') , by='month')

length(date.test2)


all.fcodes <- unique(d2$fcode)

hyper.besag =   hyper = list(prec = list(prior = "loggamma",
                                          param = c(1, 1), initial = 0.01))

hyper1 = list(prec.unstruct=list(prior='pc.prec',param=c(3, 0.01)),
              prec.spatial=list(prior='pc.prec', param=c(3, 0.01)))


# iid model 
hyper.iid = list(theta = list(prior="pc.prec", param=c(1, 0.01)))

# ar1 model
hyper.ar1 = list(theta1 = list(prior='pc.prec', param=c(0.5, 0.01)),
                 rho = list(prior='pc.cor0', param = c(0.5, 0.75)))

hyper.ar2 <- list(
  theta1 = list(prior = "loggamma", param = c(3, 2)),
  theta2 = list(prior = "loggamma", param = c(3, 2))
)

# bym model
hyper.bym = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                 theta2 = list(prior="pc.prec", param=c(1, 0.01)))

# bym2 model
# probability of SD of theta1 > 1 = 0.01
hyper.bym2 = list(theta1 = list(prior="pc.prec", param=c(1, 0.01)),
                  theta2 = list(prior="pc", param=c(0.5, 0.5)))


# (puts more or less prior probability density on more or less wiggly)
hyper1.rw = list(prec = list(prior='pc.prec', param=c(0.1, 0.01))) # strictest smoothing; sd constrained to be low
hyper2.rw = list(prec = list(prior='pc.prec', param=c(0.3, 0.01))) # medium
hyper3.rw = list(prec = list(prior='pc.prec', param=c(1, 0.01))) # weaker (suggested INLA default) 
hyper4.rw = list(prec = list(prior='pc.prec', param=c(2, 0.01))) # weakest; sd can be quite wide 


