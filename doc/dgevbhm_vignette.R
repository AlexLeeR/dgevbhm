## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(dgevbhm)

## ----data_path, eval=FALSE----------------------------------------------------
# # Replace "/path/to/directory/" with the actual path where you saved the data
# UK_data_full <- data_read(dir_path = "/path/to/directory/",
#                           missing_data = 5,
#                           durs = c(1, 3, 6, 12, 24, 48),
#                           cores = 4)

## ----full_data_length---------------------------------------------------------
length(UK_data_full$y)

## ----mean_median--------------------------------------------------------------
# mean and median hourly precipitation rate (mm)
mean(UK_data_full$amp)
median(UK_data_full$amp)

## ----selection----------------------------------------------------------------
med_amp <- median(UK_data_full$amp)
# similar stations: selection of 25 stations that are closest to the median precipitation rate
sim_stations <- which(rank(abs(UK_data_full$amp - med_amp)) <= 25)
# plotting a histogram of the 25 precipitation rates
hist(UK_data_full$amp[sim_stations], main = 'histogram', xlab = 'stations')

## ----data_names---------------------------------------------------------------
names(UK_data_full) #checking header of the GSDR data

## ----latlondim----------------------------------------------------------------
dim(UK_data_full$latlon) #dimension of the latlon 

## ----clean_data---------------------------------------------------------------
# clean the data: select the 25 rain gauges and store them in a new object 
UK_data <- list()
UK_data$y <- UK_data_full$y[sim_stations]
UK_data$amp <- UK_data_full$amp[sim_stations]
UK_data$latlon <- UK_data_full$latlon[sim_stations, ]
UK_data$stationID <- UK_data_full$stationID[sim_stations]
# don't forget the durations!
UK_data$durs <- c(1, 3, 6, 12, 24, 48)

## ----dgev_bhm_run_d-----------------------------------------------------------
# Fit the dGEV-BHM with the shape parameter fixed over durations (shp_d = 'd'). We opt for 2000 iterations of the MCMC sampler, where the first half is used as burn-in. 
# We assume that 2 cores are available for computation. 
# As you can only use one core per chain, we use 2 chains.

# WARNING: May take some time. You can also load the prepared fitted models.
fit_xi_d <- dgev_bhm(data = UK_data, chains = 2, iter = 2000, cores = 2, shp_d = 'd')

## ----dgev_bhm_run_j-----------------------------------------------------------
# Fit the dGEV-BHM with the shape parameter fixed over space (shp_d = 'j').

# WARNING: May take some time. You can also load the prepared fitted models.
fit_xi_j <- dgev_bhm(data = UK_data, chains = 2, iter = 2000, cores = 2, shp_d = 'j')

## ----fit_names----------------------------------------------------------------
# structure of the fitted models
names(fit_xi_d) 
names(fit_xi_j)

## ----mcmc_diag----------------------------------------------------------------
# MCMC diagnostics for mut
fit_xi_d$mcmc_diagnostics$rhat$mut
fit_xi_d$mcmc_diagnostics$ess$mut

## ----ic-----------------------------------------------------------------------
fit_xi_j$information_criteria
fit_xi_d$information_criteria

## ----pp_plot_d, fig.width=7, fig.height=5, out.width="100%"-------------------
pp_plot(fit = fit_xi_d, j = 5)
pp_plot(fit = fit_xi_j, j = 5)

## ----post_pred_plot_median, fig.width=7, fig.height=8, out.width="100%"-------
par(mfrow = c(2, 1))
post_pred_plot(fit = fit_xi_d, stat = 'median', j = 1, d = 6) # posterior predictive diagnostic plot for rain gauge #1 and duration 48h for the fit_xi_d model 
post_pred_plot(fit = fit_xi_j, stat = 'median', j = 1, d = 6) # posterior predictive diagnostic plot for rain gauge #1 and duration 48h for the fit_xi_j model

## ----post_pred_plot_min, fig.width=7, fig.height=8, out.width="100%"----------
par(mfrow = c(2, 1))
post_pred_plot(fit = fit_xi_d, stat = 'min', j = 6, d = 1) # posterior predictive diagnostic plot for rain gauge #6 and duration 1h for the fit_xi_d model 
post_pred_plot(fit = fit_xi_j, stat = 'min', j = 6, d = 1) # posterior predictive diagnostic plot for rain gauge #6 and duration 1h for the fit_xi_j model

## ----post_pred_plot_max, fig.width=7, fig.height=8, out.width="100%"----------
par(mfrow = c(2, 1))
post_pred_plot(fit = fit_xi_d, stat = 'max', j = 3, d = 1)
post_pred_plot(fit = fit_xi_j, stat = 'max', j = 3, d = 1)

## ----post_pred_plot_pval, fig.width=7, fig.height=8, out.width="100%"---------
post_pred_pval(fit = fit_xi_d, stat = 'median', j = 2, d = 6)
post_pred_pval(fit = fit_xi_j, stat = 'median', j = 2, d = 6)

## ----return_level, fig.width=7, fig.height=5, out.width="100%"----------------
# let's try a significance level of 0.1
return_level_plot(fit_xi_d, j = 5, d = 1, max_rp = 100, alpha = 0.1)
return_level_plot(fit_xi_j, j = 5, d = 1, max_rp = 100, alpha = 0.1)

## ----IDF1, fig.width=7, fig.height=5, out.width="100%"------------------------
idf_plot(fit = fit_xi_d,j = 5, rp = 100, alpha = 0.05)

## ----IDF3, fig.width=7, fig.height=5, out.width="100%"------------------------
idf_plot(fit = fit_xi_d,j = 5, rp = c(2,10,100), alpha = 0.05)

