#' Posterior Predictive Histogram Plot
#'
#' Generates a histogram of the posterior predictive distribution for a specific
#' statistic and compares it against the observed (empirical) value.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param stat Numeric (0-1) or Character ("median", "max", "min"). The statistic to check.
#' @param d Integer. Index of the duration to plot.
#' @param j Integer. Index of the station to plot.
#' @param cores Integer. Number of threads for torch.
#'
#' @importFrom torch torch_tensor torch_pow torch_quantile torch_rand torch_sort torch_median torch_max torch_min torch_set_num_threads
#' @importFrom graphics hist abline
#' @export
pp_plot_torch = function(fit, stat, d, j, cores = 4) {
  # torch is handled via Imports and torch:: prefix

  if (is.numeric(stat) && (stat >= 1 || stat <= 0)) {
    stop('please enter the quantile between 0 and 1.0')
  }
  if (!is.numeric(stat) && !any(stat == c('median', 'max', 'min'))) {
    stop('please enter a stat: "max","min","median"')
  }

  torch::torch_set_num_threads(cores)

  dur    = fit$durs
  dur_t  = torch::torch_tensor(dur, dtype = torch::torch_float64())
  data   = fit$data[[j]][,d]
  data_t = torch::torch_tensor(data, dtype = torch::torch_float64())

  par_names = c('mut', 'sigma0', 'theta', 'eta')

  pars_t = lapply(par_names, function(p_name) {
    pars_data = fit$pars[[p_name]][,j]
    return(torch::torch_tensor(pars_data, dtype = torch::torch_float64()))
  })

  sigma_t = pars_t[[2]] * torch::torch_pow(dur_t[d] + pars_t[[3]], -pars_t[[4]])
  mu_t    = sigma_t * pars_t[[1]]

  # define shape parameter logic
  if (fit$shp_d == 'd') {
    xi_t = torch::torch_tensor(fit$pars$xi[,d], dtype = torch::torch_float64())
  } else if (fit$shp_d == 'j') {
    xi_t = torch::torch_tensor(fit$pars$xi[,j], dtype = torch::torch_float64())
  }

  if (is.numeric(stat)) {
    stat_t  = torch::torch_tensor(stat, dtype = torch::torch_float64())
    q_pp    = dgev_torch(p = stat_t, mu = mu_t, sigma = sigma_t, xi = xi_t)
    q_emp_t = torch::torch_quantile(data_t, q = stat_t)

    val_q_pp = as.numeric(q_pp)
    val_q_emp = as.numeric(q_emp_t)

    x_max = max(val_q_pp, val_q_emp)
    x_min = min(val_q_pp, val_q_emp)
    x_min_plot = max(0, x_min - (x_max - x_min) * 0.05)

    hist(val_q_pp,
         xlim = c(x_min_plot, x_max * 1.05),
         main = paste0('Quantile = ', stat, ', xi = ', fit$shp_d, ', Dur = ', fit$durs[d], 'h',
                       '\n', fit$stationID[j]),
         cex.main = 0.85, xlab = 'Intensity [mm/h]')
    abline(v = val_q_emp, col = 'red', lty = 2, lwd = 2)

  } else if (stat == "median") {
    uni_rand     = torch::torch_rand(size = c(10001), dtype = torch::torch_float64())
    q_pp         = dgev_torch(p = uni_rand$unsqueeze(1),
                              mu = mu_t$unsqueeze(2),
                              sigma = sigma_t$unsqueeze(2),
                              xi = xi_t$unsqueeze(2))
    q_pp_median  = torch::torch_sort(q_pp, dim = 2)[[1]][.., 5000]
    q_emp_median = torch::torch_median(data_t)

    val_med_pp = as.numeric(q_pp_median)
    val_med_emp = as.numeric(q_emp_median)

    x_max = max(val_med_pp, val_med_emp)
    x_min = min(val_med_pp, val_med_emp)
    x_min_plot = max(0, x_min - (x_max - x_min) * 0.05)

    hist(val_med_pp,
         xlim = c(x_min_plot, x_max * 1.05),
         main = paste0('Median, xi = ', fit$shp_d, ', Dur = ', fit$durs[d], 'h',
                       '\n', fit$stationID[j]),
         cex.main = 0.85, xlab = 'Intensity [mm/h]')
    abline(v = val_med_emp, col = 'red', lty = 2, lwd = 2)

  } else if (stat == "max") {
    data_len = data_t$shape[1]
    p        = 1 - (1 / data_len)
    q_pp     = dgev_torch(p = torch::torch_tensor(p, dtype = torch::torch_float64()),
                          mu = mu_t, sigma = sigma_t, xi = xi_t)
    q_emp    = torch::torch_max(data_t)

    val_q_pp = as.numeric(q_pp)
    val_q_emp = as.numeric(q_emp)

    x_max = max(val_q_pp, val_q_emp)
    x_min = min(val_q_pp, val_q_emp)
    x_min_plot = max(0, x_min - (x_max - x_min) * 0.05)

    hist(val_q_pp, xlim = c(x_min_plot, x_max * 1.05),
         main = paste0('Max, xi = ', fit$shp_d, ', Dur = ', fit$durs[d], 'h',
                       '\n', fit$stationID[j]),
         cex.main = 0.85, xlab = 'Intensity [mm/h]')
    abline(v = val_q_emp, col = 'red', lty = 2, lwd = 2)

  } else if (stat == "min") {
    data_len = data_t$shape[1]
    p        = 1 / data_len
    q_pp     = dgev_torch(p = torch::torch_tensor(p, dtype = torch::torch_float64()),
                          mu = mu_t, sigma = sigma_t, xi = xi_t)
    q_emp    = torch::torch_min(data_t)

    val_q_pp = as.numeric(q_pp)
    val_q_emp = as.numeric(q_emp)

    x_max = max(val_q_pp, val_q_emp)
    x_min = min(val_q_pp, val_q_emp)
    x_min_plot = max(0, x_min - (x_max - x_min) * 0.05)

    hist(val_q_pp, xlim = c(x_min_plot, x_max * 1.05),
         main = paste0('Min, xi = ', fit$shp_d, ', Dur = ', fit$durs[d], 'h',
                       '\n', fit$stationID[j]),
         cex.main = 0.85, xlab = 'Intensity [mm/h]')
    abline(v = val_q_emp, col = 'red', lty = 2, lwd = 2)
  }
}

#' Posterior Predictive P-value (PPP)
#'
#' Calculates the posterior predictive p-value for a given statistic.
#' A value near 0.5 indicates a good model fit.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param stat Numeric (0-1) or Character ("median", "max", "min"). The statistic to check.
#' @param d Integer. Index of the duration.
#' @param j Integer. Index of the station.
#' @param cores Integer. Number of threads for torch.
#'
#' @return Numeric. The calculated PPP value.
#' @export
ppp_torch = function(fit, stat, d, j, cores = 4) {
  torch::torch_set_num_threads(cores)

  dur_t  = torch::torch_tensor(fit$durs[d], dtype = torch::torch_float64())
  data_t = torch::torch_tensor(fit$data[[j]][,d], dtype = torch::torch_float64())

  par_names = c('mut', 'sigma0', 'theta', 'eta')
  pars_t = lapply(par_names, function(p) {
    torch::torch_tensor(fit$pars[[p]][,j], dtype = torch::torch_float64())
  })

  sigma_t = pars_t[[2]] * torch::torch_pow(dur_t + pars_t[[3]], -pars_t[[4]])
  mu_t    = sigma_t * pars_t[[1]]

  xi_t = if(fit$shp_d == 'd') torch::torch_tensor(fit$pars$xi[,d], dtype = torch::torch_float64()) else
    torch::torch_tensor(fit$pars$xi[,j], dtype = torch::torch_float64())

  if (is.numeric(stat)) {
    stat_t  = torch::torch_tensor(stat, dtype = torch::torch_float64())
    q_pp    = dgev_torch(p = stat_t, mu = mu_t, sigma = sigma_t, xi = xi_t)
    q_emp_t = torch::torch_quantile(data_t, q = stat_t)
    ppp_res = torch::torch_where(q_pp > q_emp_t, 1, 0)
  } else if (stat == "median") {
    uni_rand     = torch::torch_rand(size = c(10001), dtype = torch::torch_float64())
    q_pp         = dgev_torch(p = uni_rand$unsqueeze(1),
                              mu = mu_t$unsqueeze(2),
                              sigma = sigma_t$unsqueeze(2),
                              xi = xi_t$unsqueeze(2))
    q_pp_median  = torch::torch_sort(q_pp, dim = 2)[[1]][.., 5000]
    q_emp_median = torch::torch_median(data_t)
    ppp_res = torch::torch_where(q_pp_median > q_emp_median, 1, 0)
  } else if (stat == "max") {
    p = 1 - (1 / data_t$shape[1])
    q_pp = dgev_torch(p = torch::torch_tensor(p, dtype = torch::torch_float64()), mu = mu_t, sigma = sigma_t, xi = xi_t)
    ppp_res = torch::torch_where(q_pp > torch::torch_max(data_t), 1, 0)
  } else if (stat == "min") {
    p = 1 / data_t$shape[1]
    q_pp = dgev_torch(p = torch::torch_tensor(p, dtype = torch::torch_float64()), mu = mu_t, sigma = sigma_t, xi = xi_t)
    ppp_res = torch::torch_where(q_pp > torch::torch_min(data_t), 1, 0)
  }

  return(as.numeric(torch::torch_mean(ppp_res$to(dtype = torch::torch_float64()))))
}
