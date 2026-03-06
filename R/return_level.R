#' Plot Return Level Curves for d-GEV
#'
#' Calculates and plots the relationship between return periods and precipitation
#' intensity (return levels) with Bayesian credible intervals.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param max_rp Numeric or Character. Maximum return period to plot. Use "data"
#' to cap at the record length.
#' @param j Integer. Index of the station.
#' @param d Integer. Index of the duration.
#' @param alpha Numeric. Significance level for the credible interval (default 0.05).
#'
#' @importFrom torch torch_tensor torch_pow torch_mean torch_quantile
#' @importFrom graphics plot polygon lines points
#' @export
return_level = function(fit, max_rp = 100, j, d, alpha = 0.05) {
  if (is.numeric(max_rp) && (max_rp > 1000 || max_rp < 5)) {
    stop('Return period cannot be larger 1000 or smaller than 5')
  }

  durs_t = torch::torch_tensor(fit$durs, dtype = torch::torch_float64())
  data = fit$data[[j]][,d]
  N = length(data)
  emp_T = (N + 1) / (1:N)

  if (max_rp == 'data') max_rp = emp_T[1]

  data_sorted = sort(data, decreasing = TRUE)
  T_periods = exp(seq(log(emp_T[N]), log(max_rp), length.out = 100))
  prob_t = torch::torch_tensor(1 - (1 / T_periods), dtype = torch::torch_float64())

  par_names = c('mut', 'sigma0', 'theta', 'eta')
  pars_t = lapply(par_names, function(p) {
    torch::torch_tensor(fit$pars[[p]][,j], dtype = torch::torch_float64())
  })

  sigma = pars_t[[2]] * torch::torch_pow(durs_t[d] + pars_t[[3]], -pars_t[[4]])
  mu    = pars_t[[1]] * sigma
  xi = if(fit$shp_d == 'd') torch::torch_tensor(fit$pars$xi[,d], dtype=torch::torch_float64()) else
    torch::torch_tensor(fit$pars$xi[,j], dtype=torch::torch_float64())

  rl_all = dgev_torch(prob_t$unsqueeze(1), mu$unsqueeze(2), sigma$unsqueeze(2), xi$unsqueeze(2))

  rl_mean  = torch::torch_mean(rl_all, dim = 1)
  rl_upper = torch::torch_quantile(rl_all, q = torch::torch_tensor(1-alpha/2, dtype = torch::torch_float64()), dim = 1)$squeeze()
  rl_lower = torch::torch_quantile(rl_all, q = torch::torch_tensor(alpha/2, dtype = torch::torch_float64()), dim = 1)$squeeze()

  plot(T_periods, as.numeric(rl_mean), type = 'n', log = 'xy',
       ylim = range(c(as.numeric(rl_lower), as.numeric(rl_upper))),
       xlab = 'Return period [y]', ylab = 'Return level [mm/h]',
       main = paste0('StationID: ',fit_xi_d$stationID[j],'\n','duration = ', fit$durs[d], 'h'))
  polygon(c(T_periods, rev(T_periods)), c(as.numeric(rl_upper), rev(as.numeric(rl_lower))), col = '#FFFF00', border = NA)
  lines(T_periods, as.numeric(rl_mean), col = '#D70064', lwd = 2)
  points(emp_T[emp_T <= max_rp], data_sorted[emp_T <= max_rp], col = '#0050A0')
}
