#' GEV Cumulative Distribution Function using Torch
#'
#' Calculates the CDF for the Generalized Extreme Value distribution using
#' tensor operations.
#'
#' @param x Tensor. The data points.
#' @param mu Tensor. Location parameter.
#' @param sigma Tensor. Scale parameter.
#' @param xi Tensor. Shape parameter.
#'
#' @importFrom torch torch_tensor torch_exp torch_pow torch_where torch_abs
#' @export
gev_cdf_torch = function(x, mu, sigma, xi) {
  z = (x - mu) / sigma
  one = torch::torch_tensor(1.0, dtype = torch::torch_float64())
  gumbel = torch::torch_exp(-torch::torch_exp(-z))
  gev    = torch::torch_exp(-torch::torch_pow(one + xi * z, -1/xi))
  eps = torch::torch_tensor(1e-10, dtype = torch::torch_float64())
  cdf = torch::torch_where(torch::torch_abs(xi) > eps, gev, gumbel)
  return(cdf)
}

#' Probability-Probability (P-P) Plot for d-GEV
#'
#' Creates a P-P plot comparing empirical probabilities to theoretical
#' d-GEV probabilities with credible intervals.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param j Integer. Index of the station to plot.
#' @param alpha Numeric. Significance level for credible intervals (default 0.05).
#'
#' @importFrom torch torch_tensor torch_pow torch_quantile torch_mean
#' @importFrom graphics plot points lines legend abline
#' @export
PP_plot = function(fit, j, alpha = 0.05) {
  data_j = t(apply(fit$data[[j]], 2, sort))
  data_j_t = torch::torch_tensor(data_j, dtype = torch::torch_float64())
  data_length = dim(data_j)[2]
  durs = fit$durs
  durs_t = torch::torch_tensor(durs, dtype = torch::torch_float64())

  upper_q_t = torch::torch_tensor(1 - (alpha/2), dtype = torch::torch_float64())
  lower_q_t = torch::torch_tensor(alpha/2, dtype = torch::torch_float64())

  p_i = (1:data_length) / (data_length + 1)

  par_names = c('mut', 'sigma0', 'theta', 'eta')
  pars_t = lapply(par_names, function(p_name) {
    torch::torch_tensor(fit$pars[[p_name]][,j], dtype = torch::torch_float64())$view(c(-1, 1, 1))
  })

  durs_t2 = durs_t$view(c(1, -1, 1))
  sigma = pars_t[[2]] * torch::torch_pow(durs_t2 + pars_t[[3]], -pars_t[[4]])
  mu    = pars_t[[1]] * sigma

  if (fit$shp_d == "d") {
    xi = torch::torch_tensor(fit$pars$xi, dtype = torch::torch_float64())$unsqueeze(3)
  } else {
    xi = torch::torch_tensor(fit$pars$xi[,j], dtype = torch::torch_float64())$view(c(-1, 1, 1))
  }

  cdf_t = gev_cdf_torch(data_j_t$unsqueeze(1), mu, sigma, xi)
  cdf_mean_t = torch::torch_mean(cdf_t, dim = 1)
  shape_dim = c(prod(cdf_t$shape[1:2]), cdf_t$shape[3])

  upper_ci = torch::torch_quantile(cdf_t$view(shape_dim), q = upper_q_t, dim = 1)
  lower_ci = torch::torch_quantile(cdf_t$view(shape_dim), q = lower_q_t, dim = 1)

  n_durs = length(fit$durs)
  plot_cols = (1:n_durs) + 1
  plot_pchs = 1:n_durs

  plot(0, xlim = c(0, 1), ylim = c(0, 1), type = 'n',
       main = fit$stationID[j], xlab = 'empirical probability', ylab = 'theoretical probability')
  for (dd in 1:length(durs)) {
    points(p_i, as.numeric(cdf_mean_t[dd,]),
           col = plot_cols[dd], pch = plot_pchs[dd])
  }
  lines(as.numeric(upper_ci), p_i, lty = 2, col = '#80808080')
  lines(as.numeric(lower_ci), p_i, lty = 2, col = '#80808080')
  legend('topleft',
         legend = paste0(fit$durs, "h"),
         col = plot_cols, pch = plot_pchs,
         ncol = 2,bty = 'n', cex = 0.8,
         pt.cex = 1.2, text.width = 0.1)
  abline(0, 1)
}
