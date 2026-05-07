#' GEV Cumulative Distribution Function using Base R
#'
#' Calculates the CDF for the Generalized Extreme Value distribution using
#' highly optimized base R logical subsetting.
#'
#' @param x Numeric vector, matrix, or array. The data points.
#' @param mu Numeric vector, matrix, or array. Location parameter.
#' @param sigma Numeric vector, matrix, or array. Scale parameter.
#' @param xi Numeric vector, matrix, or array. Shape parameter.
#'
#' @export
gev_cdf <- function(x, mu, sigma, xi) {
  z <- (x - mu) / sigma

  gumbel <- exp(-exp(-z))

  suppressWarnings({
    gev <- exp(-((1 + xi * z)^(-1 / xi)))
  })

  bound_mask <- (1 + xi * z) <= 0
  gev[bound_mask & xi > 0] <- 0
  gev[bound_mask & xi <= 0] <- 1

  cdf <- gumbel
  gev_mask <- abs(xi) > 1e-10
  cdf[gev_mask] <- gev[gev_mask]

  return(cdf)
}

#' Probability-Probability (P-P) Plot for d-GEV
#'
#' Creates a P-P plot comparing empirical probabilities to theoretical
#' d-GEV probabilities with credible intervals using base R.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param j Integer. Index of the station to plot.
#' @param alpha Numeric. Significance level for credible intervals (default 0.05).
#'
#' @importFrom graphics plot points lines legend abline
#' @importFrom stats quantile
#' @export
pp_plot <- function(fit, j, alpha = 0.05) {

  data_j <- t(apply(fit$data[[j]], 2, sort, na.last = TRUE))
  D <- nrow(data_j)
  L <- ncol(data_j)

  durs <- fit$durs
  n_durs <- length(durs)

  upper_q <- 1 - (alpha / 2)
  lower_q <- alpha / 2

  p_i <- (1:L) / (L + 1)

  mut    <- fit$pars$mut[, j]
  sigma0 <- fit$pars$sigma0[, j]
  theta  <- fit$pars$theta[, j]
  eta    <- fit$pars$eta[, j]

  N <- length(mut)

  sigma_mat <- matrix(0, nrow = N, ncol = D)
  for (dd in 1:D) {
    sigma_mat[, dd] <- sigma0 * (durs[dd] + theta)^(-eta)
  }
  mu_mat <- mut * sigma_mat

  if (fit$shp_d == "d") {
    xi_mat <- fit$pars$xi
  } else {
    xi_mat <- matrix(fit$pars$xi[, j], nrow = N, ncol = D)
  }

  x_array     <- array(rep(data_j, each = N), dim = c(N, D, L))
  mu_array    <- array(rep(mu_mat, times = L), dim = c(N, D, L))
  sigma_array <- array(rep(sigma_mat, times = L), dim = c(N, D, L))
  xi_array    <- array(rep(xi_mat, times = L), dim = c(N, D, L))

  cdf_t <- gev_cdf(x_array, mu_array, sigma_array, xi_array)

  cdf_mean <- apply(cdf_t, MARGIN = c(2, 3), mean, na.rm = TRUE)

  dim(cdf_t) <- c(N * D, L)
  upper_ci <- apply(cdf_t, MARGIN = 2, quantile, probs = upper_q, na.rm = TRUE)
  lower_ci <- apply(cdf_t, MARGIN = 2, quantile, probs = lower_q, na.rm = TRUE)

  plot_cols <- (1:n_durs) + 1
  plot_pchs <- 1:n_durs

  plot(0, xlim = c(0, 1), ylim = c(0, 1), type = 'n',
       main = fit$stationID[j], xlab = 'empirical probability', ylab = 'theoretical probability')

  for (dd in 1:n_durs) {
    points(p_i, cdf_mean[dd, ], col = plot_cols[dd], pch = plot_pchs[dd])
  }

  lines(p_i, upper_ci, lty = 2, col = '#80808080')
  lines(p_i, lower_ci, lty = 2, col = '#80808080')

  legend('topleft',
         legend = paste0(durs, "h"),
         col = plot_cols, pch = plot_pchs,
         ncol = 2, bty = 'n', cex = 0.8,
         pt.cex = 1.2, text.width = 0.1)

  abline(0, 1)
}
