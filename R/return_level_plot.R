#' Plot Return Level Curves for d-GEV (Base R)
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
#' @importFrom graphics plot polygon lines points
#' @importFrom matrixStats colQuantiles
#' @export
return_level_plot <- function(fit, max_rp = 100, j, d, alpha = 0.05) {

  if (is.numeric(max_rp) && (max_rp > 1000 || max_rp < 5)) {
    stop('Return period cannot be larger than 1000 or smaller than 5')
  }

  dur  <- fit$durs[d]
  data <- fit$data[[j]][, d]
  N    <- length(data)

  # Correct Weibull plotting positions
  # For sorted data (largest to smallest): i=1 is largest, i=N is smallest
  data_sorted <- sort(data, decreasing = TRUE)
  emp_T <- (N + 1) / (1:N)  # Return periods for sorted data

  if (max_rp == 'data') {
    max_rp <- max(emp_T)  # Use the maximum empirical return period
  }

  # Generate evenly-spaced return periods on log scale
  T_periods <- exp(seq(log(min(emp_T)), log(max_rp), length.out = 100))
  prob_vec  <- 1 - (1 / T_periods)

  mut    <- fit$pars$mut[, j]
  sigma0 <- fit$pars$sigma0[, j]
  theta  <- fit$pars$theta[, j]
  eta    <- fit$pars$eta[, j]

  sigma <- sigma0 * (dur + theta)^(-eta)
  mu    <- mut * sigma

  if (fit$shp_d == 'd') {
    xi <- fit$pars$xi[, d]
  } else if (fit$shp_d == 'j') {
    xi <- fit$pars$xi[, j]
  }

  rl_mat <- sapply(prob_vec, function(p) {
    eps <- 1e-10
    log_p <- -log(p)
    multiplier <- ifelse(abs(xi) > eps,
                         (1 - log_p^(-xi)) / xi,
                         log(log_p))
    return(mu - sigma * multiplier)
  })

  rl_mean  <- colMeans(rl_mat)
  rl_upper <- matrixStats::colQuantiles(rl_mat, probs = 1 - alpha/2)
  rl_lower <- matrixStats::colQuantiles(rl_mat, probs = alpha/2)

  # Filter empirical points to plot only those within the display range
  plot_idx <- emp_T <= max_rp

  plot(T_periods, rl_mean, type = 'n', log = 'xy',
       ylim = range(c(rl_lower, rl_upper, data_sorted[plot_idx])),
       xlab = 'Return period [y]', ylab = 'Return level [mm/h]',
       main = paste0('StationID: ', fit$stationID[j], '\nduration = ', dur, 'h'))

  # Add credible interval
  polygon(c(T_periods, rev(T_periods)),
          c(rl_upper, rev(rl_lower)),
          col = '#FFFF00', border = NA)  # Light gray with transparency

  # Add mean return level curve
  lines(T_periods, rl_mean, col = '#D70064', lwd = 2)

  # Add empirical points
  points(emp_T[plot_idx], data_sorted[plot_idx], col = '#0050A0')

  # Add legend
  legend('bottomright',
         legend = c('Posterior mean', paste0((1-alpha)*100, '% Credible interval'), 'Empirical'),
         col = c('#D70064', '#FFFF00', '#0050A0'),
         lty = c(1, NA, NA),
         lwd = c(2, NA, NA),
         pch = c(NA, 15, 1),
         pt.cex = c(NA, 2, 1),
         bty = 'n')
}
