#' Posterior Predictive Histogram Plot (Base R)
#'
#' Generates a histogram of the posterior predictive distribution for a specific
#' statistic and compares it against the observed (empirical) value.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param stat Numeric (0-1) or Character ("median", "max", "min"). The statistic to check.
#' @param d Integer. Index of the duration to plot.
#' @param j Integer. Index of the station to plot.
#'
#' @importFrom graphics hist abline
#' @export
post_pred_plot <- function(fit, stat, d, j) {

  if (is.numeric(stat) && (stat <= 0 || stat >= 1)) {
    stop('please enter the quantile between 0 and 1.0')
  }
  if (!is.numeric(stat) && !any(stat == c('median', 'max', 'min'))) {
    stop('please enter a stat: "max","min","median"')
  }

  dur  <- fit$durs[d]
  data <- fit$data[[j]][, d]
  N    <- length(data)

  mut    <- fit$pars$mut[, j]
  sigma0 <- fit$pars$sigma0[, j]
  theta  <- fit$pars$theta[, j]
  eta    <- fit$pars$eta[, j]

  sigma <- sigma0 * (dur + theta)^(-eta)
  mu    <- sigma * mut

  if (fit$shp_d == 'd') {
    xi <- fit$pars$xi[, d]
  } else if (fit$shp_d == 'j') {
    xi <- fit$pars$xi[, j]
  }

  qgev_scalar_p <- function(p_val) {
    eps <- 1e-10
    log_p <- -log(p_val)
    multiplier <- ifelse(abs(xi) > eps,
                         (1 - log_p^(-xi)) / xi,
                         log(log_p))
    return(mu - sigma * multiplier)
  }

  if (is.numeric(stat)) {
    q_pp  <- qgev_scalar_p(stat)
    q_emp <- quantile(data, probs = stat, names = FALSE, na.rm = TRUE)

  } else if (stat == "median") {
    # Median is the 0.5 quantile
    q_pp  <- qgev_scalar_p(0.5)
    q_emp <- median(data, na.rm = TRUE)

  } else if (stat == "max") {
    p <- 1 - (1 / N)
    q_pp  <- qgev_scalar_p(p)
    q_emp <- max(data, na.rm = TRUE)

  } else if (stat == "min") {
    p <- 1 / N
    q_pp  <- qgev_scalar_p(p)
    q_emp <- min(data, na.rm = TRUE)
  }

  x_max <- max(q_pp, q_emp)
  x_min <- min(q_pp, q_emp)
  x_min_plot <- max(0, x_min - (x_max - x_min) * 0.05)

  stat_title <- if (is.numeric(stat)) paste0("Quantile = ", stat) else tools::toTitleCase(stat)

  hist(q_pp, xlim = c(x_min_plot, x_max * 1.05),
       main = paste0(stat_title, ', xi = ', fit$shp_d, ', Dur = ', dur, 'h\n', fit$stationID[j]),
       cex.main = 0.85, xlab = 'Intensity [mm/h]')
  abline(v = q_emp, col = 'red', lty = 2, lwd = 2)
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
#'
#' @return Numeric. The calculated PPP.
#' @export
post_pred_pval <- function(fit, stat, d, j) {

  dur  <- fit$durs[d]
  data <- fit$data[[j]][, d]
  N    <- length(data)

  mut    <- fit$pars$mut[, j]
  sigma0 <- fit$pars$sigma0[, j]
  theta  <- fit$pars$theta[, j]
  eta    <- fit$pars$eta[, j]

  sigma <- sigma0 * (dur + theta)^(-eta)
  mu    <- sigma * mut

  if (fit$shp_d == 'd') {
    xi <- fit$pars$xi[, d]
  } else if (fit$shp_d == 'j') {
    xi <- fit$pars$xi[, j]
  }

  qgev_scalar_p <- function(p_val) {
    eps <- 1e-10
    log_p <- -log(p_val)
    multiplier <- ifelse(abs(xi) > eps,
                         (1 - log_p^(-xi)) / xi,
                         log(log_p))
    return(mu - sigma * multiplier)
  }

  if (is.numeric(stat)) {
    q_pp  <- qgev_scalar_p(stat)
    q_emp <- quantile(data, probs = stat, names = FALSE, na.rm = TRUE)

  } else if (stat == "median") {
    # Median is the 0.5 quantile
    q_pp  <- qgev_scalar_p(0.5)
    q_emp <- median(data, na.rm = TRUE)

  } else if (stat == "max") {
    p <- 1 - (1 / N)
    q_pp  <- qgev_scalar_p(p)
    q_emp <- max(data, na.rm = TRUE)

  } else if (stat == "min") {
    p <- 1 / N
    q_pp  <- qgev_scalar_p(p)
    q_emp <- min(data, na.rm = TRUE)
  }

  ppp_res <- mean(q_pp > q_emp)

  return(ppp_res)
}
