#' Quantile function for d-GEV using Base R
#'
#' @param p Probability (scalar)
#' @param mu_mat Location parameter matrix (Durations x Samples)
#' @param sigma_mat Scale parameter matrix (Durations x Samples)
#' @param xi_vec Shape parameter vector
#' @param xi_dim Character. "duration" for shp_d='d', "sample" for shp_d='j'
#' @export
gev_mat <- function(p, mu_mat, sigma_mat, xi_vec, xi_dim) {
  eps <- 1e-10
  log_p <- -log(p)

  multiplier <- ifelse(abs(xi_vec) > eps,
                       (1 - log_p^(-xi_vec)) / xi_vec,
                       log(log_p))

  if (xi_dim == "duration") {
    sigma_scaled <- sigma_mat * multiplier
  } else if (xi_dim == "sample") {
    sigma_scaled <- sweep(sigma_mat, 2, multiplier, `*`)
  }

  return(mu_mat - sigma_scaled)
}

#' Plot IDF Curves using Base R
#'
#' @param fit List. The output object from the dgev_bhm function.
#' @param j Integer. Station index.
#' @param rp Numeric vector. Return periods (max 3).
#' @param alpha Numeric. Significance level.
#'
#' @importFrom graphics plot lines polygon legend
#' @importFrom matrixStats rowQuantiles
#' @importFrom grDevices adjustcolor
#' @export
idf_plot <- function(fit, j = NULL, rp, alpha = 0.05) {

  if (length(rp) > 3) stop('No more than 3 return periods allowed!')
  if (alpha < 0 | alpha > 1) stop("alpha should be less than 1 and greater than 0")

  durs <- fit$durs
  if (length(durs) <= 3) stop("more than 3 durations needed")

  prob_vec <- 1 - (1 / rp)
  si_col <- c("#D81B60", "#1E88E5", "#004D40")
  quant_probs <- c(0 + (alpha / 2), 1 - (alpha / 2))

  interpolate_vector <- function(v, steps = 101) {
    res <- numeric(0)
    for (i in 2:length(v)) {
      seq_val <- seq(v[i - 1], v[i], length.out = steps)
      if (i == 2) {
        res <- c(res, seq_val)
      } else {
        res <- c(res, seq_val[-1])
      }
    }
    return(res)
  }

  durs_ip <- interpolate_vector(durs, steps = 101)

  dur_indices <- seq(1, length(durs_ip), by = 100)

  calc_mu_sigma <- function(durs_vec, mut, sigma0, theta, eta) {
    term_inside <- outer(durs_vec, theta, `+`)
    term_pow    <- sweep(term_inside, 2, -eta, `^`)
    sigma_mat   <- sweep(term_pow, 2, sigma0, `*`)
    mu_mat      <- sweep(sigma_mat, 2, mut, `*`)
    return(list(mu = mu_mat, sigma = sigma_mat))
  }

  mut    <- fit$pars$mut[, j]
  sigma0 <- fit$pars$sigma0[, j]
  theta  <- fit$pars$theta[, j]
  eta    <- fit$pars$eta[, j]

  # plotting objects
  res_list <- list()
  all_lower <- numeric(0)
  all_upper <- numeric(0)

  if (fit$shp_d == "d") {
    cat("Computing d-dimensional shape IDF curve...\n")

    xi_mean    <- colMeans(fit$pars$xi)
    xi_mean_ip <- interpolate_vector(xi_mean, steps = 101)

    params_dv     <- calc_mu_sigma(durs_ip, mut, sigma0, theta, eta)
    params_actual <- calc_mu_sigma(durs, mut, sigma0, theta, eta)

    for (si in 1:length(prob_vec)) {
      p <- prob_vec[si]

      idf_samp_dv <- gev_mat(p, params_dv$mu, params_dv$sigma, xi_mean_ip, "duration")
      idf_mean_dv <- rowMeans(idf_samp_dv)
      quants_dv   <- matrixStats::rowQuantiles(idf_samp_dv, probs = quant_probs)

      idf_samp_actual <- gev_mat(p, params_actual$mu, params_actual$sigma, xi_mean, "duration")
      idf_mean_actual <- rowMeans(idf_samp_actual)
      quants_actual   <- matrixStats::rowQuantiles(idf_samp_actual, probs = quant_probs)

      diff_mean  <- idf_mean_actual - idf_mean_dv[dur_indices]
      diff_lower <- quants_actual[, 1] - quants_dv[dur_indices, 1]
      diff_upper <- quants_actual[, 2] - quants_dv[dur_indices, 2]

      corr_mean  <- idf_mean_dv + interpolate_vector(diff_mean, steps = 101)
      corr_lower <- quants_dv[, 1] + interpolate_vector(diff_lower, steps = 101)
      corr_upper <- quants_dv[, 2] + interpolate_vector(diff_upper, steps = 101)

      res_list[[si]] <- list(mean = corr_mean, lower = corr_lower, upper = corr_upper)
      all_lower <- c(all_lower, corr_lower)
      all_upper <- c(all_upper, corr_upper)
    }


  } else if (fit$shp_d == "j") {
    cat("Computing j-dimensional shape IDF curve...\n")

    xi_j <- fit$pars$xi[, j]
    params_dv <- calc_mu_sigma(durs_ip, mut, sigma0, theta, eta)

    for (si in 1:length(prob_vec)) {
      p <- prob_vec[si]

      idf_samp <- gev_mat(p, params_dv$mu, params_dv$sigma, xi_j, "sample")

      corr_mean  <- rowMeans(idf_samp)
      quants     <- matrixStats::rowQuantiles(idf_samp, probs = quant_probs)
      corr_lower <- quants[, 1]
      corr_upper <- quants[, 2]

      res_list[[si]] <- list(mean = corr_mean, lower = corr_lower, upper = corr_upper)
      all_lower <- c(all_lower, corr_lower)
      all_upper <- c(all_upper, corr_upper)
    }
  }

  idf_xlims <- range(durs_ip)
  idf_ylims <- c(min(all_lower), max(all_upper) + 10)
  rp_paste  <- paste(rp, collapse = ",")

  plot(1, type = "n", xlim = idf_xlims, ylim = idf_ylims,
       main = paste0(fit$stationID[j], "\nfrequency = ", rp_paste, " years, xi shape = ", fit$shp_d),
       xaxt = "n",
       ylab = "Intensity [mm/h]", xlab = "Duration [h]", cex.main = 0.75, log = "xy")

  for (si in 1:length(prob_vec)) {
    lines(durs_ip, res_list[[si]]$mean, col = si_col[si], lwd = 2)

    if (length(prob_vec) == 1) {

      lines(durs_ip, res_list[[si]]$upper, lty = 2, col = "grey")
      lines(durs_ip, res_list[[si]]$lower, lty = 2, col = "grey")
      legend("topright", legend = c("credible \nintervals", "dGEV"),
             lwd = c(1, 2), lty = c(2, 1), col = c("grey", "black"))
    } else {
      polygon(x = c(durs_ip, rev(durs_ip)),
              y = c(res_list[[si]]$upper, rev(res_list[[si]]$lower)),
              border = NA, col = adjustcolor(si_col[si], alpha.f = 0.5))
    }
  }

  axis(side = 1, at = durs, labels = durs)

  if (length(prob_vec) > 1) {
    legend("topright", legend = paste("T =", rp, "years"),
           lwd = 1, lty = 1, col = si_col[1:length(rp)],
           pt.bg = adjustcolor(si_col[1:length(rp)], alpha.f = 0.5), pch = 15)
  }
}
