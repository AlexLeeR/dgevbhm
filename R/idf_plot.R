#' Quantile function for d-GEV using Torch
#'
#' @param p Probability
#' @param mu Location parameter
#' @param sigma Scale parameter
#' @param xi Shape parameter
#' @importFrom torch torch_tensor torch_where torch_log
#' @export
dgev_torch <- function(p, mu, sigma, xi) {
  p <- torch::torch_tensor(p, dtype = torch::torch_float64())
  one <- torch::torch_tensor(1.0, dtype = torch::torch_float64())

  q <- torch::torch_where(
    xi != torch::torch_tensor(0, dtype = torch::torch_float64()),
    mu - (sigma / xi) * (one - (-torch::torch_log(p))^(-xi)),
    mu - sigma * torch::torch_log(-torch::torch_log(p))
  )
  return(q)
}

#' Plot IDF Curves using Torch
#'
#' @param fit List. The output object from the dgev_bhm function.
#' @param j Integer. Station index.
#' @param rp Numeric vector. Return periods (max 3).
#' @param alpha Numeric. Significance level.
#' @param cores Integer. Threads for torch.
#'
#' @importFrom torch torch_tensor torch_where torch_log torch_empty torch_linspace torch_cat torch_pow torch_mean torch_quantile torch_stack torch_min torch_max torch_set_num_threads
#' @importFrom graphics plot lines polygon legend
#' @export
idf_plot_torch <- function(fit, j=NULL, rp, alpha = 0.05, cores = 4) {

  if (length(rp) > 3) {
    stop('No than 3 return periods allowed!')
  }

  if (alpha < 0 | alpha > 1) {
    stop("alpha should be less than 1 and greater than 0")
  }

  torch::torch_set_num_threads(cores)

  prob <- 1 - (1 / rp)
  prob_t <- torch::torch_tensor(prob, dtype = torch::torch_float64())

  durs <- fit$durs
  if (length(durs) <= 3) {
    stop("more than 3 durations needed")
  }

  si_col <- c("#D81B60", "#1E88E5", "#004D40")

  durs_ip <- torch::torch_empty(0, dtype = torch::torch_float64())
  durs_t <- torch::torch_tensor(durs, dtype = torch::torch_float64())

  for (dv in 2:length(durs)) {
    if (dv == 2) {
      durs_sq <- torch::torch_linspace(durs_t[dv - 1], durs_t[dv], steps = 101)
      durs_ip <- torch::torch_cat(list(durs_ip, durs_sq))
    } else {
      durs_sq <- torch::torch_linspace(durs_t[dv - 1], durs_t[dv], steps = 101)
      durs_ip <- torch::torch_cat(list(durs_ip, durs_sq[2:101]))
    }
  }

  mut <- torch::torch_tensor(fit$pars$mut[, j], dtype = torch::torch_float64())
  sigma0 <- torch::torch_tensor(fit$pars$sigma0[, j], dtype = torch::torch_float64())
  theta <- torch::torch_tensor(fit$pars$theta[, j], dtype = torch::torch_float64())
  eta <- torch::torch_tensor(fit$pars$eta[, j], dtype = torch::torch_float64())

  q_upper <- torch::torch_tensor(1 - (alpha / 2), dtype = torch::torch_float64())
  q_lower <- torch::torch_tensor(0 + (alpha / 2), dtype = torch::torch_float64())

  if (fit$shp_d == "d") {
    cat("Computing d-dimensional shape IDF curve...\n")

    xi <- torch::torch_tensor(fit$pars$xi, dtype = torch::torch_float64())
    xi_mean <- torch::torch_mean(xi, dim = 1)

    xi_mean_ip <- torch::torch_empty(0, dtype = torch::torch_float64())

    for (dx in 2:length(xi_mean)) {
      if (dx == 2) {
        seq_points <- torch::torch_linspace(xi_mean[dx - 1], xi_mean[dx], steps = 101)
        xi_mean_ip <- torch::torch_cat(list(xi_mean_ip, seq_points))
      } else {
        seq_points <- torch::torch_linspace(xi_mean[dx - 1], xi_mean[dx], steps = 101)
        xi_mean_ip <- torch::torch_cat(list(xi_mean_ip, seq_points[2:101]))
      }
    }

    sigma_d_dv <- sigma0$unsqueeze(2) * torch::torch_pow(durs_ip$unsqueeze(1) + theta$unsqueeze(2), -eta$unsqueeze(2))
    mu_d_dv <- sigma_d_dv * mut$unsqueeze(2)

    sigma_d_actual <- sigma0$unsqueeze(2) * torch::torch_pow(durs_t$unsqueeze(1) + theta$unsqueeze(2), -eta$unsqueeze(2))
    mu_d_actual <- sigma_d_actual * mut$unsqueeze(2)

    if (prob_t$shape == 1) {

      idf_samp_dv <- dgev_torch(prob_t$unsqueeze(1),
                                mu_d_dv, sigma_d_dv,
                                xi_mean_ip$unsqueeze(1))

      idf_mean_dv <- torch::torch_mean(idf_samp_dv, dim = 1)
      idf_q_upper_dv <- torch::torch_quantile(idf_samp_dv, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower_dv <- torch::torch_quantile(idf_samp_dv, q = q_lower, dim = 1)$squeeze(1)

      idf_samp_actual <- dgev_torch(prob_t$unsqueeze(1), mu_d_actual, sigma_d_actual, xi)

      idf_mean_actual <- torch::torch_mean(idf_samp_actual, dim = 1)
      idf_q_upper_actual <- torch::torch_quantile(idf_samp_actual, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower_actual <- torch::torch_quantile(idf_samp_actual, q = q_lower, dim = 1)$squeeze(1)

      dur_indices <- which(as.array(durs_ip) %in% as.array(durs_t))

      idf_mean_dv_di <- idf_mean_dv[dur_indices]
      idf_q_upper_dv_di <- idf_q_upper_dv[dur_indices]
      idf_q_lower_dv_di <- idf_q_lower_dv[dur_indices]

      diff_idf_mean <- idf_mean_actual - idf_mean_dv_di
      diff_idf_q_upper <- idf_q_upper_actual - idf_q_upper_dv_di
      diff_idf_q_lower <- idf_q_lower_actual - idf_q_lower_dv_di

      diff_idf_mean_ip <- torch::torch_empty(0, dtype = torch::torch_float64())
      diff_idf_q_upper_ip <- torch::torch_empty(0, dtype = torch::torch_float64())
      diff_idf_q_lower_ip <- torch::torch_empty(0, dtype = torch::torch_float64())

      for (dx in 2:length(diff_idf_mean)) {
        if (dx == 2) {
          diff_idf_mean_sp <- torch::torch_linspace(diff_idf_mean[dx - 1], diff_idf_mean[dx], steps = 101)
          diff_idf_mean_ip <- torch::torch_cat(list(diff_idf_mean_ip, diff_idf_mean_sp))

          diff_idf_q_upper_sp <- torch::torch_linspace(diff_idf_q_upper[dx - 1], diff_idf_q_upper[dx], steps = 101)
          diff_idf_q_upper_ip <- torch::torch_cat(list(diff_idf_q_upper_ip, diff_idf_q_upper_sp))

          diff_idf_q_lower_sp <- torch::torch_linspace(diff_idf_q_lower[dx - 1], diff_idf_q_lower[dx], steps = 101)
          diff_idf_q_lower_ip <- torch::torch_cat(list(diff_idf_q_lower_ip, diff_idf_q_lower_sp))
        } else {
          diff_idf_mean_sp <- torch::torch_linspace(diff_idf_mean[dx - 1], diff_idf_mean[dx], steps = 101)
          diff_idf_mean_ip <- torch::torch_cat(list(diff_idf_mean_ip, diff_idf_mean_sp[2:101]))

          diff_idf_q_upper_sp <- torch::torch_linspace(diff_idf_q_upper[dx - 1], diff_idf_q_upper[dx], steps = 101)
          diff_idf_q_upper_ip <- torch::torch_cat(list(diff_idf_q_upper_ip, diff_idf_q_upper_sp[2:101]))

          diff_idf_q_lower_sp <- torch::torch_linspace(diff_idf_q_lower[dx - 1], diff_idf_q_lower[dx], steps = 101)
          diff_idf_q_lower_ip <- torch::torch_cat(list(diff_idf_q_lower_ip, diff_idf_q_lower_sp[2:101]))
        }
      }

      idf_mean_corrected <- idf_mean_dv + diff_idf_mean_ip
      idf_q_upper_corrected <- idf_q_upper_dv + diff_idf_q_upper_ip
      idf_q_lower_corrected <- idf_q_lower_dv + diff_idf_q_lower_ip

      idf_ylims <- c(min(as.numeric(idf_q_lower_corrected)), max(as.numeric(idf_q_upper_corrected) + 10))

      plot(as.numeric(durs_ip), as.numeric(idf_mean_corrected), ylim = idf_ylims, type = "l", lwd = 2, log = "xy",
           main = paste0(fit$stationID[j], "\nfrequency = ", rp, " years, xi shape = ", fit$shp_d),
           ylab = "Intensity [mm/h]", xlab = "Duration [h]", cex.main = 0.75)
      lines(as.numeric(durs_ip),
            as.numeric(idf_q_upper_corrected),
            lty = 2, col = "grey")
      lines(as.numeric(durs_ip),
            as.numeric(idf_q_lower_corrected),
            lty = 2, col = "grey")
      legend("topright",
             legend = c("confidence \nintervals", "dGEV"),
             lwd = c(1, 2), lty = c(2, 1), col = c("grey", "black"))
    } else if (prob_t$shape > 1 & prob_t$shape <= 3) {

      idf_samp_dv <- dgev_torch(prob_t$unsqueeze(1)$unsqueeze(1),
                                mu_d_dv$unsqueeze(3),
                                sigma_d_dv$unsqueeze(3),
                                xi_mean_ip$unsqueeze(2)$unsqueeze(1))

      idf_mean_dv <- torch::torch_mean(idf_samp_dv, dim = 1)
      idf_q_upper_dv <- torch::torch_quantile(idf_samp_dv, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower_dv <- torch::torch_quantile(idf_samp_dv, q = q_lower, dim = 1)$squeeze(1)

      idf_samp_actual <- dgev_torch(prob_t$unsqueeze(1)$unsqueeze(1),
                                    mu_d_actual$unsqueeze(3),
                                    sigma_d_actual$unsqueeze(3),
                                    xi$unsqueeze(3))

      idf_mean_actual <- torch::torch_mean(idf_samp_actual, dim = 1)
      idf_q_upper_actual <- torch::torch_quantile(idf_samp_actual, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower_actual <- torch::torch_quantile(idf_samp_actual, q = q_lower, dim = 1)$squeeze(1)

      dur_indices <- which(as.array(durs_ip) %in% as.array(durs_t))

      idf_mean_dv_di <- idf_mean_dv[dur_indices, ]
      idf_q_upper_dv_di <- idf_q_upper_dv[dur_indices, ]
      idf_q_lower_dv_di <- idf_q_lower_dv[dur_indices, ]

      diff_idf_mean <- idf_mean_actual - idf_mean_dv_di
      diff_idf_q_upper <- idf_q_upper_actual - idf_q_upper_dv_di
      diff_idf_q_lower <- idf_q_lower_actual - idf_q_lower_dv_di

      diff_idf_mean_ip_si_list <- list()
      diff_idf_q_upper_ip_si_list <- list()
      diff_idf_q_lower_ip_si_list <- list()

      for (si in 1:prob_t$shape) {

        diff_idf_mean_ip <- torch::torch_empty(0, dtype = torch::torch_float64())
        diff_idf_q_upper_ip <- torch::torch_empty(0, dtype = torch::torch_float64())
        diff_idf_q_lower_ip <- torch::torch_empty(0, dtype = torch::torch_float64())

        for (dx in 2:length(diff_idf_mean[, 1])) {
          if (dx == 2) {

            diff_idf_mean_sp <- torch::torch_linspace(diff_idf_mean[dx - 1, si], diff_idf_mean[dx, si], steps = 101)
            diff_idf_mean_ip <- torch::torch_cat(list(diff_idf_mean_ip, diff_idf_mean_sp))

            diff_idf_q_upper_sp <- torch::torch_linspace(diff_idf_q_upper[dx - 1, si], diff_idf_q_upper[dx, si], steps = 101)
            diff_idf_q_upper_ip <- torch::torch_cat(list(diff_idf_q_upper_ip, diff_idf_q_upper_sp))

            diff_idf_q_lower_sp <- torch::torch_linspace(diff_idf_q_lower[dx - 1, si], diff_idf_q_lower[dx, si], steps = 101)
            diff_idf_q_lower_ip <- torch::torch_cat(list(diff_idf_q_lower_ip, diff_idf_q_lower_sp))
          } else {
            diff_idf_mean_sp <- torch::torch_linspace(diff_idf_mean[dx - 1, si], diff_idf_mean[dx, si], steps = 101)
            diff_idf_mean_ip <- torch::torch_cat(list(diff_idf_mean_ip, diff_idf_mean_sp[2:101]))

            diff_idf_q_upper_sp <- torch::torch_linspace(diff_idf_q_upper[dx - 1, si], diff_idf_q_upper[dx, si], steps = 101)
            diff_idf_q_upper_ip <- torch::torch_cat(list(diff_idf_q_upper_ip, diff_idf_q_upper_sp[2:101]))

            diff_idf_q_lower_sp <- torch::torch_linspace(diff_idf_q_lower[dx - 1, si], diff_idf_q_lower[dx, si], steps = 101)
            diff_idf_q_lower_ip <- torch::torch_cat(list(diff_idf_q_lower_ip, diff_idf_q_lower_sp[2:101]))
          }
        }
        diff_idf_mean_ip_si_list[[si]] <- diff_idf_mean_ip
        diff_idf_q_upper_ip_si_list[[si]] <- diff_idf_q_upper_ip
        diff_idf_q_lower_ip_si_list[[si]] <- diff_idf_q_lower_ip
      }
      diff_idf_mean_ip_si <- torch::torch_stack(diff_idf_mean_ip_si_list, dim = 2)
      diff_idf_q_upper_ip_si <- torch::torch_stack(diff_idf_q_upper_ip_si_list, dim = 2)
      diff_idf_q_lower_ip_si <- torch::torch_stack(diff_idf_q_lower_ip_si_list, dim = 2)

      idf_mean_corrected <- idf_mean_dv + diff_idf_mean_ip_si
      idf_q_upper_corrected <- idf_q_upper_dv + diff_idf_q_upper_ip_si
      idf_q_lower_corrected <- idf_q_lower_dv + diff_idf_q_lower_ip_si

      idf_xlims <- c(as.numeric(torch::torch_min(durs_ip)), as.numeric(torch::torch_max(durs_ip)))
      idf_ylims <- c(as.numeric(torch::torch_min(idf_q_lower_corrected)), as.numeric(torch::torch_max(idf_q_upper_corrected)))

      rp_paste <- paste(rp, collapse = ",")

      plot(1, type = "n", xlim = idf_xlims, ylim = idf_ylims,
           main = paste0(fit$stationID[j], "\nfrequency = ",
                         rp_paste, " years, xi shape = ", fit$shp_d),
           ylab = "Intensity [mm/h]", xlab = "Duration [h]", cex.main = 0.75, log = "xy")
      for (si in 1:prob_t$shape) {
        lines(as.numeric(durs_ip), as.numeric(idf_mean_corrected[, si]), col = si_col[si], lwd = 2)
        polygon(x = c(as.numeric(durs_ip), rev(as.numeric(durs_ip))),
                y = c(as.numeric(idf_q_upper_corrected[, si]), rev(as.numeric(idf_q_lower_corrected[, si]))),
                border = NA, col = paste0(si_col[si], "80"))

      }
      l_si <- rep(1, prob_t$shape)
      legend("topright", legend = paste("T =", rp, "years"), lwd = l_si,
             lty = l_si, col = si_col, pt.bg = paste0(si_col, "80"), pch = 15)
    }
  } else if (fit$shp_d == "j") {
    cat("Computing j-dimensional shape IDF curve...\n")

    xi_j <- torch::torch_tensor(fit$pars$xi[, j], dtype = torch::torch_float64())

    sigma_d_dv <- sigma0$unsqueeze(2) * torch::torch_pow(durs_ip$unsqueeze(1) + theta$unsqueeze(2), -eta$unsqueeze(2))
    mu_d_dv <- sigma_d_dv * mut$unsqueeze(2)

    if (prob_t$shape == 1) {

      idf_samp <- dgev_torch(prob_t$unsqueeze(1), mu_d_dv, sigma_d_dv, xi_j$unsqueeze(2))

      idf_mean <- torch::torch_mean(idf_samp, dim = 1)
      idf_q_upper <- torch::torch_quantile(idf_samp, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower <- torch::torch_quantile(idf_samp, q = q_lower, dim = 1)$squeeze(1)

      idf_xlims <- c(as.numeric(torch::torch_min(durs_ip)),
                     as.numeric(torch::torch_max(durs_ip)))
      idf_ylims <- c(as.numeric(torch::torch_min(idf_q_lower)),
                     as.numeric(torch::torch_max(idf_q_upper)))

      plot(as.numeric(durs_ip), as.numeric(idf_mean),
           ylim = idf_ylims, type = "l", lwd = 2, log = "xy",
           main = paste0(fit$stationID[j], "\nfrequency = ",
                         rp, " years, xi shape = ", fit$shp_d),
           ylab = "Intensity [mm/h]", xlab = "Duration [h]", cex.main = 0.75)
      lines(as.numeric(durs_ip),
            as.numeric(idf_q_upper),
            lty = 2, col = "grey")
      lines(as.numeric(durs_ip),
            as.numeric(idf_q_lower),
            lty = 2, col = "grey")
      legend("topright",
             legend = c("confidence \nintervals", "dGEV"),
             lwd = c(1, 2), lty = c(2, 1), col = c("grey", "black"))
    } else if (prob_t$shape > 1 & prob_t$shape <= 3) {

      rs <- fit$mcsamp

      idf_samp <- dgev_torch(prob_t$unsqueeze(1)$unsqueeze(1),
                             mu_d_dv$unsqueeze(3), sigma_d_dv$unsqueeze(3),
                             xi_j$view(c(rs, 1, 1)))

      idf_mean <- torch::torch_mean(idf_samp, dim = 1)
      idf_q_upper <- torch::torch_quantile(idf_samp, q = q_upper, dim = 1)$squeeze(1)
      idf_q_lower <- torch::torch_quantile(idf_samp, q = q_lower, dim = 1)$squeeze(1)

      idf_xlims <- c(as.numeric(torch::torch_min(durs_ip)),
                     as.numeric(torch::torch_max(durs_ip)))
      idf_ylims <- c(as.numeric(torch::torch_min(idf_q_lower)),
                     as.numeric(torch::torch_max(idf_q_upper)))

      rp_paste <- paste(rp, collapse = ",")

      plot(1, type = "n", xlim = idf_xlims, ylim = idf_ylims,
           main = paste0(fit$stationID[j], "\nfrequency = ",
                         rp_paste, " years, xi shape = ", fit$shp_d),
           ylab = "Intensity [mm/h]", xlab = "Duration [h]", cex.main = 0.75, log = "xy")
      for (si in 1:prob_t$shape) {
        lines(as.numeric(durs_ip), as.numeric(idf_mean[, si]), col = si_col[si], lwd = 2)
        polygon(x = c(as.numeric(durs_ip), rev(as.numeric(durs_ip))),
                y = c(as.numeric(idf_q_upper[, si]), rev(as.numeric(idf_q_lower[, si]))),
                border = NA, col = paste0(si_col[si], "80"))

      }
      l_si <- rep(1, prob_t$shape)
      legend("topright", legend = paste("T =", rp, "years"),
             lwd = l_si, lty = l_si, col = si_col, pt.bg = paste0(si_col, "80"), pch = 15)
    }
  }
}
