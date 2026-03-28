#' Robust GEV Log-PDF using Torch
#'
#' Calculates the log-probability density function for the GEV distribution.
#' Handles the Gumbel limit (xi near 0) and invalid support areas robustly.
#'
#' @param x Tensor. Observed data.
#' @param mu Tensor. Location parameter.
#' @param sigma Tensor. Scale parameter.
#' @param xi Tensor. Shape parameter.
#'
#' @importFrom torch torch_tensor torch_log torch_exp torch_reciprocal torch_where torch_abs
#' @export
gev_lpdf_robust <- function(x, mu, sigma, xi) {
  one     <- torch::torch_tensor(1,      dtype = torch::torch_float64())
  neg_one <- torch::torch_tensor(-1,     dtype = torch::torch_float64())
  epsilon <- torch::torch_tensor(1e-10, dtype = torch::torch_float64())

  z <- (x - mu) / sigma
  T_term <- one + xi * z
  lpdf_gumbel <- -torch::torch_log(sigma) - torch::torch_exp(-z) - z
  xi_reciprocal <- torch::torch_reciprocal(xi)

  lpdf_non_gumbel <- -torch::torch_log(sigma) +
    (neg_one * T_term ^ (neg_one * xi_reciprocal)) +
    (neg_one * (one + xi_reciprocal) * torch::torch_log(T_term))

  neg_inf <- torch::torch_tensor(-Inf, dtype = torch::torch_float64())

  lpdf_non_gumbel <- torch::torch_where(
    T_term <= torch::torch_tensor(0, dtype = torch::torch_float64()),
    neg_inf,
    lpdf_non_gumbel)

  result <- torch::torch_where(torch::torch_abs(xi) < epsilon,
                               lpdf_gumbel,
                               lpdf_non_gumbel)
  return(result)
}

#' Calculate WAIC for d-GEV BHM
#'
#' Computes the Watanabe-Akaike Information Criterion (WAIC) for the fitted
#' Bayesian Hierarchical Model.
#'
#' @param fit List. Output from \code{dgev_bhm}.
#' @param cores Integer. Number of threads for torch (default 4).
#'
#' @return Numeric. The calculated WAIC value.
#'
#' @importFrom stats var
#' @importFrom torch torch_tensor torch_pow torch_stack torch_flatten torch_sum
#' @export
waic = function(fit, cores = 4) {
  # dependencies are in DESCRIPTION

  torch::torch_set_num_threads(cores)

  lpdf_list = lapply(1:length(fit$data), function(i) {
    mut_t    = torch::torch_tensor(fit$pars$mut[,i],    dtype = torch::torch_float64())
    sigma0_t = torch::torch_tensor(fit$pars$sigma0[,i], dtype = torch::torch_float64())
    eta_t    = torch::torch_tensor(fit$pars$eta[,i],    dtype = torch::torch_float64())
    theta_t  = torch::torch_tensor(fit$pars$theta[,i],  dtype = torch::torch_float64())
    data_t   = torch::torch_tensor(t(fit$data[[i]]),    dtype = torch::torch_float64())
    durs_t   = torch::torch_tensor(fit$durs,            dtype = torch::torch_float64())

    mut_t1    = mut_t$unsqueeze(2)$unsqueeze(3)
    sigma0_t1 = sigma0_t$unsqueeze(2)$unsqueeze(3)
    eta_t1    = eta_t$unsqueeze(2)$unsqueeze(3)
    theta_t1  = theta_t$unsqueeze(2)$unsqueeze(3)

    data_t1   = data_t$unsqueeze(1)
    durs_t1   = durs_t$unsqueeze(1)$unsqueeze(3)

    if (fit$shp_d == 'd') {
      xi_t  = torch::torch_tensor(fit$pars$xi, dtype = torch::torch_float64())
      xi_t1 = xi_t$unsqueeze(3)
    } else if (fit$shp_d == 'j') {
      xi_t  = torch::torch_tensor(fit$pars$xi[,i], dtype = torch::torch_float64())
      xi_t1 = xi_t$unsqueeze(2)$unsqueeze(3)
    }

    sigma_d_t  = sigma0_t1 * torch::torch_pow(durs_t1 + theta_t1, -eta_t1)
    mu_d_t     = mut_t1 * sigma_d_t

    lpdf_t = gev_lpdf_robust(data_t1, mu_d_t, sigma_d_t, xi_t1)
    return(torch::torch_flatten(lpdf_t, start_dim = 2, end_dim = 3))
  })

  log_lik_mat <- torch::torch_cat(lpdf_list, dim = 2)
  lse_i       <- torch::torch_logsumexp(log_lik_mat, dim = 1)

  S           <- dim(log_lik_mat)[1]
  S_tensor    <- torch::torch_tensor(S, dtype = torch::torch_float64())

  lppd_i      <- lse_i - torch::torch_log(S_tensor)
  lppd        <- torch::torch_sum(lppd_i)

  p_waic_i    <- torch::torch_var(log_lik_mat, dim = 1, unbiased = TRUE)
  p_waic      <- torch::torch_sum(p_waic_i)

  negtwo      <- torch::torch_tensor(-2, dtype = torch::torch_float64())
  waic        <- negtwo * (lppd - p_waic)
  waic_val    <- as.numeric(waic)

  return(waic_val)
}
