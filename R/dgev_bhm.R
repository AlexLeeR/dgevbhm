#' Bayesian Hierarchical Model for d-GEV
#'
#' Fits a Bayesian Hierarchical Model using Stan for duration-dependent
#' Generalized Extreme Value distributions.
#'
#' @param data A list containing processed station data (output from \code{data_read}).
#' @param chains Integer. Number of independent MCMC chains.
#' @param iter Integer. Number of iterations per chain (including warmup).
#' @param cores Integer. Number of cores for parallel sampling.
#' @param shp_d Character. Dimension for shape parameter: "d" for durational, "j" for spatial.
#'
#' @return A list containing MCMC samples (pars), original data, and metadata.
#'
#' @importFrom rstan stan_model sampling extract
#' @importFrom IDF gev.d.fit
#' @export
dgev_bhm = function(data, chains = 4, iter = 2000, cores = 4, shp_d = NULL) {
  # No require() needed here; we define dependencies in the DESCRIPTION file

  if (is.null(shp_d) || (shp_d != "d" && shp_d != "j")) {
    stop("please enter \"d\" durational or \"j\" spatial for shape parameter dimension")
  }

  if (cores > chains) {
    warning(sprintf("cores = %d are specified, but chains = %d. Stan can only use one core per chain; %d cores will remain idle.",
                    cores, chains, cores - chains))
  }

  cat("\nComputing initial parameters... ")

  data_list = data$y
  durs      = data$durs

  unique_dims = unique(sapply(data_list, ncol))
  if (length(unique_dims) > 1) {
    stop("Duration length mismatch across stations")
  }

  dgev_params_list = lapply(data_list, function(x) {
    if (length(unique(sapply(data_list, function(x) dim(x)[2]))) > 1) {
      stop("Duration length mismatch")
    }
    ds_durs    = rep(durs, each = dim(x)[1])
    dat        = as.vector(x)

    # Use IDF:: to call functions from the IDF package
    fit        = IDF::gev.d.fit(xdat = dat, ds = ds_durs, show = FALSE)
    fit_params = fit$mle
    return(fit_params)
  })

  dgev_params = array(unlist(dgev_params_list),
                      dim = c(length(dgev_params_list[[1]]),
                              length(dgev_params_list)))

  cat("Complete!\n")

  # Common data preparation for both cases
  y_flat = unlist(data_list)
  d_len = length(durs)
  s_dj = matrix(rep(sapply(data_list, nrow), d_len), nrow = d_len, byrow = TRUE)

  stan_data = list(J = length(data_list), D = ncol(data_list[[1]]),
                   N_total = length(y_flat), s = s_dj, y = y_flat,
                   d = durs)

  if (shp_d == "d") {
    init_fun <- function(...) {
      list(mut    = dgev_params[1, ],
           sigma0 = dgev_params[2, ],
           theta  = dgev_params[4, ],
           eta    = dgev_params[5, ],
           xi     = rep(mean(dgev_params[3, ]), length(durs)))
    }

    cat("\nReading Stan file (durational)... ")

    stan_model_path <- system.file("stan", "dgev_bhm_xi_d.stan", package = "dgevbhm")
    if (stan_model_path == "") {
      stop("Stan model file not found. Ensure it is in inst/stan/ and re-install.")
    }
    mod = rstan::stan_model(stan_model_path)
    cat("Complete!\n")

  } else if (shp_d == "j") {
    init_fun <- function(...) {
      list(mut    = dgev_params[1, ],
           sigma0 = dgev_params[2, ],
           theta  = dgev_params[4, ],
           eta    = dgev_params[5, ],
           xi     = dgev_params[3, ])
    }

    cat("\nReading Stan file (spatial)... ")

    stan_model_path <- system.file("stan", "dgev_bhm_xi_j.stan", package = "dgevbhm")
    if (stan_model_path == "") {
      stop("Stan model file not found. Ensure it is in inst/stan/ and re-install.")
    }
    mod <- rstan::stan_model(stan_model_path)

    cat("Complete!\n")
  }

  fit <- rstan::sampling(mod,
                         data = stan_data,
                         init = init_fun,
                         chains = chains,
                         iter = iter,
                         warmup = iter/2,
                         cores = cores)

  warmup <- iter / 2
  samples_per_chain <- iter - warmup

  return(list(pars = rstan::extract(fit),
              data = data_list,
              durs = durs,
              j = length(data_list),
              d = length(durs),
              shp_d = shp_d,
              mcsamp = chains * samples_per_chain,
              stationID = data$stationID))
}
