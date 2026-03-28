#' UK Extreme Rainfall Data
#'
#' A dataset containing annual maxima for 47 stations.
#' The full list of processed data.
#'
#' @format A list with elements:
#' \describe{
#'   \item{y}{A list of matrices duration dependent annual maxima for each station. year x duration}
#'   \item{amp}{Annual mean precipitation for each station.}
#'   \item{latlon}{Latitude and longitude for each station.}
#'   \item{stationID}{station ID for each station.}
#'   \item{durs}{Durations corresponding to the columns of the y matrices.}
#' }
"UK_data_full"

#' UK Extreme Rainfall Data
#'
#' A dataset containing annual maxima for 25 stations.
#' Filtered using amp median.
#'
#' @format A list with elements:
#' \describe{
#'   \item{y}{A list of matrices duration dependent annual maxima for each station. year x duration}
#'   \item{amp}{Annual mean precipitation for each station.}
#'   \item{latlon}{Latitude and longitude for each station.}
#'   \item{stationID}{station ID for each station.}
#'   \item{durs}{Durations corresponding to the columns of the y matrices.}
#' }
"UK_data"

#' Result of xi_d Model Fit
#'
#' Pre-computed BHM fit using the shape parameter hierarchical with respect to duration.
"fit_xi_d"

#' Result of xi_d Model Fit
#'
#' Pre-computed BHM fit using the shape parameter hierarchical with respect to space.
"fit_xi_j"
