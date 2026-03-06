#' Read and Process Precipitation Station Data
#'
#' This function reads .txt files from a directory, extracts metadata,
#' and calculates annual maximum precipitation for specified durations.
#'
#' @param dir_path Character. Path to the directory containing .txt files.
#' @param missing_data Numeric. Percentage threshold for missing data (default is 5).
#' @param durs Numeric vector. Durations (in hours) to calculate max precipitation.
#' @param cores Integer. Number of CPU cores for parallel processing.
#' @param min_year Integer. Minimum number of valid years required to keep a station.
#'
#' @return A list containing processed precipitation arrays, annual mean precipitation,
#' coordinates, and station IDs.
#'
#' @importFrom stats setNames
#' @importFrom parallel mclapply
#' @export
data_read = function(dir_path,
                     missing_data = 5,
                     durs = c(1, 3, 6, 12, 24, 48),
                     cores = 4,
                     min_year = 25) {
  # dependencies are handled in DESCRIPTION, not via require()

  file_list = list.files(dir_path, pattern = "\\.txt$", full.names = TRUE)

  # Pre-calculate common values
  missing_threshold = 1 - (missing_data / 100)

  message("Processing ", length(file_list), " files...")

  data_list = parallel::mclapply(seq_along(file_list), function(fp) {
    file_contents = readLines(file_list[fp], warn = FALSE)

    # Extract all metadata at once
    metadata_lines = grep("^(Start datetime|End datetime|Latitude|Longitude|Station ID):", file_contents, value = TRUE)
    metadata = strsplit(metadata_lines, ": ")
    metadata_dict = setNames(sapply(metadata, `[`, 2), sapply(metadata, `[`, 1))

    # Time processing using lubridate::
    start_date_ut = as.numeric(lubridate::ymd_h(metadata_dict[["Start datetime"]]))
    end_date_ut = as.numeric(lubridate::ymd_h(metadata_dict[["End datetime"]]))

    start_year = lubridate::year(lubridate::as_datetime(start_date_ut))
    end_year = lubridate::year(lubridate::as_datetime(end_date_ut))
    year_range = start_year:end_year

    # Precip data processing
    pos = grep("^Other:", file_contents, ignore.case = TRUE)[1]
    precip = as.numeric(file_contents[(pos + 1):length(file_contents)])
    precip[precip == -999] = NA

    # Create time sequence
    actual_hours = seq(start_date_ut, end_date_ut, by = 3600)
    actual_years = lubridate::year(lubridate::as_datetime(actual_hours))

    # Process each year
    year_data = lapply(year_range, function(yr) {
      year_mask = actual_years == yr
      yearly_precip = precip[year_mask]

      valid_count = sum(!is.na(yearly_precip))
      rain_perc = valid_count / length(yearly_precip)

      if (rain_perc > missing_threshold && valid_count > 0) {
        yearly_precip[is.na(yearly_precip)] = 0

        precip_durs_list = lapply(durs, function(k) {
          if (length(yearly_precip) < k) return(NA)

          if (k == 1) {
            return(max(yearly_precip, na.rm = TRUE))
          } else {
            # Use zoo:: for rolling operations
            roll_means = zoo::rollmean(yearly_precip, k = k, align = "left", fill = NA)
            roll_means_clean = roll_means[!is.na(roll_means)]
            if (length(roll_means_clean) > 0) {
              return(max(roll_means_clean, na.rm = TRUE))
            } else {
              return(NA)
            }
          }
        })

        return(list(
          precip_durs1 = unlist(precip_durs_list),
          mean_precip1 = mean(yearly_precip, na.rm = TRUE)
        ))
      } else {
        return(list(precip_durs1 = NULL, mean_precip1 = NULL))
      }
    })

    precip_durs = lapply(year_data, "[[", "precip_durs1")
    mean_precip = lapply(year_data, "[[", "mean_precip1")
    valid_years = !sapply(precip_durs, is.null)

    if (any(valid_years)) {
      filtered_precip = precip_durs[valid_years]
      complete_years = sapply(filtered_precip, function(x) {
        length(x) == length(durs) && all(!is.na(x))
      })

      if (any(complete_years)) {
        filtered_precip_clean = filtered_precip[complete_years]
        final_array = t(matrix(unlist(filtered_precip_clean), nrow = length(durs)))
        annual_mean_precip = mean(unlist(mean_precip[valid_years][complete_years]), na.rm = TRUE)
        lat = as.numeric(metadata_dict[["Latitude"]])
        lon = as.numeric(metadata_dict[["Longitude"]])
        stationID_read = metadata_dict[["Station ID"]]
      } else {
        final_array = NA
        annual_mean_precip = NA
        lat = as.numeric(metadata_dict[["Latitude"]])
        lon = as.numeric(metadata_dict[["Longitude"]])
        stationID_read = metadata_dict[["Station ID"]]
      }
    } else {
      final_array = NA; annual_mean_precip = NA; lat = NA; lon = NA; stationID_read = NA
    }

    return(list(y = final_array, amp = annual_mean_precip, lat = lat, lon = lon, stationID = stationID_read))

  }, mc.cores = cores)

  # Filter and combine results
  valid_stations = sapply(data_list, function(x) {
    is.array(x$y) && dim(x$y)[1] >= min_year && all(!is.na(x$y))
  })

  if (sum(valid_stations) == 0) {
    warning("No valid stations found with >= ", min_year, " years of complete data")
    return(list(y = list(), amp = numeric(), latlon = matrix(ncol=2), stationID = character(), durs = durs))
  }

  filtered_data = data_list[valid_stations]

  final_list = list(
    y = lapply(filtered_data, `[[`, "y"),
    amp = sapply(filtered_data, `[[`, "amp"),
    latlon = do.call(rbind, lapply(filtered_data, function(x) c(x$lat, x$lon))),
    stationID = sapply(filtered_data, `[[`, "stationID"),
    durs = durs
  )

  message("Done! Processed ", sum(valid_stations), " valid stations.")
  return(final_list)
}
