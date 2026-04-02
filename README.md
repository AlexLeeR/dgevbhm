The `dgevbhm` package is an R and Stan-based package for Bayesian hierarchical regional modeling of extreme precipitation. 
It has the ability to provide robust intensity-duration-frequency (IDF) curve estimations in data-scarce scenarios (e.g. records under 50 years).

The requirements are R (>=4.4.0), `rstan` and `torch`.

To install please run: `pak::pak("AlexLeeR/dgevbhm")`

Basic workflow:
1. **Data pre-processing:** Use `data_read()` to format the sub-daily precipitation TXTs
2. **Model Fitting:** Use `dgev_bhm()` to fit the Bayesian hierarchical model with processed data
3. **Model evaluation and application tools:** Utilize diagnostics tools such as `pp_plot_torch()` and `waic`, and application tools like `idf_plot_torch()`.

```r
library(dgevbhm)
data('UK_data')

# Fit a model with a shape parameter hierarchical over durations
fit <- dgev_bhm(UK_data, iter=2000, cores=8, shp_d='d')

# Generate an IDF plot for the second locality
idf_plot_torch(fit, j=2, rp=c(2, 10, 100))
```

The `data_read()` function...
```r
data_processed <- data_read(
  dir_path = 'path/to/directory',  # path to the directory containing TXTs
  missing_data = 5,                # minimum missing data per year (in percent)
  durs = c(1,3,6,12,24,24,48),     # durations e.g. in hours
  cores = 4,                       # cores for parallel processing
  min_year = 25                    # minimum years for station to be selected
)
```

This project is licensed under the MIT License.

For questions or support, please contact: Alexander Lee Rischmuller - alex.rischmuller@gmail.com
