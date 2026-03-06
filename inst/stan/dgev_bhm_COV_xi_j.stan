functions {
  real gev_lpdf(real y, real mu, real sigma, real xi) {
    real z;
    if (fabs(xi) > 1e-10) {
      z = 1 + xi * (y - mu) / sigma;
      if (z <= 0) return negative_infinity();
      return -log(sigma) - (1 + 1/xi) * log(z) - pow(z, -1/xi);
    } else {
      z = (y - mu) / sigma;
      return -log(sigma) - z - exp(-z);
    }
  }
}

data {
  int<lower=0> J;                   // # of groups
  int<lower=0> D;                   // # of subgroups per group (fixed)
  int<lower=0> N_total;             // total amount of observations
  array[D,J] int<lower=0> s;        // Subgroup size for each group (s)
  vector[N_total] y;                // Flattened observations
  vector[D] d;
  vector[J] COV;
}

parameters {
  real mut_0;  // Fixed: consistent naming
  real mut_1;  // Fixed: consistent naming
  vector<lower=0>[J] sigma0;
  vector<lower=-1, upper=1>[J] xi;
  vector<lower=0>[J] theta;
  vector<lower=0, upper=1>[J] eta;
  real<lower=0> beta;
  real<lower=-0.5, upper=0.5> delta;
}

transformed parameters {
  vector[J] mut;
  for (j in 1:J) {
    mut[j] = mut_0 + mut_1 * COV[j];  // Fixed: use mut_0 and mut_1
  }
}

model {
  int pos = 1;
  
  mut_0 ~ normal(0, 10);
  mut_1 ~ normal(0, 5);
  sigma0 ~ gamma(beta * 10, 10);
  xi ~ normal(delta, 0.1);

  for (j in 1:J) {
    for (dd in 1:D) {
      real mu = mut[j] * sigma0[j] * pow(d[dd] + theta[j], -eta[j]);
      real sigma = sigma0[j] * pow(d[dd] + theta[j], -eta[j]);
      
      for (i in 1:s[dd,j]) {
        target += gev_lpdf(y[pos] | mu, sigma, xi[j]);
        pos += 1;  // Fixed: increment position correctly
      }
    }
  }
}