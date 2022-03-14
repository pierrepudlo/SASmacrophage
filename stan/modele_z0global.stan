functions {
  vector path(int N, vector tps, vector z0, vector a, vector log_b) {
    vector[N] y = log_b + (z0 - log_b) .* exp(-a .* tps);
    return y;
  }

  real pos_normal_rng(real mu, real sigma){
    real u = uniform_rng(0, 1);
    real z0 = - mu/sigma;
    real z = inv_Phi((1 - Phi(z0)) * u + Phi(z0));
    real y = mu + sigma * z;
    return y;
  }
}

data {
  // Hierarchical data
  int<lower=0> n_groups;
  int<lower=0> n_ind;
  int<lower=1, upper=n_groups> g[n_ind];
  int<lower=0> N;
  int<lower=1, upper=n_ind> id[N];
  vector[N] z;
  vector[N] tps;
  vector[N] lambda;

  // Hyperparameters
  real sigma2_max;
  vector[n_groups] mub;
  vector[n_groups] taub2priormean;
  vector[n_groups] sigmab;
  vector[n_groups] mua;
  vector[n_groups] taua2priormean;
  vector[n_groups] sigmaa;
  real mu0;
  real tau02priormean;
  real sigma0;

  // Dates to use to sample from the predictive
  int<lower=0> T;
  vector[T] tpspred;
}

transformed data {
  vector[N] sqrt_lambda = sqrt(lambda);

  int<lower=1, upper=n_groups> gpred[n_groups];
  int<lower=1, upper=n_groups> idpred[T * n_groups];
  vector[T*n_groups] realtpspred;
  for (i in 1:n_groups) {
    gpred[i] = i;
    idpred[((i-1)*T + 1):(i*T)] = rep_array(i, T);
    realtpspred[((i-1)*T + 1):(i*T)] = tpspred;
  }
}
parameters {
  // Global level
  real z0bar;
  real<lower=0> tau02;
  real<lower=0, upper=sigma2_max> sigma2;

  // Group level
  vector[n_groups] log_bbar;
  vector<lower=0>[n_groups] abar;
  vector<lower=0>[n_groups] taub2;
  vector<lower=0>[n_groups] taua2;

  // Individual level
  vector[n_ind] log_b;
  vector[n_ind] z0;
  vector<lower=0>[n_ind] a;
}

transformed parameters {
  real tau0 = sqrt(tau02);
  real sigma = sqrt(sigma2);
  vector[n_groups] taub = sqrt(taub2);
  vector[n_groups] taua = sqrt(taua2);
}

model {
  vector[N] theo;
  vector[N] residual;

  sigma2 ~ normal(0, sigma2_max);

  // plateaux
  log_bbar ~ normal(mub, sigmab);
  taub2 ~ inv_gamma(2., 1 * taub2priormean);
  log_b ~ normal(log_bbar[g], taub[g]);

  // vitesse
  abar ~ normal(mua, sigmaa);
  taua2 ~ inv_gamma(2., 1. * taua2priormean);
  a ~ normal(abar[g], taua[g]);

  // valeurs initiale
  z0bar ~ normal(mu0, sigma0);
  tau02 ~ inv_gamma(2., 1. * tau02priormean);
  z0 ~ normal(z0bar, tau0);

  // bruit
  theo = path(N, tps, z0[id], a[id], log_b[id]);
  z ~ normal(theo, sigma * sqrt_lambda);
}
generated quantities {
  vector<lower=0>[n_groups] apred;
  vector[n_groups] log_bpred;
  vector[n_groups] z0pred;
  vector[n_groups * T] ypred;
  for (i in 1:n_groups) {
    apred[i] = pos_normal_rng(abar[i], taua[i]);
    log_bpred[i] = pos_normal_rng(log_bbar[i], taub[i]);
    z0pred[i] = pos_normal_rng(z0bar, tau0);
  }

  ypred = path(n_groups*T, realtpspred, z0pred[idpred],
    apred[idpred], log_bpred[idpred]) +
    to_vector(normal_rng(rep_array(0., n_groups*T),
    rep_array(sigma, n_groups*T)));
}
