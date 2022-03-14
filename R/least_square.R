#' Fit Gompertz curves with a least square estimate
#'
#' @param data The data.set, in long format
#' @param lambda A vector to change the weight of a few time points
#' @param parinit The initial parameter value
#'
#' @return The least square estimates for each individual in the dataset
least_square <- function(data, lambda = NULL,
                            parinit = c(20, 10, 0.5)) {
  min.RSS <- function(v, tps, lambda, par) {
    log_b <- par[1]
    z0 <- par[2]
    a <- par[3]
    v_theo <- log_b + exp(-a * tps) * (z0 - log_b)
    return(sum((v - v_theo)^2/lambda))
  }
  par_LS <- data %>% group_by(id, type) %>%
    summarise(log_b = as.numeric(NA),
              z0 = log_b,
              a = log_b,
              .groups = 'drop')
  tps <- data$tps[data$id == 1]
  if (is.null(lambda)) {
    lambda <- rep(1, times = length(tps))
  }
  for (i in par_LS$id) {
    result <- optim(par = parinit, fn = min.RSS,
                    v = data$log_y[data$id == i],
                    tps = tps,
                    lambda = lambda)
    par_LS[i, -(1:2)] <- as.list(result$par)
  }
  return(par_LS)
}
