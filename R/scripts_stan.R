#' Prepare a list of data for Stan
#'
#' @param dat The dataset
#' @param z0.global If FALSE, assume a group hierachical prior on z0
#' @param groups.explain If TRUE, write how the type were coded by numbers.
#'
#' @return A list ready to be passed to stan
set_data <- function(dat, z0.global = TRUE, groups.explain = TRUE) {
  rtn <- list()

  # Code groups with numbers instead of labels
  ll <- levels(factor(dat$type, levels = unique(dat$type)))
  dat$type <- as.integer(factor(dat$type,
                                    levels = unique(dat$type)))
  dat$id <- as.integer(as.factor(dat$id))

  individus <- dat %>%
    group_by(id, type) %>%
    summarise(T = n(), .groups = "drop")

  if (groups.explain) {
    cat("Groups have been coded as follow:\n")
    for (i in 1:length(ll)) {
      cat("Group number", i, "is", ll[i], "\n")
    }
  }

  # Hierarchical data
  rtn$n_groups <- length(ll)
  rtn$n_ind <- nrow(individus)
  rtn$g <- individus$type
  rtn$N <- nrow(dat)
  rtn$id <- dat$id
  rtn$z <- dat$log_y
  rtn$tps <- dat$tps
  rtn$lambda <- dat$lambda

  # Hyperparameters
  rtn$sigma2_max <- 0.3 # 0.4
  rtn$mub <- rep(20, rtn$n_groups)
  rtn$sigmab <- rep(5, rtn$n_groups)
  rtn$taub2priormean <- rep(1.6, rtn$n_groups)^2
  rtn$mua <- rep(0.1, rtn$n_groups)
  rtn$sigmaa <- rep(0.5, rtn$n_groups)
  rtn$taua2priormean <- rep(0.05, rtn$n_groups)^2
  if (z0.global) {
    rtn$mu0 <- 10
    rtn$tau02priormean <- 2
    rtn$sigma0 <- 5
  } else {
    rtn$mu0 <- rep(10, rtn$n_groups)
    rtn$tau02priormean <- rep(2, rtn$n_groups)^2
    rtn$sigma0 <- rep(5, rtn$n_groups)
  }
  # Sample from the predictive at those dates
  rtn$tpspred <- seq(from = 0, to = 12, by = 1)
  rtn$T <- length(rtn$tpspred)

  return(rtn)
}

#' Propose initial values for the MCMC
#'
#' @param par_LS the least square estimates
#' @param data_stan the list, as outputed by set_data
#' @param z0.global If FALSE, assume a group hierachical prior on z0
#'
#' @return A list to be used by stan
init_fun <- function(par_LS, data_stan, z0.global = TRUE,...) {
  rtn <- list()
  g <- par_LS$type <- as.integer(factor(par_LS$type,
                                        levels = unique(par_LS$type)))
  par_LS$id <- as.integer(as.factor(par_LS$id))
  par_LS_mean <- par_LS %>%
    group_by(type) %>%
    summarise(log_bbar = median(log_b),
              abar = median(a),
              z0bar = median(z0),
              taub2 = var(log_b),
              taua2 = var(a),
              tau02 = var(z0),
              .groups = "drop")

  # Global level
  rtn$sigma2 <- data_stan$sigma2_max / 4
  if (z0.global) {
    rtn$z0bar <- median(par_LS$z0)
    rtn$tau02 <- var(par_LS$z0)
  }

  # Group level
  rtn$log_bbar <- par_LS_mean$log_bbar
  rtn$abar <- par_LS_mean$abar
  rtn$taub2 <- par_LS_mean$taub2
  rtn$taua2 <- par_LS_mean$taua2
  if (!z0.global) {
    rtn$z0bar <- par_LS_mean$z0bar
    rtn$tau02 <- par_LS_mean$tau02
  }

  # Individual level
  rtn$log_b <- par_LS$log_b
  rtn$z0 <- par_LS$z0
  rtn$a <- par_LS$a

  return(rtn)
}

mysummary <- function(x, pars){
  rtn <- rstan::summary(x, pars = pars)$summary[, 1:8]
  as_tibble(round(rtn, digits = 3), rownames = "par")
}
