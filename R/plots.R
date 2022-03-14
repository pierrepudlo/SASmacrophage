#' Plot the dataset and the fitted Gompertz curves
#'
#' @param dat The dataset
#' @param par_LS The least square estimates
#' @param n Number of points to draw the Gompertz curves
#'
#' @return A ggplot object
plot_LS <- function(dat, par_LS, n = 101) {
  tps.grid <- seq(from = min(dat$tps),
                  to = max(dat$tps),
                  length.out = n)
  ff <- function(id){
    log_b <- par_LS$log_b[par_LS$id == id]
    z0 <- par_LS$z0[par_LS$id == id]
    a <- par_LS$a[par_LS$id == id]
    log_b + exp(-a * tps.grid) * (z0 - log_b)
  }
  par_LS <- par_LS %>%
    mutate(tps = map(id, function(id) tps.grid),
           log_y = map(id, ff)) %>%
    unnest(cols = c(tps, log_y))
  g <- ggplot(dat, aes(x = tps, y = log_y))
  if ("lambda" %in% colnames(dat)) {
    g <- g + geom_point(aes(size = 1/sqrt(lambda))) +
      guides(size = FALSE)
  } else {
    g <- g + geom_point()
  }
  g <- g +
    geom_line() +
    geom_line(data = par_LS, color = "red") +
    scale_x_continuous(breaks = sort(unique(dat$tps))) +
    labs(x = "t", y = "log y(t)") +
    facet_wrap(~ id, scales = "free_y")
}

#' Draw posterior Gompertz curve envelop
#'
#' @param fit A stanfint object
#' @param data The dataset
#' @param id A vector of length 1, containing the id of the individual to plot
#' @param t The time grid
#' @param ci_level Posterior probability of the ribbon of darkest color
#' @param color Color of the ribbon
#'
#' @return A ggplot object
envelop <- function(fit, data, id, t = seq(0, 11, length.out = 101),
                        ci_level = 0.5, color = "#4B8EC7") {
  nt <- length(t)
  pars <- paste0(c("log_b[", "a[", "z0["),
                 id, "]")
  pars <- c(pars, "lp__")
  postdata <- as.data.frame(fit, pars = pars)
  colnames(postdata) <- c("log_b","a", "z0", "lp__")
  env1 <- tibble(
    t = t,
    mean = NA,
    min.95 = NA,
    max.95 = NA,
    min.ci = NA,
    max.ci = NA
  )

  for (i in 1:nt) {
    y <- (postdata$log_b) + ((postdata$z) - (postdata$log_b)) *
      exp(-postdata$a * t[i])
    env1$mean[i] <- mean(y)
    env1$min.95[i] <- quantile(y, probs = 0.025)
    env1$max.95[i] <- quantile(y, probs = 0.975)
    env1$min.ci[i] <- quantile(y, probs = (1 - ci_level)/2)
    env1$max.ci[i] <- quantile(y, probs = 1 - (1 - ci_level)/2)
  }
  ggplot(env1, aes(x = t)) +
    geom_path(aes(y = mean), color = color) +
    geom_ribbon(aes(ymin = min.ci, ymax = max.ci),
                fill = color, color = color, alpha = .5) +
    geom_ribbon(aes(ymin = min.95, ymax = max.95),
                fill = color, color = color, alpha = .1) +
    geom_point(mapping = aes(x = tps, y = log_y, size = 1/lambda),
               data = data[data$id == id, ]) +
    geom_line(mapping = aes(x = tps, y = log_y),
              data = data[data$id == id, ]) +
    guides(size = FALSE) +
    labs(x = "t (days)", y = "log y(t)", title = id)
}

#' Draw the posterior distribution at the group level
#'
#' @param fit A stanfit object
#' @param xlim Range of the x-axis
#' @param ylim Range of the y-axis
#'
#' @return An object of class ggExtraPlot
compare_type <- function(fit, xlim = NULL, ylim = NULL) {
  getLevel <- function(x,y,prob=c(0.95, 0.75, 0.5, 0.25, 0.05)) {
    kk <- MASS::kde2d(x, y)
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    #approx(c1, sz, xout = 1 - prob)$y
    dimnames(kk$z) <- list(kk$x, kk$y)
    dc <- reshape2::melt(kk$z)
    dc$prob <- approx(sz, 1 - c1, dc$value)$y
    approx(c1, sz, xout = 1 - prob)$y
  }
  postdata <- as.data.frame(fit, pars = c("log_bbar", "abar", "z0bar", "lp__"))
  colnames(postdata) <- c("b1", "b2", "a1", "a2", "z0", "lp__")
  n <- nrow(postdata)
  longpostdata <- tibble(
    type = rep(c("AA", "6-OHDA"), each = n),
    log_b = 0,
    a = 0,
    lp__ = 0)
  longpostdata$log_b[1:n] <- postdata$b1
  longpostdata$log_b[n + 1:n] <- postdata$b2
  longpostdata$a[1:n] <- postdata$a1
  longpostdata$a[n + 1:n] <- postdata$a2

  lv1 <- getLevel(postdata$b1, postdata$a1)
  lv2 <- getLevel(postdata$b2, postdata$a2)
  if (is.null(xlim))
    xlim <- range(longpostdata$log_b)
  if (is.null(ylim))
    ylim <- range(longpostdata$a)

  j2 <- ggplot(longpostdata, aes(x = log_b, y = a, colour = type)) +
    stat_density2d(data = longpostdata %>% filter(`type` == "AA"),
                   geom = "polygon",
                   aes(alpha = (..level..)^2, fill = type), colour = 0,
                   breaks = lv1) +
    stat_density2d(data = longpostdata %>% filter(`type` == "6-OHDA"),
                   geom = "polygon",
                   aes(alpha = (..level..)^2, fill = type), colour = 0,
                   breaks = lv2) +
    scale_fill_manual(name = "Group:",
                      values = c("#4B8EC7", "#000000"),
                      breaks = c("6-OHDA", "AA")) +
    scale_colour_manual(name = "Group:",
                        values = c("#4B8EC7", "#000000"),
                        breaks = c("6-OHDA", "AA")) +
    # scale_alpha_continuous(name = "Probability:",
    #                        range = c(0.1, 1),
    #                        labels = c("5%", "25%", "50%", "75%", "95%"),
    #                        breaks = c((lv1^2)) ) +
    geom_point(alpha = 0) +  xlim(xlim) + ylim(ylim)  +
    guides(alpha = FALSE) +
    labs(x = "log(b)", y = "a")
  rtn <- ggMarginal(j2, type = "boxplot",
                    groupColour = TRUE, groupFill = TRUE)
  return(rtn)
}
