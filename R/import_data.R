#' Transform dataset to a long format
#'
#' @details Assume the 1st column of the input is a type, and all other columns
#' are observed values at different times. Names of thoses columns should be
#' the time, prefixed by "T"
#'
#' @param dat A data.frame
#' @return A data.frame in long format.
long_log_format <- function(dat) {
  n <- nrow(dat)
  dat <- cbind(tibble(id = 1:n), dat)
  p <- ncol(dat)
  dat[, 3:p] <- log(dat[, 3:p])
  dat <- dat %>%
    mutate(log_y0 = T0) %>%
    pivot_longer(cols = 3:p, names_to = "tps", names_prefix = "T",
                 values_to = "log_y") %>%
    mutate(tps = as.double(tps))
  colnames(dat)[2] <- "type"
  return(dat)
}

#' Add a column named lambda to the data.frame
#'
#' @param dat The data.frame to where we will add the lambda column
#' @param change_lambda The data.frame explaining where lambda should not
#' be equal to 1
#'
#' @return The input dat completed with a new column named lambda.
correct_lambda <- function(dat, change_lambda) {
  dat <- dat %>%
    mutate(lambda = 1)
  for (i in 1:nrow(change_lambda)) {
    ligne <- which(dat$id == change_lambda$id[i] &
                     dat$tps == change_lambda$tps[i])
    if (!length(ligne) == 1)
      stop("bug")
    dat$lambda[ligne] <- change_lambda$lambda[i]
  }
  return(dat)
}
