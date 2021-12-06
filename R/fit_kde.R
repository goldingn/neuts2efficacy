#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param params
#' @return
#' @author Nick Golding
#' @export # fit a 2D kernel density estimate to the posterior samples of the
#'   key variant parameters
fit_kde <- function(params, smooth = c(0.25, 0.25)) {

  expansion <- 0.5 * c(diff(range(params[[1]])), diff(range(params[[2]])))
  expand <- c(
    c(-1, 1) * expansion[1],
    c(-1, 1) * expansion[2]
  )

  limits <- c(
    range(params[[1]]) + expand[1:2],
    range(params[[2]]) + expand[3:4]
  )

  kde <- MASS::kde2d(
    x = params[[1]],
    y = params[[2]],
    lims = limits,
    n = 500,
    h = smooth
  )

  # crop the KDE to the observed region
  keep_x <- kde$x >= min(params[[1]]) & kde$x <= max(params[[1]])
  keep_y <- kde$y >= min(params[[2]]) & kde$y <= max(params[[2]])

  kde$x <- kde$x[keep_x]
  kde$y <- kde$y[keep_y]
  kde$z <- kde$z[keep_x, keep_y]

  kde

}

