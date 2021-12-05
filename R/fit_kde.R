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
fit_kde <- function(params) {
  MASS::kde2d(
    x = params[[1]],
    y = params[[2]],
    n = 300
    # h = c(3, 15)
  )
}