#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param kd
#' @param nameme1
#' @return
#' @author Nick Golding
#' @export
get_kde_contour <- function(kde, level = 0.95) {

  vals <- sort(kde$z)
  densities <- cumsum(vals) / sum(sum(kde$z))
  cutoff <- vals[which(densities >= (1 - level))[1]]

  contour <- contourLines(
    x = kde$x,
    y = kde$y,
    z = kde$z,
    levels = cutoff
  )

  as_tibble(contour[[1]])

}
