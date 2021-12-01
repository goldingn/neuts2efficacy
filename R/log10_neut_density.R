#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param x
#' @param mean_log10_neut
#' @param sd_log10_neut
#' @return
#' @author Nick Golding
#' @export
# compute the density of the log10 neut distribution, with greta or R
log10_neut_density <- function(x, mean, sd) {

  with_greta <- inherits(x, "greta_array") |
    inherits(mean, "greta_array") |
    inherits(sd, "greta_array")

  if (with_greta) {
    normal_density(x, mean, sd)
  } else {
    dnorm(x, mean, sd)
  }
}
