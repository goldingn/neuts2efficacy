#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param x
#' @param mean
#' @param sd
#' @return
#' @author Nick Golding
#' @export
# tensorflow function for the normal PDF
tf_normal_density <- function(x, mean, sd) {
  dist <- greta:::tfp$distributions$Normal(loc = mean, scale = sd)
  dist$prob(x)
}
