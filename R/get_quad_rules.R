#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @return
#' @author Nick Golding
#' @export
# get quadrature rules (values and weights) for Legendre Gaussian quadrature
get_quad_rules <- function(n_observations, lower = -3, upper = 3) {

  n_quads <- round(5 * (upper - lower))

  # get quadrature rules on (-1, 1)
  quads <- gaussquad::legendre.quadrature.rules(n_quads)[[n_quads]]

  # transform them to (lower, upper) and return
  lambda <- (upper - lower) / 2
  mu <- (lower + upper) / 2

  list(
    values = lambda * quads$x + mu,
    weights = lambda * quads$w
  )

}
