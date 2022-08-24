#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ngm_with_immunity
#' @return
#' @author Nick Golding
#' @export

# given a square matrix, compute the reproduction number (dominant eigenvalue)
get_reproduction_number <- function(ngm_with_immunity) {

  eigs <- eigen(ngm_with_immunity)
  eigs$values[1]

}
