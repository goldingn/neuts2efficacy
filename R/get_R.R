#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @return
#' @author Nick Golding
#' @export
get_R <- function(matrix) {
  Re(eigen(matrix)$values[1])
}
