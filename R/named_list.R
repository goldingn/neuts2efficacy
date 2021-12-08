#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param object
#' @param name
#' @return
#' @author Nick Golding
#' @export
named_list <- function(object, name) {
  x <- list(object)
  names(x) <- name
  x
}
