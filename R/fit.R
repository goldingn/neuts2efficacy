#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @return
#' @author Nick Golding
#' @export
fit <- function(neut_model) {

  mcmc(neut_model$greta_model, chains = 10)

}
