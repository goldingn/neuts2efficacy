#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param sar_omicron
#' @return
#' @author Nick Golding
#' @export
far_to_sar <- function(sar, coefs = get_far_sar_coefs()) {
  coefs["(Intercept)"] + sar * coefs["sar"]
}
