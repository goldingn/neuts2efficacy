#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_estimates <- function() {

  bind_rows(
    get_ve_transmission_estimates(),
    get_ve_clinical_waning_estimates()
  )

}
