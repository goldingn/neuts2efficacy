#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# load and prepare vaccine scenarios
get_vaccine_scenarios <- function() {

  # restore if already created
  if (file.exists("outputs/vaccine_scenarios.RDS")) {
    vaccine_scenarios <- readRDS("outputs/vaccine_scenarios.RDS")
  } else {
    vaccine_scenarios <- make_vaccine_scenarios()
    saveRDS(vaccine_scenarios, "outputs/vaccine_scenarios.RDS")
  }

  vaccine_scenarios

}
