#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param vaccine_cohorts_now
#' @return
#' @author Nick Golding
#' @export
get_coverage <- function(vaccine_cohorts) {

  # get current coverage with any dose in each age band, for each scenario
  coverage <- vaccine_cohorts %>%
    mutate(
      immune = !is.na(immunity)
    ) %>%
    group_by(
      scenario, age_band
    ) %>%
    summarise(
      coverage = weighted.mean(immune, num_people),
      .groups = "drop"
    )

}
