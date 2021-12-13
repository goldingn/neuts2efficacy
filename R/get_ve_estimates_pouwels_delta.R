#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_estimates_pouwels_delta <- function() {

  tibble::tribble(
    ~variant, ~outcome, ~product, ~dose, ~ve, ~ve_lower, ~ve_upper, ~source,
    "delta", "acquisition", "AZ", 1, 0.46, 0.35, 0.55, "pouwels",
    "delta", "acquisition", "AZ", 2, 0.67, 0.62, 0.71, "pouwels",
    "delta", "acquisition", "Pfizer", 1, 0.57, 0.50, 0.63, "pouwels",
    "delta", "acquisition", "Pfizer", 2, 0.8, 0.77, 0.83, "pouwels"
  ) %>%
    mutate(
      days_earliest = 0,
      # set the upper bound to 8 weeks
      days_latest = 8 * 7
    )

}
