#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_estimates_eyre_delta <- function() {

  tibble::tribble(
    ~variant, ~outcome, ~product, ~dose, ~ve, ~ve_lower, ~ve_upper, ~source,
    "delta", "transmission", "AZ", 1, 0.05, 0.01, 0.09, "eyre",
    "delta", "transmission", "AZ", 2, 0.24, 0.18, 0.30, "eyre",
    "delta", "transmission", "Pfizer", 1, 0.17, 0.14, 0.19, "eyre",
    "delta", "transmission", "Pfizer", 2, 0.50, 0.35, 0.61, "eyre"
  ) %>%
    mutate(
      days_earliest = 0,
      # for Eyre, set the upper bound to twice the median, for Pouwels set it to 8 weeks
      days_latest = case_when(
        product == "AZ" ~ 2 * 51,
        product == "Pfizer" ~ 2 * 90,
        TRUE ~ 8 * 7
      )
    )

}
