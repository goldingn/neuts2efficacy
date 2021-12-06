#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# list VEs based on national plan parameters doc - but with Andrews estiamtes
# for symptomatic disease
get_ve_transmission_estimates <- function() {

  tibble::tribble(
    ~outcome, ~product, ~dose, ~ve, ~ve_lower, ~ve_upper, ~source,
    "acquisition", "AZ", 1, 0.46, 0.35, 0.55, "Pouwels",
    "acquisition", "AZ", 2, 0.67, 0.62, 0.71, "Pouwels",
    "acquisition", "Pfizer", 1, 0.57, 0.50, 0.63, "Pouwels",
    "acquisition", "Pfizer", 2, 0.8, 0.77, 0.83, "Pouwels",
    "transmission", "AZ", 1, 0.05, 0.01, 0.09, "Eyre",
    "transmission", "AZ", 2, 0.24, 0.18, 0.30, "Eyre",
    "transmission", "Pfizer", 1, 0.17, 0.14, 0.19, "Eyre",
    "transmission", "Pfizer", 2, 0.50, 0.35, 0.61, "Eyre"
  ) %>%
    mutate(
      days_earliest = 0,
      # for Eyre, set the upper bound to twice the median, for Pouwels set it to 8 weeks
      days_latest = case_when(
        outcome == "transmission" & product == "AZ" ~ 2 * 51,
        outcome == "transmission" & product == "Pfizer" ~ 2 * 90,
        TRUE ~ 8 * 7
      )
    ) %>%
    rowwise() %>%
    mutate(
      days = mean(days_earliest:days_latest),
      .after = dose
    ) %>%
    ungroup()

}
