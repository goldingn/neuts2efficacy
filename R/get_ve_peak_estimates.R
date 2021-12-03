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
get_ve_peak_estimates <- function() {

  tibble::tribble(
    ~outcome, ~product, ~dose, ~ve, ~ve_lower, ~ve_upper, ~source,
    "acquisition", "AZ", 1, 0.46, 0.35, 0.55, "Pouwels",
    "acquisition", "AZ", 2, 0.67, 0.62, 0.71, "Pouwels",
    "acquisition", "Pfizer", 1, 0.57, 0.50, 0.63, "Pouwels",
    "acquisition", "Pfizer", 2, 0.8, 0.77, 0.83, "Pouwels",
    "transmission", "AZ", 1, 0.05, 0.01, 0.09, "Eyre",
    "transmission", "AZ", 2, 0.24, 0.18, 0.30, "Eyre",
    "transmission", "Pfizer", 1, 0.17, 0.14, 0.19, "Eyre",
    "transmission", "Pfizer", 2, 0.50, 0.35, 0.61, "Eyre",
    "symptoms", "AZ", 1, 0.433, 0.423, 0.442, "Andrews",
    "symptoms", "AZ", 2, 0.652, 0.649, 0.656, "Andrews",
    "symptoms", "Pfizer", 1, 0.519, 0.514, 0.524, "Andrews",
    "symptoms", "Pfizer", 2, 0.835, 0.833, 0.836, "Andrews",
    "hospitalisation", "AZ", 1, 0.814, 0.787, 0.837, "Andrews",
    "hospitalisation", "AZ", 2, 0.930, 0.924, 0.935, "Andrews",
    "hospitalisation", "Pfizer", 1, 0.918, 0.904, 0.930, "Andrews",
    "hospitalisation", "Pfizer", 2, 0.967, 0.963, 0.970, "Andrews",
    "death", "AZ", 1, 0.884, 0.782, 0.938, "Andrews",
    "death", "AZ", 2, 0.927, 0.907, 0.943, "Andrews",
    "death", "Pfizer", 1, 0.886, 0.773, 0.943, "Andrews",
    "death", "Pfizer", 2, 0.952, 0.937, 0.964, "Andrews"
  ) %>%
    mutate(
      days = 0,
      .before = ve
    ) %>%
    mutate(
      days_earliest = days,
      days_latest = days
    )

}
