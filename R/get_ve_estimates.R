#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# list VE against symptomatic disease for different products and doses
# from national plan parameters doc
get_ve_estimates <- function() {

  ve_estimates <- tibble::tribble(
    ~outcome, ~product, ~dose, ~ve,
    "acquisition", "AZ", 1, 0.46,
    "acquisition", "AZ", 2, 0.67,
    "acquisition", "Pfizer", 1, 0.57,
    "acquisition", "Pfizer", 2, 0.8,
    "transmission", "AZ", 1, 0.05,  # 0.02
    "transmission", "AZ", 2, 0.24,  # 0.36
    "transmission", "Pfizer", 1, 0.17, # 0.13
    "transmission", "Pfizer", 2, 0.50, # 0.65
    "symptoms", "AZ", 1, 0.4,
    "symptoms", "AZ", 2, 0.71,
    "symptoms", "Pfizer", 1, 0.58,
    "symptoms", "Pfizer", 2, 0.84,
    "hospitalisation", "AZ", 1, 0.81,
    "hospitalisation", "AZ", 2, 0.93,
    "hospitalisation", "Pfizer", 1, 0.92,
    "hospitalisation", "Pfizer", 2, 0.97,
    "death", "AZ", 1, 0.88,
    "death", "AZ", 2, 0.93,
    "death", "Pfizer", 1, 0.89,
    "death", "Pfizer", 2, 0.95
  )

}
