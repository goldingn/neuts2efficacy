#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_clinical_waning_estimates <- function() {

  tibble::tribble(
    ~outcome, ~product, ~weeks, ~ve, ~ve_lower, ~ve_upper,
    "symptoms", "AZ", "1", 0.627, 0.617, 0.638,
    "symptoms", "AZ", "2-9", 0.667, 0.663, 0.670,
    "symptoms", "AZ", "10-14", 0.593, 0.588, 0.599,
    "symptoms", "AZ", "15-19", 0.526, 0.517, 0.535,
    "symptoms", "AZ", "20+", 0.473, 0.450, 0.496,
    "symptoms", "Pfizer", "1", 0.924, 0.921, 0.927,
    "symptoms", "Pfizer", "2-9", 0.898, 0.896, 0.900,
    "symptoms", "Pfizer", "10-14", 0.803, 0.799, 0.806,
    "symptoms", "Pfizer", "15-19", 0.734, 0.729, 0.739,
    "symptoms", "Pfizer", "20+", 0.697, 0.687, 0.705,
    "hospitalisation", "AZ", "1", 0.939, 0.913, 0.957,
    "hospitalisation", "AZ", "2-9", 0.952, 0.946, 0.956,
    "hospitalisation", "AZ", "10-14", 0.914, 0.905, 0.922,
    "hospitalisation", "AZ", "15-19", 0.868, 0.851, 0.884,
    "hospitalisation", "AZ", "20+", 0.770, 0.703, 0.823,
    "hospitalisation", "Pfizer", "1", 0.997, 0.967, 0.100,
    "hospitalisation", "Pfizer", "2-9", 0.984, 0.979, 0.988,
    "hospitalisation", "Pfizer", "10-14", 0.965, 0.959, 0.971,
    "hospitalisation", "Pfizer", "15-19", 0.944, 0.934, 0.952,
    "hospitalisation", "Pfizer", "20+", 0.927, 0.903, 0.946,
    "death", "AZ", "2-9", 0.941, 0.897, 0.944,
    "death", "AZ", "10-14", 0.924, 0.897, 0.944,
    "death", "AZ", "15-19", 0.891, 0.842, 0.925,
    "death", "AZ", "20+", 0.787, 0.527, 0.904,
    "death", "Pfizer", "2-9", 0.982, 0.959, 0.992,
    "death", "Pfizer", "10-14", 0.952, 0.930, 0.967,
    "death", "Pfizer", "15-19", 0.939, 0.911, 0.958,
    "death", "Pfizer", "20+", 0.904, 0.851, 0.938
  ) %>%
    mutate(
      # add dose and days info
      dose = 2,
      .after = product
    ) %>%
    mutate(
      source = "Andrews",
      # convert 'weeks' to days post second dose text says week 1 is 7-13 days,
      # so presumably week 0 is 0-6 days. Compute midpoint in days of other
      # periods accordingly. For 20+ weeks, assume the midpoint of 20-30 weeks.
      # Data includes events up to 3rd September 2021, and most vaccinations
      # starting Jan 4 2021 (supp 1), so subtracting 4 weeks for time to
      # immunity, longest post-vaccination period could be is 214 days = ~30
      # weeks
      days_earliest = case_when(
        weeks == "1" ~ 7,
        weeks == "2-9" ~ 7 * 2,
        weeks == "10-14" ~ 7 * 10,
        weeks == "15-19" ~ 7 * 15,
        weeks == "20+" ~ 7 * 20
      ),
      days_latest = case_when(
        weeks == "1" ~ 7 * 2 - 1,
        weeks == "2-9" ~ 7 * 10 - 1,
        weeks == "10-14" ~ 7 * 15 - 1,
        weeks == "15-19" ~ 7 * 20 - 1,
        weeks == "20+" ~ 7 * 31 - 1,
      )
    ) %>%
    rowwise() %>%
    mutate(
      days = mean(days_earliest:days_latest),
      .after = dose
    ) %>%
    # days are expressed as time since second dose, reset to time since peak (14 days post dose) and filter out any before then
    mutate(
      across(
        starts_with("days"),
        ~ .x - 14
      )
    ) %>%
    filter(
      days >= 0
    ) %>%
    ungroup() %>%
    select(
      -weeks
    )

}
