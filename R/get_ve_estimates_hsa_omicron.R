#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1046853/technical-briefing-34-14-january-2022.pdf
get_ve_estimates_hsa_omicron <- function() {

  tibble::tribble(

    ~variant, ~outcome, ~product, ~dose, ~weeks, ~ve, ~ve_lower, ~ve_upper,

    # table 2:
    "omicron", "hospitalisation", "Pfizer", 3, "2-4", 0.92, 0.89, 0.94,
    "omicron", "hospitalisation", "Pfizer", 3, "5-9", 0.88, 0.84, 0.91,
    "omicron", "hospitalisation", "Pfizer", 3, "10+", 0.83, 0.78, 0.87,

    # table 3 - SIREN study (boosters started 16 Sept, status recorded 30 Nov, data reported Jan 16:
    "omicron", "acquisition", "Pfizer", 3, "6-12", 0.71, 0.56, 0.82,

    # Figure 11 only include Pfizer boosters:
    "omicron", "symptoms", "AZ", 2, "2-4", 0.43, 0.30, 0.55,
    "omicron", "symptoms", "AZ", 2, "5-9", 0.32, 0.21, 0.41,
    "omicron", "symptoms", "AZ", 2, "10-14", 0.27, 0.16, 0.35,
    "omicron", "symptoms", "AZ", 2, "15-19", 0.18, 0.12, 0.22,
    "omicron", "symptoms", "AZ", 2, "20-24", 0.05, 0.03, 0.07,
    "omicron", "symptoms", "AZ", 2, "25+", 0, 0, 0.01,

    "omicron", "symptoms", "Pfizer", 2, "2-4", 0.64, 0.61, 0.67,
    "omicron", "symptoms", "Pfizer", 2, "5-9", 0.47, 0.44, 0.50,
    "omicron", "symptoms", "Pfizer", 2, "10-14", 0.29, 0.27, 0.31,
    "omicron", "symptoms", "Pfizer", 2, "15-19", 0.18, 0.16, 0.20,
    "omicron", "symptoms", "Pfizer", 2, "20-24", 0.13, 0.11, 0.15,
    "omicron", "symptoms", "Pfizer", 2, "25+", 0.14, 0.11, 0.17,

    # "omicron", "symptoms", "AZ", 3, "1", 0.58, 0.56, 0.60,
    # "omicron", "symptoms", "AZ", 3, "2-4", 0.62, 0.60, 0.64,
    # "omicron", "symptoms", "AZ", 3, "5-9", 0.55, 0.53, 0.57,
    # "omicron", "symptoms", "AZ", 3, "10+", 0.48, 0.46, 0.50,

    "omicron", "symptoms", "Pfizer", 3, "1", 0.66, 0.64, 0.68,
    "omicron", "symptoms", "Pfizer", 3, "2-4", 0.68, 0.66, 0.70,
    "omicron", "symptoms", "Pfizer", 3, "5-9", 0.58, 0.56, 0.60,
    "omicron", "symptoms", "Pfizer", 3, "10+", 0.55, 0.53, 0.57,

  ) %>%
    mutate(
      source = "hsa_omicron",
      # convert 'weeks' to days post second dose text says week 1 is 7-13 days,
      # so presumably week 0 is 0-6 days. Compute midpoint in days of other
      # periods accordingly. For 25+ weeks, assume the midpoint of 25-30 weeks.
      days_earliest = case_when(
        weeks == "2-4" ~ 7 * 2,
        weeks == "5-9" ~ 7 * 5,
        weeks == "10-14" ~ 7 * 10,
        weeks == "15-19" ~ 7 * 15,
        weeks == "20-24" ~ 7 * 20,
        weeks == "10+" ~ 7 * 10,
        weeks == "25+" ~ 7 * 25,
      ),
      days_latest = case_when(
        weeks == "2-4" ~ 7 * (4 + 1) - 1,
        weeks == "5-9" ~ 7 * (9 + 1) - 1,
        weeks == "10-14" ~ 7 * (14 + 1) - 1,
        weeks == "15-19" ~ 7 * (19 + 1) - 1,
        weeks == "20-24" ~ 7 * (24 + 1) - 1,
        weeks == "25+" ~ 7 * (30 + 1) - 1,
        # UK boosters started Sept 16, so max 4 months
        dose == 3 & weeks == "10+" ~ 7 * 16
      )
    ) %>%
    # days are expressed as time since second dose, reset to time since peak (14 days post dose) and filter out any before then
    mutate(
      across(
        starts_with("days"),
        ~ .x - 14
      )
    ) %>%
    filter(
      days_earliest >= 0
    ) %>%
    select(
      -weeks
    )

}
