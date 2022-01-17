#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1046853/technical-briefing-34-14-january-2022.pdf
get_ve_estimates_hsa_delta <- function() {

  tibble::tribble(

    ~variant, ~outcome, ~product, ~dose, ~weeks, ~ve, ~ve_lower, ~ve_upper,

    # Figure 11 only include Pfizer boosters:
    "delta", "symptoms", "AZ", 2, "2-4", 0.82, 0.74, 0.89,
    "delta", "symptoms", "AZ", 2, "5-9", 0.76, 0.69, 0.81,
    "delta", "symptoms", "AZ", 2, "10-14", 0.68, 0.65, 0.73,
    "delta", "symptoms", "AZ", 2, "15-19", 0.54, 0.52, 0.56,
    "delta", "symptoms", "AZ", 2, "20-24", 0.46, 0.44, 0.48,
    "delta", "symptoms", "AZ", 2, "25+", 0.42, 0.40, 0.44,

    "delta", "symptoms", "Pfizer", 2, "2-4", 0.92, 0.90, 0.94,
    "delta", "symptoms", "Pfizer", 2, "5-9", 0.85, 0.83, 0.87,
    "delta", "symptoms", "Pfizer", 2, "10-14", 0.79, 0.77, 0.81,
    "delta", "symptoms", "Pfizer", 2, "15-19", 0.76, 0.74, 0.78,
    "delta", "symptoms", "Pfizer", 2, "20-24", 0.67, 0.65, 0.69,
    "delta", "symptoms", "Pfizer", 2, "25+", 0.63, 0.61, 0.65,

    "delta", "symptoms", "Pfizer", 3, "1", 0.93, 0.91, 0.95,
    "delta", "symptoms", "Pfizer", 3, "2-4", 0.96, 0.94, 0.98,
    "delta", "symptoms", "Pfizer", 3, "5-9", 0.93, 0.91, 0.95,
    "delta", "symptoms", "Pfizer", 3, "10+", 0.91, 0.89, 0.93,

  ) %>%
    mutate(
      source = "hsa_delta",
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
