#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_ve_estimates_andrews_omicron <- function() {

  tibble::tribble(

    ~variant, ~outcome, ~product, ~dose, ~weeks, ~ve, ~ve_lower, ~ve_upper, ~cases,

    "omicron", "symptoms", "AZ", 2, "15-19", -0.547, -1.74, 0.126, 17,
    "omicron", "symptoms", "AZ", 2, "20-24", -0.132, -0.602, 0.201, 76,
    "omicron", "symptoms", "AZ", 2, "25+", 0.059, -0.297, 0.317, 96,
    "omicron", "symptoms", "AZ", 3, "1-2", 0.719, 0.091, 0.913, 3,
    "omicron", "symptoms", "AZ", 3, "2+", 0.714, 0.418, 0.860, 10,

    "omicron", "symptoms", "Pfizer", 1, "4+", 0.342, -0.035, 0.581, 28,
    "omicron", "symptoms", "Pfizer", 2, "2-9", 0.88, 0.659, 0.958, 4,
    "omicron", "symptoms", "Pfizer", 2, "10-14", 0.485, 0.243, 0.650, 39,
    "omicron", "symptoms", "Pfizer", 2, "15-19", 0.341, 0.097, 0.520, 83,
    "omicron", "symptoms", "Pfizer", 2, "20-24", 0.366, 0.004, 0.596, 27,
    "omicron", "symptoms", "Pfizer", 2, "25+", 0.342, -0.050, 0.587, 25,
    "omicron", "symptoms", "Pfizer", 3, "2+", 0.755, 0.561, 0.863, 16,

    "delta", "symptoms", "AZ", 1, "4+", 0.454, 0.389, 0.511, 553,
    "delta", "symptoms", "AZ", 2, "2-9", 0.762, 0.637, 0.844, 29,
    "delta", "symptoms", "AZ", 2, "10-14", 0.649, 0.552, 0.724, 97,
    "delta", "symptoms", "AZ", 2, "15-19", 0.485, 0.447, 0.520, 1751,
    "delta", "symptoms", "AZ", 2, "20-24", 0.454, 0.430, 0.467, 10728,
    "delta", "symptoms", "AZ", 2, "25+", 0.418, 0.394, 0.441, 13376,
    "delta", "symptoms", "AZ", 3, "1-2", 0.87, 0.855, 0.884, 430,
    "delta", "symptoms", "AZ", 3, "2+", 0.938, 0.932, 0.943, 669,

    "delta", "symptoms", "Pfizer", 1, "4+", 0.361, 0.323, 0.398, 2715,
    "delta", "symptoms", "Pfizer", 2, "2-9", 0.882, 0.867, 0.895, 336,
    "delta", "symptoms", "Pfizer", 2, "10-14", 0.777, 0.763, 0.790, 1818,
    "delta", "symptoms", "Pfizer", 2, "15-19", 0.722, 0.710, 0.734, 4746,
    "delta", "symptoms", "Pfizer", 2, "20-24", 0.648, 0.626, 0.669, 1877,
    "delta", "symptoms", "Pfizer", 2, "25+", 0.635, 0.614, 0.655, 2528,
    "delta", "symptoms", "Pfizer", 3, "1-2", 0.922, 0.907, 0.934, 142,
    "delta", "symptoms", "Pfizer", 3, "2+", 0.926, 0.920, 0.931, 1135

  ) %>%
    mutate(
      source = "andrews_omicron",
      # convert 'weeks' to days post second dose text says week 1 is 7-13 days,
      # so presumably week 0 is 0-6 days. Compute midpoint in days of other
      # periods accordingly. For 20+ weeks, assume the midpoint of 20-30 weeks.
      # Data includes events up to 3rd September 2021, and most vaccinations
      # starting Jan 4 2021 (supp 1), so subtracting 4 weeks for time to
      # immunity, longest post-vaccination period could be is 214 days = ~30
      # weeks
      days_earliest = case_when(
        weeks == "1-2" ~ 7,
        weeks == "2-9" ~ 7 * 2,
        weeks == "10-14" ~ 7 * 10,
        weeks == "15-19" ~ 7 * 15,
        weeks == "20-24" ~ 7 * 20,
        weeks == "25+" ~ 7 * 25,
        weeks == "2+" ~ 7 * 2,
        weeks == "4+" ~ 7 * 4
      ),
      days_latest = case_when(
        weeks == "1-2" ~ 7 * 2 - 1,
        weeks == "2-9" ~ 7 * 10 - 1,
        weeks == "10-14" ~ 7 * 15 - 1,
        weeks == "15-19" ~ 7 * 20 - 1,
        weeks == "20-24" ~ 7 * 25 - 1,
        weeks == "25+" ~ 7 * 31 - 1,
        weeks == "4+" ~ 7 * 31 - 1,
        # UK boosters started Sept 16, so max 3 months
        dose == 3 & weeks == "2+" ~ 7 * 13 - 1
      )
    ) %>%
    # days are expressed as time since second dose, reset to time since peak (14 days post dose) and filter out any before then
    mutate(
      across(
        starts_with("days"),
        ~ .x - 14
      )
    ) %>%
    # subset to those with at least 10 cases to compute the VE
    filter(
      cases > 10
    ) %>%
    select(
      -cases,
      -weeks
    )

}
