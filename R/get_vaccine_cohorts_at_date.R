#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param vaccine_scenarios
#' @param nameme1
#' @return
#' @author Nick Golding
#' @export
get_vaccine_cohorts_at_date <- function(vaccine_scenarios, target_date) {

  # set future times to NA and collapse to get contemporary data
  vaccine_scenarios %>%
    mutate(
      across(
        starts_with("date"),
        ~if_else(.x > target_date, as.Date(NA), .x)
      )
    ) %>%
    group_by(
      across(
        -num_people
      )
    ) %>%
    summarise(
      num_people = sum(num_people),
      .groups = "drop"
    ) %>%
    # compute most recent vaccines and how long ago they were for each cohort
    mutate(
      most_recent_dose = pmax(date_dose_1, date_dose_2, date_booster, na.rm = TRUE)
    ) %>%
    pivot_longer(
      cols = starts_with("date"),
      names_to = "dose",
      values_to = "date",
      names_prefix = "date_"
    ) %>%
    # keep only one of these entries
    mutate(
      keep = case_when(
        is.na(date) ~ dose == "dose_1",
        TRUE ~ date == most_recent_dose
      )
    ) %>%
    filter(
      keep
    ) %>%
    select(
      -most_recent_dose,
      -keep
    ) %>%
    mutate(
      dose = if_else(is.na(date), NA_character_, dose),
      product = case_when(
        dose == "booster" ~ vaccine_booster,
        TRUE ~ vaccine
      ),
      # rename products to match VE model, recoding Moderna as Pfizer for now
      # since there is not enough evidence on efficacy to distringuish it from
      # Pfizer
      product = case_when(
        product == "Pfizer" ~ "Pf",
        product == "AstraZeneca" ~ "AZ",
        product == "Moderna" ~ "Pf",
        product == "Booster" ~ "mRNA"
      ),
      immunity = case_when(
        !is.na(date) ~ paste(product, dose, sep = "_"),
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(
      days_ago = as.numeric(target_date - date)
    ) %>%
    select(
      -starts_with("vaccine"),
      -product,
      -dose,
      -date
    )

}
