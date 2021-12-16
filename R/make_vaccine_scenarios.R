#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
make_vaccine_scenarios <- function() {

  # get quantium lookup tables
  lookups <- get_quantium_lookups()

  # load all vaccination scenario data (v. large)
  vaccines <- read_csv(
    "~/not_synced/vaccination/quantium_wp4/vaccines.csv",
    col_types = cols(
      .default = col_double()
    )
  )

  # collapse over space
  vaccine_all <- vaccines %>%
    group_by(
      across(c(-sa4_code16, -num_people))
    ) %>%
    summarise(
      num_people = sum(num_people)
    )

  # add on times, age bands, and products, and return
  vaccine_all %>%
    # join on the dates
    ungroup() %>%
    left_join(
      lookups$date,
      by = c(
        "time_dose_1" = "time"
      )
    ) %>%
    rename(
      date_dose_1 = week_starting
    ) %>%
    left_join(
      lookups$date,
      by = c(
        "time_dose_2" = "time"
      )
    ) %>%
    rename(
      date_dose_2 = week_starting
    ) %>%
    left_join(
      lookups$date,
      by = c(
        "time_booster" = "time"
      )
    ) %>%
    rename(
      date_booster = week_starting
    ) %>%
    select(
      -starts_with("time")
    ) %>%
    # join on the products
    left_join(
      lookups$product,
      by = "vaccine"
    ) %>%
    select(
      -vaccine,
      -short_name
    ) %>%
    rename(
      vaccine = name
    ) %>%
    left_join(
      lookups$product,
      by = c("vaccine_booster" = "vaccine")
    ) %>%
    select(
      -vaccine_booster,
      -short_name
    ) %>%
    rename(
      vaccine_booster = name
    ) %>%
    # join on the age bands and convert to a factor
    left_join(
      lookups$age,
      by = "age_band_id"
    ) %>%
    select(
      -age_band_id
    ) %>%
    mutate(
      age_band = factor(
        age_band,
        levels = sort_age_groups(unique(age_band))
      )
    )

}
