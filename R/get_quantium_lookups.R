#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
get_quantium_lookups <- function() {

  list(
    age = read_csv(
      "~/not_synced/vaccination/quantium_wp4/dim_age_band.csv",
      col_types = cols(
        age_band_id = col_double(),
        age_band = col_character()
      )
    ),
    product = read_csv(
      "~/not_synced/vaccination/quantium_wp4/dim_vaccine.csv",
      col_types = cols(
        vaccine = col_double(),
        name = col_character(),
        short_name = col_character()
      )
    ),
    date = read_csv(
      "~/not_synced/vaccination/quantium_wp4/dim_time.csv",
      col_types = cols(
        time = col_double(),
        week_starting = col_date("%d/%m/%Y")
      )
    ),
    scenario = read_csv(
      "~/not_synced/vaccination/quantium_wp4/dim_scenario.csv",
      col_types = cols(
        scenario = col_double(),
        `5-11 uptake curve` = col_character(),
        booster_shape = col_character(),
        booster_scale_by_age = col_character(),
        booster_uptake_terminal = col_double(),
        booster_uptake_months = col_double()
      )
    )
  )

}
