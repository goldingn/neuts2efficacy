#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
# load omicron scenario parameters and rearrange to wide format
get_omicron_params_wide <- function() {
  read_csv(
    "outputs/scenario_parameters_omicron.csv",
    col_types = cols(
      parameter = col_character(),
      intermediate = col_double(),
      optimistic = col_double(),
      pessimistic = col_double()
    )
  ) %>%
    pivot_longer(
      cols = c(-parameter),
      names_to = "omicron_scenario",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = parameter,
      values_from = value
    )
}
