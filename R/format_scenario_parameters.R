#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @param draws
#' @param scenarios
#' @return
#' @author Nick Golding
#' @export
# get parameters for IBM in for these scenarios
format_scenario_parameters <- function(neut_model, draws, scenarios, tolerance = 0.05) {

  # simulate parameter values, including the immune evasion parameter
  sims <- format_parameters(
    neut_model = neut_model,
    draws = draws,
    summarise = FALSE,
    include_immune_evasion = TRUE
  )

  sims_tibble <- sims %>%
    lapply(c) %>%
    as_tibble()

  # simulate them all, join to scenarios, subset to those matching the scenarios
  # and summarise
  scenario_params <- sims_tibble %>%
    full_join(
      scenarios,
      by = character()
    ) %>%
    filter(
      near(R0_ratio, R0_ratio_target, tol = tolerance),
      near(za_immune_evasion, immune_evasion_target, tol = tolerance)
    ) %>%
    select(
      -ends_with("target"),
      -za_immune_evasion
    ) %>%
    group_by(
      scenario
    ) %>%
    summarise(
      across(
        everything(),
        mean
      )
    ) %>%
    # need to turn these into named lists now
    pivot_longer(
      cols = -scenario,
      names_to = "parameter",
      values_to = "value"
    ) %>%
    pivot_wider(
      names_from = scenario,
      values_from = value
    )

  # calculate posterior means and add those on too
  means <- sims_tibble %>%
    select(
      -za_immune_evasion
    ) %>%
    summarise(
      across(
        everything(),
        mean
      )
    ) %>%
    # need to turn these into named lists now
    pivot_longer(
      cols = everything(),
      names_to = "parameter",
      values_to = "estimate"
    )

  # return both
  scenario_params %>%
    full_join(
      means,
      by = "parameter"
    )

}
