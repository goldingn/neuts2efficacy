#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ves_now
#' @return
#' @author Nick Golding
#' @export
get_vaccine_transmission_effects <- function(ves, coverage) {

  lookups <- get_quantium_lookups()

  # get a conmat NGM for Australia
  australia_ngm <- get_ngm()

  # combine coverage and VEs to get transmission reduction for each rollout
  # scenario, and omicron scenario
  ves %>%
    # add back in the younger age_groups
    complete(
      scenario,
      omicron_scenario,
      variant,
      outcome,
      age_band = unique(lookups$age$age_band),
      fill = list(
        ve = 0
      )
    ) %>%
    # get the two transmission VEs as columns
    filter(
      outcome %in% c("acquisition", "transmission")
    ) %>%
    pivot_wider(
      names_from = outcome,
      values_from = ve
    ) %>%
    group_by(
      scenario,
      omicron_scenario,
      variant
    ) %>%
    arrange(
      age_band,
      .by_group = TRUE
    ) %>%
    # join on coverages
    left_join(
      coverage,
      by = c("scenario", "age_band")
    ) %>%
    # compute percentage reduction in acquisition and transmission in each age group
    mutate(
      acquisition_multiplier = 1 - acquisition * coverage,
      transmission_multiplier = 1 - transmission * coverage,
    ) %>%
    select(
      -acquisition,
      -transmission,
      -coverage
    ) %>%
    # transform these into matrices of reduction in transmission, matching the NGM
    summarise(
      transmission_reduction_matrix =
        list(
          outer(
            # 'to' groups are on the rows in conmat, and first element in outer is rows,
            # so acquisition first
            acquisition_multiplier,
            transmission_multiplier,
            FUN = "*"
          )
        ),
      .groups = "drop"
    ) %>%
    group_by(
      scenario,
      omicron_scenario,
      variant
    ) %>%
    mutate(
      ngm_unvaccinated = list(australia_ngm),
      ngm_vaccinated = list(ngm_unvaccinated[[1]] * transmission_reduction_matrix[[1]]),
      vaccination_effect_multiplier = get_R(ngm_vaccinated[[1]]) / get_R(ngm_unvaccinated[[1]])
    ) %>%
    select(
      -transmission_reduction_matrix,
      -ngm_unvaccinated,
      -ngm_vaccinated
    )

}
