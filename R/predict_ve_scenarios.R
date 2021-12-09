#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param scenarios
#' @param neut_model
#' @param draws
#' @param tolerance
#' @return
#' @author Nick Golding
#' @export

# predict a single deterministic set of VEs for a each of a set of scenarios of
# R0 and immune evasion, by finding the posterior means of parameter sets that
# yield R0 and immune evasion in the vicinity of those scenario values
predict_ve_scenarios <- function(scenarios, neut_model, draws, omicron = FALSE, tolerance = 0.05) {

  greta_arrays <- c(
    neut_model$model_objects,
    list(za_immune_evasion = immune_evasion(neut_model))
  )

  sims <- with(
    greta_arrays,
    calculate(
      R0_ratio,
      za_immune_evasion,
      omicron_log10_neut_fold,
      za_baseline_immunity_log10_neut_fold,
      za_distancing_effect,
      za_waning_days,
      reinfection_correction,
      fraction_immune,
      c50s,
      log_k,
      neut_decay,
      dose_1_mean_log10_neuts,
      dose_2_mean_log10_neuts,
      values = draws,
      nsim = 1e4
    )
  )

  params <- sims %>%
    as_tibble() %>%
    mutate(
      booster_multiplier = array(5, dim = c(n(), 1, 1))
    ) %>%
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
        ~list(t(colMeans(.x)[, 1]))
      )
    ) %>%
    # need to turn these into named lists now
    pivot_longer(
      cols = -scenario,
      names_to = "parameter",
      values_to = "value"
    ) %>%
    rowwise() %>%
    mutate(
      param = named_list(value, parameter)
    ) %>%
    ungroup() %>%
    group_by(
      scenario
    ) %>%
    summarise(
      param = list(param),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = scenario,
      values_from = param
    ) %>%
    as.list() %>%
    lapply(`[[`, 1)

  # predict for each scenario and combine
  preds <- list()
  for (i in seq_along(params)) {
    preds[[i]] <- predict_ves(
      neut_model = neut_model,
      draws = params[[i]],
      omicron = omicron,
      nsim = 3
    )
  }
  names(preds) <- names(params)

  pred_df <- do.call(
    bind_rows,
    c(preds, list(.id = "scenario"))
  )

  # format the output
  pred_df %>%
    select(
      -ve_predict_sd,
      -contains("upper"),
      -contains("lower")
    ) %>%
    rename(
      ve = ve_predict_mean
    ) %>%
    relocate(
      ve,
      .after = everything()
    )

}
