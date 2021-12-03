#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param neut_model
#' @param draws
#' @return
#' @author Nick Golding
#' @export
predict_ves <- function(neut_model, draws) {

  # uncertain multiplier for booster doses (mean 5-fold, 95% CI 3 to 7-fold)
  booster_multiplier <- normal(5, 1, truncation = c(0, Inf))
  log10_booster_multipler <- log10(booster_multiplier)

  az_idx <- neut_model$lookups$product == "AZ"
  pfizer_idx <- neut_model$lookups$product == "Pfizer"

  # prepare for prediction
  log10_neuts_list <- list(
    az_dose_1 = neut_model$model_objects$dose_1_mean_log10_neuts[az_idx],
    az_dose_2 = neut_model$model_objects$dose_2_mean_log10_neuts[az_idx],
    pfizer_dose_1 = neut_model$model_objects$dose_1_mean_log10_neuts[pfizer_idx],
    pfizer_dose_2 = neut_model$model_objects$dose_2_mean_log10_neuts[pfizer_idx],
    booster = neut_model$model_objects$dose_2_mean_log10_neuts[pfizer_idx] +
      log10_booster_multipler,
    infection = 0
  )

  lookups <- neut_model$lookups
  lookups$immunity <- names(log10_neuts_list)
  peak_mean_log10_neuts_all <- do.call(c, log10_neuts_list)

  # create a table of all the outcomes, immunity types, and days post dose, to compute the VEs
  ve_prediction_data <- expand_grid(
    outcome = lookups$outcome,
    immunity = lookups$immunity,
    days = 0:365
  )

  indices <- list(
    outcome_idx = match(ve_prediction_data$outcome, lookups$outcome),
    immunity_idx = match(ve_prediction_data$immunity, lookups$immunity)
  )

  mean_log10_neuts_all <- log10_neut_over_time(
    time = ve_prediction_data$days,
    maximum_log10_neut = peak_mean_log10_neuts_all[indices$immunity_idx],
    decay = neut_model$model_objects$neut_decay
  )

  # predict VEs
  ve_predict <- ve_from_mean_log10_neut(
    mean_log10_neut_vec = mean_log10_neuts_all,
    c50_vec = neut_model$model_objects$c50s[indices$outcome_idx],
    log_k = neut_model$model_objects$log_k,
    sd_log10_neut = neut_model$model_objects$sd_log10_neut_titres,
    method = "gaussian"
  )

  ve_predict_sims <- calculate(
    ve_predict,
    values = draws,
    nsim = 1000
  )[[1]][, , 1]

  ve_prediction_data %>%
    mutate(
      ve_predict_mean = colMeans(ve_predict_sims),
      ve_predict_sd = apply(ve_predict_sims, 2, stats::sd),
      ve_predict_lower_50 = apply(ve_predict_sims, 2, quantile, 0.25),
      ve_predict_upper_50 = apply(ve_predict_sims, 2, quantile, 0.75),
      ve_predict_lower_90 = apply(ve_predict_sims, 2, quantile, 0.05),
      ve_predict_upper_90 = apply(ve_predict_sims, 2, quantile, 0.95)
    ) %>%
    mutate(
      immunity_type = case_when(
        immunity == "az_dose_2" ~ "AZ vaccine dose 2",
        immunity == "az_dose_1" ~ "AZ vaccine dose 1",
        immunity == "pfizer_dose_2" ~ "Pfizer vaccine dose 2",
        immunity == "pfizer_dose_1" ~ "Pfizer vaccine dose 1",
        immunity == "booster" ~ "mRNA booster",
        immunity == "infection" ~ "Infection"
      )
    )


}
