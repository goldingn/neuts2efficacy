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
predict_ves <- function(neut_model, draws, omicron = FALSE, nsim = 1000, omicron_infection_multiplier = 1) {

  log10_neuts_list <- list(
    neut_model$model_objects$peak_mean_log10_neuts[1],
    neut_model$model_objects$peak_mean_log10_neuts[2],
    neut_model$model_objects$peak_mean_log10_neuts[3],
    neut_model$model_objects$peak_mean_log10_neuts[4],
    neut_model$model_objects$peak_mean_log10_neuts[5]
  )
  names(log10_neuts_list) <- neut_model$lookups$immunity

  # add WT infection - just 0 because this is the baseline
  log10_neuts_list$wt_infection <- log10_neuts_list$Pfizer_dose_2 * 0

  # increase neuts when paired with Omicron infection
  for(old_name in names(log10_neuts_list)) {

    new_name <- paste0(old_name, "_plus_omicron_infection")
    old_neut <- log10_neuts_list[[old_name]]
    new_neut <- old_neut + log10(omicron_infection_multiplier)
    log10_neuts_list[[new_name]] <- new_neut

  }

  # prepare for prediction

  # if we are predicting for omicron, add the omicron adjustment to the log10
  # neuts, and add an omicron-only infection
  if (omicron) {

    omicron_adjustment <- neut_model$model_objects$omicron_log10_neut_fold

    log10_neuts_list <- lapply(
      log10_neuts_list,
      "+",
      omicron_adjustment
    )

    # add omicron infection protection against omicron (same as WT infection against WT)
    log10_neuts_list$omicron_infection <- log10_neuts_list$Pfizer_dose_2 * 0

  }

  # overwrite the lookup, with these new ones added on
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

  ve_predict_sims <- with(
    neut_model$model_objects,
    calculate(
      ve_predict,
      values = draws,
      nsim = nsim
    )[[1]][, , 1]
  )

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
        immunity == "AZ_dose_2" ~ "AZ vaccine dose 2",
        immunity == "AZ_dose_1" ~ "AZ vaccine dose 1",
        immunity == "Pfizer_dose_2" ~ "Pfizer vaccine dose 2",
        immunity == "Pfizer_dose_1" ~ "Pfizer vaccine dose 1",
        immunity == "mRNA_booster" ~ "mRNA booster",
        immunity == "wt_infection" ~ "WT Infection",
        immunity == "omicron_infection" ~ "Omicron Infection",
        immunity == "wt_infection_plus_omicron_infection" ~ "WT Infection + Omicron infection",
        immunity == "AZ_dose_2_plus_omicron_infection" ~ "AZ vaccine dose 2 + Omicron infection",
        immunity == "AZ_dose_1_plus_omicron_infection" ~ "AZ vaccine dose 1 + Omicron_infection",
        immunity == "Pfizer_dose_2_plus_omicron_infection" ~ "Pfizer vaccine dose 2 + Omicron infection",
        immunity == "Pfizer_dose_1_plus_omicron_infection" ~ "Pfizer vaccine dose 1 + Omicron infection",
        immunity == "mRNA_booster_plus_omicron_infection" ~ "mRNA booster + Omicron infection",
      )
    )

}
