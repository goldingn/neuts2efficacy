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
# summarise parameters for further modelling
format_parameters <- function(neut_model, draws) {

  c50_params <- list()
  for(i in seq_along(neut_model$lookups$outcome)) {
    name <- paste0("c50_", neut_model$lookups$outcome[i])
    c50_params[[name]] <- neut_model$model_objects$c50s[i]
  }

  neut_params <- list(
    neut_model$model_objects$peak_mean_log10_neuts[1],
    neut_model$model_objects$peak_mean_log10_neuts[2],
    neut_model$model_objects$peak_mean_log10_neuts[3],
    neut_model$model_objects$peak_mean_log10_neuts[4]
  )
  names(neut_params) <- paste0("log10_mean_neut_", neut_model$lookups$immunity)
  neut_params$log10_mean_neut_Pfizer_dose_3 <- neut_params$log10_mean_neut_Pfizer_dose_2 * 5
  neut_params$log10_mean_neut_infection <- neut_params$log10_mean_neut_Pfizer_dose_2 * 0

  other_params <- neut_model$model_objects[c("log_k", "neut_decay")]

  params <- c(c50_params, neut_params, other_params)
  param_names <- as.list(names(params))

  param_sims <- do.call(
    calculate,
    c(
      params,
      list(
        values = draws,
        nsim = 1e4
      )
    )
  )

  # get means
  params <- vapply(param_sims, mean, FUN.VALUE = numeric(1))

  # add on the SD of neut titres
  c(
    params,
    sd_log10_neut_titres = neut_model$model_objects$sd_log10_neut_titres
  )

}
