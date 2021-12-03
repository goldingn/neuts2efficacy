#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ve_data_modelling
#' @param ve_data_modelling
#' @param ve_fitted_sims
#' @return
#' @author Nick Golding
#' @export
plot_ve_fit <- function(ve_data_modelling, ve_fitted_sims) {

  ve_fitted <- colMeans(ve_fitted_sims)
  ve_fitted_quants <- apply(ve_fitted_sims, 2, quantile, c(0.025, 0.975))

  ve_data_validation <- ve_data_modelling %>%
    mutate(
      predicted_ve = ve_fitted,
      predicted_ve_lower = ve_fitted_quants[1, ],
      predicted_ve_upper = ve_fitted_quants[2, ],
      .after = ve
    ) %>%
    mutate(
      difference = predicted_ve - ve,
      percent_difference = 100 * (1 - predicted_ve / ve)
    )

  ve_data_validation %>%
    ggplot(
      aes(
        x = ve,
        y = predicted_ve,
        colour = days
      )
    ) +
    geom_point() +
    geom_errorbar(
      aes(
        ymin = predicted_ve_lower,
        ymax = predicted_ve_upper
      )
    ) +
    geom_errorbar(
      aes(
        xmin = ve_lower,
        xmax = ve_upper
      )
    ) +
    geom_abline(
      intercept = 0,
      slope = 1
    ) +
    xlab("Published VE estimate") +
    ylab("Predicted VE") +
    theme_minimal()

}
