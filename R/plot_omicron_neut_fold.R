#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param omicron_params
#' @return
#' @author Nick Golding
#' @export
plot_omicron_neut_fold <- function(omicron_params) {

  omicron_params %>%
    ggplot(
      aes(
        x = exp(titre_fold)
      )
    ) +
    geom_density(
      fill = lighten("seagreen", 0.6)
    ) +
    coord_cartesian(
      xlim = c(1, 10)
    ) +
    xlab("Fold increase in antibodies required to neutralise virus") +
    ylab("Posterior density") +
    ggtitle(
      "Posterior distribution over predicted neutralisation reduction for Omicron",
      "Predicted multiplicative increase in antibodies required to neutralise Omicron\nin vitro, relative to Delta"
    ) +
    theme_minimal()
}
