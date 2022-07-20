#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param sim_omicron_params
#' @return
#' @author Nick Golding
#' @export
plot_omicron_params <- function(omicron_params) {

  # plot in LHSTM style
  omicron_kde <- omicron_params %>%
    relocate(
      immune_evasion,
      R0_ratio,
      .before = everything()
    ) %>%
    fit_kde()

  ggplot() +
    geom_polygon(
      aes(x, y),
      data = get_kde_contour(omicron_kde, 0.95),
      fill = lighten("seagreen", 0.6)
    ) +
    geom_polygon(
      aes(x, y),
      data = get_kde_contour(omicron_kde, 0.5),
      fill = lighten("seagreen", 0.1)
    ) +
    scale_x_continuous(
      labels = scales::percent,
      limits = c(0, 1)
    ) +
    geom_hline(
      yintercept = 1,
      linetype = 2
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ylim(0.2, 1.5) +
    ggtitle("Immune escape and relative transmissibility of Omicron vs. Delta") +
    xlab("Omicron immune evasion (reduction in overall VE against transmission, compared to Delta)") +
    ylab("Omicron relative transmissibility (R0 ratio)") +
    theme_minimal()

}
