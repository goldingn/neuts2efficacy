#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param sim_omicron_params
#' @return
#' @author Nick Golding
#' @export
plot_omicron_params <- function(omicron_params, neut_fold_lines = c()) {

  # plot in LHSTM style
  omicron_kde <- omicron_params %>%
    relocate(
      immune_evasion,
      R0_ratio,
      .before = everything()
    ) %>%
    fit_kde()

  plot <- ggplot() +
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
      limits = c(0, 1),
      labels = scales::percent
    ) +
    geom_hline(
      yintercept = 1,
      linetype = 2
    ) +
    geom_vline(
      xintercept = 0,
      linetype = 2
    ) +
    ggtitle("Immune escape and relative transmissibility of Omicron vs. Delta") +
    xlab("Omicron immune evasion (reduction in overall VE against transmission, compared to Delta)") +
    ylab("Omicron relative transmissibility (R0 ratio)") +
    theme_minimal()


  if (length(neut_fold_lines > 0)) {


    labels <- paste0(neut_fold_lines, "x")
    evasions <- NA * neut_fold_lines

    for (i in seq_along(neut_fold_lines)) {

      fold <- neut_fold_lines[i]
      evasions[i] <- omicron_params %>%
        filter(
          titre_fold > log10(fold - 0.2) & titre_fold < log10(fold + 0.2)
        ) %>%
        pull(
          immune_evasion
        ) %>%
        mean()

    }

    line_df <- tibble(
      fold = labels,
      immune_evasion = evasions,
    )

    plot <- plot +
      geom_vline(
        aes(xintercept = immune_evasion),
        data = line_df,
        alpha = 0.2
      ) +
      geom_text(
        aes(
          x = immune_evasion,
          y = 0.4,
          label = fold
        ),
        data = line_df,
        angle = 90,
        vjust = -0.2,
        hjust = 0,
        size = 3,
        alpha = 0.4
      )

  }

  plot

}
