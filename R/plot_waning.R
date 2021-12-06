#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param ve_predictions
#' @param ve_data_modelling
#' @return
#' @author Nick Golding
#' @export
plot_waning <- function(ve_predictions,
                        ve_data_modelling = NULL,
                        immunity_levels = c(
                          "mRNA booster",
                          "Pfizer vaccine dose 2",
                          "Pfizer vaccine dose 1",
                          "AZ vaccine dose 2",
                          "AZ vaccine dose 1",
                          "Infection"
                        )
                        ) {

  immunity_colours <- c(
    "mRNA booster" = lighten("darkorchid4", 0.1),
    "Pfizer vaccine dose 2" = lighten("darkorchid1", 0.1),
    "Infection" = grey(0.9),
    "AZ vaccine dose 2" = lighten("firebrick1", 0.1),
    "Pfizer vaccine dose 1" = lighten("darkorchid1", 0.8),
    "AZ vaccine dose 1" = lighten("firebrick1", 0.8)
  )

  immunity_colours <- immunity_colours[immunity_levels]

  plot <- ve_predictions %>%
    filter(
      immunity_type %in% immunity_levels
    ) %>%
    mutate(
      outcome = str_to_sentence(outcome),
      outcome = factor(
        outcome,
        levels = c(
          "Death",
          "Hospitalisation",
          "Symptoms",
          "Acquisition",
          "Transmission"
        )
      ),
      immunity_type = factor(
        immunity_type,
        levels = immunity_levels
      )
    ) %>%
    rename(
      `Type of immunity` = immunity_type
    ) %>%
    ggplot(
      aes(
        x = days,
        fill = `Type of immunity`
      )
    ) +
    facet_wrap(~outcome) +
    geom_ribbon(
      aes(
        ymin = ve_predict_lower_50,
        ymax = ve_predict_upper_50
      ),
      size = 0.25,
      alpha = 0.8,
      colour = grey(0.6)
    ) +
    scale_y_continuous(
      labels = scales::percent
    ) +
    scale_fill_manual(
      values = immunity_colours
    ) +
    scale_alpha_manual(
      values = c("two doses" = 0.8, "one dose" = 0.1)
    ) +
    coord_cartesian(
      xlim = c(0, 200)
    ) +
    ylab("Efficacy") +
    xlab("Days since peak immunity") +
    theme_minimal()


  plot_data <- !is.null(ve_data_modelling)

  if (plot_data) {

    ve_data_plotting <- ve_data_modelling %>%
      mutate(
        `Type of immunity` = case_when(
          product == "AZ" & dose == 2 ~ "AZ vaccine dose 2",
          product == "AZ" & dose == 1 ~ "AZ vaccine dose 1",
          product == "Pfizer" & dose == 2 ~ "Pfizer vaccine dose 2",
          product == "Pfizer" & dose == 1 ~ "Pfizer vaccine dose 1"
        )
      ) %>%
      mutate(
        outcome = str_to_sentence(outcome),
        outcome = factor(
          outcome,
          levels = c(
            "Death",
            "Hospitalisation",
            "Symptoms",
            "Acquisition",
            "Transmission"
          )
        )
      )

    plot <- plot +
      geom_errorbarh(
        aes(
          xmin = days_earliest,
          xmax = days_latest,
          y = ve
        ),
        height = 0,
        alpha = 0.25,
        size = 1,
        data = ve_data_plotting
      ) +
      geom_errorbar(
        aes(
          ymin = ve_lower,
          ymax = ve_upper,
          x = days
        ),
        width = 7,
        alpha = 0.5,
        data = ve_data_plotting
      ) +
      geom_point(
        aes(
          y = ve
        ),
        shape = 21,
        data = ve_data_plotting
      )
  }

  plot

}
