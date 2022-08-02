library(targets)

# Use tar_script() to create _targets.R and tar_edit()
# to open it again for editing.
# Then, run tar_make() to run the pipeline
# and tar_read(summary) to view the results.

# Set target-specific options such as packages.
tar_option_set(
  packages = c(
    "tidyverse",
    "gaussquad",
    "greta",
    "bayesplot",
    "testthat",
    "colorspace"
  )
)

# load all functions
. <- lapply(list.files("R", full.names = TRUE), source)

# End this file with a list of target objects.
list(

  # distributions over ratios of neut titres from Khoury et al, expressed as fold of convalescent
  tar_target(
    neut_ratios_vaccine,
    get_neut_ratios_vaccine()
  ),

  # VEs against different outcomes for products, doses at different times since vaccination.
  tar_target(
    ve_estimates,
    get_ve_estimates()
  ),

  # greta model of relationship between neuts and efficacies
  tar_target(
    neut_model,
    build_neut_model(
      ve_estimates,
      neut_ratios_vaccine
    )
  ),

  # fit model and check convergence
  tar_target(
    draws,
    fit(neut_model)
  ),

  tar_target(
    trace_plot,
    bayesplot::mcmc_trace(draws)
  ),

  tar_target(
    save_trace_plot,
    ggsave("figures/mcmc_trace.png",
       plot = trace_plot,
       width = 10,
       height = 7,
       bg = "white")
  ),

  # posterior samples of fitted VEs:
  tar_target(
    ve_fitted_sims,
    sim_fitted_ves(
      neut_model,
      draws
    )
  ),

  # plot of fit:
  tar_target(
    fit_plot,
    plot_ve_fit(
      neut_model$ve_data_modelling,
      ve_fitted_sims
    )
  ),

  tar_target(
    save_fit_plot,
    ggsave("figures/ve_fit.png",
       plot = fit_plot,
       width = 7,
       height = 7,
       bg = "white")
  ),

  # posterior predictive check on ECDF of observed VEs
  tar_target(
    ppc_plot,
    bayesplot::ppc_ecdf_overlay(
      neut_model$ve_data_modelling$ve,
      ve_fitted_sims)
  ),

  tar_target(
    save_ppc_plot,
    ggsave("figures/ecdf_ppc.png",
       plot = ppc_plot,
       width = 6,
       height = 5,
       bg = "white")
  ),

  # VEs against Delta over time
  tar_target(
    ve_predictions_delta,
    predict_ves(
      neut_model,
      draws,
      omicron_infection_multiplier = 4.4
    )
  ),
  tar_target(
    save_ve_delta,
    write_csv(
      ve_predictions_delta,
      "outputs/ve_waning_predictions_delta.csv"
    )
  ),

  tar_target(
    waning_plot_delta,
    plot_waning(
      ve_predictions_delta
    ) +
    ggtitle("Predicted waning in vaccine efficacy",
          "against the Delta variant")
  ),
  tar_target(
    save_waning_plot_delta,
    ggsave("figures/ve_waning_delta.png",
       plot = waning_plot_delta,
       width = 9,
       height = 6,
       bg = "white")
  ),

  tar_target(
    waning_plot_delta_data,
    plot_waning(
      ve_predictions_delta,
      neut_model$ve_data_modelling %>%
        filter(variant == "delta")
    ) +
    ggtitle("Predicted waning in vaccine efficacy",
          "against the Delta variant")
  ),

  tar_target(
    save_waning_plot_delta_data,
    ggsave("figures/ve_waning_delta_with_data.png",
       plot = waning_plot_delta_data,
       width = 9,
       height = 6,
       bg = "white")
  ),

  # VE predictions for omicron
  tar_target(
    ve_predictions_omicron,
    predict_ves(
      neut_model,
      draws,
      omicron = TRUE,
      omicron_infection_multiplier = 14.4
    )
  ),

  tar_target(
    save_ve_predictions_omicron,
    write_csv(
      ve_predictions_omicron,
      "outputs/ve_waning_predictions_omicron.csv"
    )
  ),

  tar_target(
    waning_plot_omicron,
    plot_waning(
      ve_predictions_omicron,
      immunity_levels = c(
        "Omicron Infection",
        "Pfizer vaccine dose 2 + Omicron infection",
        "mRNA booster + Omicron infection",
        "mRNA booster",
        "Pfizer vaccine dose 2",
        "AZ vaccine dose 2",
        "mRNA dose 4",
        "mRNA dose 4 + Omicron infection")
      ) +
      ggtitle("Predicted waning in vaccine efficacy",
              "against the Omicron variant")
  ),

  tar_target(
    save_waning_plot_omicron,
    ggsave("figures/ve_waning_omicron.png",
       plot = waning_plot_omicron,
       width = 9,
       height = 6,
       bg = "white")
  ),

  tar_target(
    omicron_ve_data,
    neut_model$ve_data_modelling %>%
      filter(
        variant == "omicron",
        dose > 1
      )
  ),

  tar_target(
    waning_plot_omicron_data,
    plot_waning(
      ve_predictions_omicron,
      omicron_ve_data,
      immunity_levels = c(
        "Omicron Infection",
        "Pfizer vaccine dose 2 + Omicron infection",
        "mRNA booster + Omicron infection",
        "mRNA booster",
        "Pfizer vaccine dose 2",
        "AZ vaccine dose 2")
    ) +
      ggtitle("Predicted waning in vaccine efficacy",
              "against the Omicron variant")
  ),

  tar_target(
    save_waning_plot_omicron_data,
    ggsave("figures/ve_waning_omicron_with_data.png",
       plot = waning_plot_omicron_data,
       width = 9,
       height = 6,
       bg = "white")
  )

)
