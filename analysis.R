# calibrate Deb's model against VE waning estimates

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# get distributions over ratios of neutralising antibody titres to convalescents
# from Khoury et al. (absolute neut titres are not comparable between studies,
# so these are expressed as fold of convalescent)
neut_ratios_vaccine <- get_neut_ratios_vaccine()

# VEs against different outcomes for different products and doses, on different
# days since peak vaccination impact
ve_estimates <- get_ve_estimates()

# build a greta model of relationship between neuts and efficacies
neut_model_initial <- build_neut_model(ve_estimates, neut_ratios_vaccine)

# extend the model to add omicron parameters
neut_model <- add_omicron_model(neut_model_initial)

# fit model by MCMC and check convergence
draws <- fit(neut_model)

summarise_fit(draws) %>%
  print(n = Inf)

# print parameters for simulation model
format_parameters(neut_model, draws) %>%
  as.data.frame()

trace_plot <- bayesplot::mcmc_trace(draws)
ggsave("figures/mcmc_trace.png",
       plot = trace_plot,
       width = 10,
       height = 7,
       bg = "white")

# get posterior samples of fitted VEs to check fit
ve_fitted_sims <- sim_fitted_ves(neut_model, draws)

# plot visualisation of fit to training data
fit_plot <- plot_ve_fit(neut_model$ve_data_modelling, ve_fitted_sims)
ggsave("figures/ve_fit.png",
       plot = fit_plot,
       width = 7,
       height = 7,
       bg = "white")

# plot posterior predictive check on ECDF of observed VEs
ppc_plot <- bayesplot::ppc_ecdf_overlay(neut_model$ve_data_modelling$ve, ve_fitted_sims)
ggsave("figures/ecdf_ppc.png",
       plot = ppc_plot,
       width = 6,
       height = 5,
       bg = "white")

# predict VEs against Delta over time for different types of immunity
ve_predictions_delta <- predict_ves(neut_model, draws)
write_csv(
  ve_predictions_delta,
  "outputs/ve_waning_predictions_delta.csv"
)

# plot model predictions waning, with data overlaid
waning_plot_delta <- plot_waning(ve_predictions_delta) +
  ggtitle("Predicted waning in vaccine efficacy",
          "against the Delta variant")

ggsave("figures/ve_waning_delta.png",
       plot = waning_plot_delta,
       width = 9,
       height = 6,
       bg = "white")

delta_ve_data <- neut_model$ve_data_modelling %>%
  filter(variant == "delta")

waning_plot_delta_data <- plot_waning(ve_predictions_delta, delta_ve_data) +
  ggtitle("Predicted waning in vaccine efficacy",
          "against the Delta variant")

ggsave("figures/ve_waning_delta_with_data.png",
       plot = waning_plot_delta_data,
       width = 9,
       height = 6,
       bg = "white")

# make ve predictions for Omicron
ve_predictions_omicron <- predict_ves(neut_model, draws, omicron = TRUE)
write_csv(
  ve_predictions_omicron,
  "outputs/ve_waning_predictions_omicron.csv"
)

waning_plot_omicron <- plot_waning(
  ve_predictions_omicron,
  immunity_levels = c(
    "mRNA booster",
    "Pfizer vaccine dose 2",
    "AZ vaccine dose 2")
) +
  ggtitle("Predicted waning in vaccine efficacy",
          "against the Omicron variant")

ggsave("figures/ve_waning_omicron.png",
       plot = waning_plot_omicron,
       width = 9,
       height = 6,
       bg = "white")

# get preliminary VEs for Omicron from Andrews
omicron_ve_data <- neut_model$ve_data_modelling %>%
  filter(
    variant == "omicron",
    # drop dose 1 as we won't plot the ribbon
    dose > 1
  )

waning_plot_omicron_data <- plot_waning(
  ve_predictions_omicron,
  omicron_ve_data,
  immunity_levels = c(
    "mRNA booster",
    "Pfizer vaccine dose 2",
    "AZ vaccine dose 2")
) +
  ggtitle("Predicted waning in vaccine efficacy",
          "against the Omicron variant")

ggsave("figures/ve_waning_omicron_with_data.png",
       plot = waning_plot_omicron_data,
       width = 9,
       height = 6,
       bg = "white")

# plot posterior density of R0 and immune escape
omicron_params <- sim_omicron_params(neut_model, draws)

omicron_params_plot <- omicron_params %>%
  plot_omicron_params()

ggsave("figures/omicron_params.png",
       plot = omicron_params_plot,
       width = 7,
       height = 7,
       bg = "white")

omicron_params_plot_lines <- omicron_params %>%
  plot_omicron_params(
    neut_fold_lines = c(3, 5, 10, 15)
  )

ggsave("figures/omicron_params_with_lines.png",
       plot = omicron_params_plot_lines,
       width = 7,
       height = 7,
       bg = "white")


omicron_neut_fold_plot <- plot_omicron_neut_fold(omicron_params)

ggsave("figures/omicron_neut_fold.png",
       plot = omicron_neut_fold_plot,
       width = 7,
       height = 7,
       bg = "white")



# given point estimates on the avocado, compute VEs for vaccines and waning for
# Omicron in Australia
omicron_scenarios <- tibble::tribble(
  ~scenario, ~R0_ratio_target, ~immune_evasion_target,
  "high_R0", 1.1, 0.25,
  "intermediate", 0.65, 0.5,
  "high_evasion", 0.5, 0.8
)

ve_predictions_omicron_scenarios <- predict_ve_scenarios(
  scenarios = omicron_scenarios,
  neut_model = neut_model,
  draws = draws,
  omicron = TRUE
)

write_csv(
  ve_predictions_omicron_scenarios,
  "outputs/ve_waning_predictions_omicron_scenarios.csv"
)


ve_predictions_omicron_scenarios %>%
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "high_R0",
        "intermediate",
        "high_evasion"
      )
    ),
    outcome = factor(
      outcome,
      levels = c(
        "death",
        "hospitalisation",
        "symptoms",
        "acquisition",
        "transmission"
      )
    ),
    immunity = factor(
      immunity,
      levels = c(
        "booster",
        "pfizer_dose_2",
        "infection",
        "az_dose_2",
        "pfizer_dose_1",
        "az_dose_1"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = days,
      y = ve,
      colour = immunity
    )
  ) +
  geom_line(
    size = 1,
    alpha = 0.8
  ) +
  facet_grid(
    outcome ~ scenario
  ) +
  ylab("Efficacy") +
  xlab("Days since peak immunity") +
  theme_minimal()


