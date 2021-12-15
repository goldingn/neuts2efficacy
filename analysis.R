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


# wrap up these waned immune evasions in a plot

# do the same but for a cohort with the same level of immunity as AZ dose 1 (or
# AZ dose 2 waned by 6 months)
omicron_params_plot_AZ_2_180d <- sim_omicron_params(
  neut_model,
  draws,
  baseline_immunity = "AZ_dose_2",
  days_waning = 180
) %>%
  plot_omicron_params()

ggsave("figures/omicron_params_az2_180d.png",
       plot = omicron_params_plot_AZ_2_180d,
       width = 7,
       height = 7,
       bg = "white")

# for Pfizer dose 1
omicron_params_plot_Pfizer_2_180d <- sim_omicron_params(
  neut_model,
  draws,
  baseline_immunity = "Pfizer_dose_2",
  days_waning = 180
) %>%
  plot_omicron_params()

ggsave("figures/omicron_params_pf_2_180d.png",
       plot = omicron_params_plot_Pfizer_2_180d,
       width = 7,
       height = 7,
       bg = "white")

# and for booster
omicron_params_plot_booster_0d <- sim_omicron_params(
  neut_model,
  draws,
  baseline_immunity = "mRNA_booster",
  days_waning = 0
) %>%
  plot_omicron_params()

ggsave("figures/omicron_params_booster_0d.png",
       plot = omicron_params_plot_booster_0d,
       width = 7,
       height = 7,
       bg = "white")

omicron_params_plot_booster_180d <- sim_omicron_params(
  neut_model,
  draws,
  baseline_immunity = "mRNA_booster",
  days_waning = 180
) %>%
  plot_omicron_params()

ggsave("figures/omicron_params_booster_180d.png",
       plot = omicron_params_plot_booster_180d,
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
  "pessimistic", 1.1, 0.5,
  "intermediate", 0.85, 0.3,
  "optimistic", 0.6, 0.1
)

omicron_params_plot_points <- omicron_params_plot +
  geom_point(
    aes(
      x = immune_evasion_target,
      y = R0_ratio_target
    ),
    data = omicron_scenarios
  ) +
  geom_text(
    aes(
      x = immune_evasion_target,
      y = R0_ratio_target,
      label = scenario
    ),
    data = omicron_scenarios,
    vjust = 0,
    nudge_y = 0.02
  ) +
  coord_cartesian(
    xlim = c(0, 0.6)
  )

ggsave("figures/omicron_params_with_points.png",
       plot = omicron_params_plot_points,
       width = 7,
       height = 7,
       bg = "white")

# get IBM parameter vallues for these scenarios
scenario_parameters <- format_scenario_parameters(
  neut_model = neut_model,
  draws = draws,
  scenarios = omicron_scenarios
)

write_csv(
  scenario_parameters,
  "outputs/scenario_parameters_omicron.csv"
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
        "pessimistic",
        "intermediate",
        "optimistic"
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
      ),
      labels = c(
        "Death",
        "Hospitalisation",
        "Symptoms",
        "Acquisition",
        "Onward transmission"
      ),
    ),
    immunity = factor(
      immunity,
      levels = c(
        "mRNA_booster",
        "Pfizer_dose_2",
        "infection",
        "AZ_dose_2",
        "Pfizer_dose_1",
        "AZ_dose_1"
      ),
      labels = c(
        "mRNA booster",
        "Pfizer dose 2",
        "Infection",
        "AZ dose 2",
        "Pfizer dose 1",
        "AZ dose 1"
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


# show how immune evasion to VE for overall transmission increases with waning

omicron_ves <- read_csv(
  "outputs/ve_waning_predictions_omicron_scenarios.csv",
  col_types = cols(
    scenario = col_character(),
    outcome = col_character(),
    immunity = col_character(),
    days = col_double(),
    immunity_type = col_character(),
    ve = col_double()
  )
) %>%
  filter(
    scenario == "intermediate"
  ) %>%
  select(
    outcome,
    immunity,
    days,
    ve
  )

delta_ves <- read_csv(
  "outputs/ve_waning_predictions_delta.csv",
  col_types = cols(
    outcome = col_character(),
    immunity = col_character(),
    days = col_double(),
    ve_predict_mean = col_double(),
    ve_predict_sd = col_double(),
    ve_predict_lower_50 = col_double(),
    ve_predict_upper_50 = col_double(),
    ve_predict_lower_90 = col_double(),
    ve_predict_upper_90 = col_double(),
    immunity_type = col_character()
  )
) %>%
  select(
    outcome,
    immunity,
    days,
    ve = ve_predict_mean
  )

evasion_plot <- bind_rows(
  omicron = omicron_ves,
  delta = delta_ves,
  .id = "variant"
) %>%
  filter(
    outcome %in% c("acquisition", "transmission")
  ) %>%
  pivot_wider(
    names_from = outcome,
    values_from = ve
  ) %>%
  mutate(
    overall = 1 - ((1 - acquisition) * (1 - transmission))
  ) %>%
  select(
    -acquisition,
    -transmission
  ) %>%
  pivot_wider(
    names_from = variant,
    values_from = overall
  ) %>%
  mutate(
    evasion = 1 - omicron / delta
  ) %>%
  mutate(
    immunity = factor(
      immunity,
      levels = c(
        "AZ_dose_1",
        "Pfizer_dose_1",
        "AZ_dose_2",
        "infection",
        "Pfizer_dose_2",
        "mRNA_booster"
      ),
      labels = c(
        "AZ dose 1",
        "Pfizer dose 1",
        "AZ dose 2",
        "Infection",
        "Pfizer dose 2",
        "mRNA booster"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = days,
      y = evasion,
      colour = immunity
    )
  ) +
  geom_line(
    size = 2
  ) +
  scale_y_continuous(
    labels = scales::percent,
    limits = c(0, 1)
  ) +
  ylab("Reduction in overall VE against transmission") +
  xlab("Days since peak immunity") +
  ggtitle("Omicron immune evasion as a function of immunity") +
  theme_minimal()

ggsave(
  "figures/evasion_vs_immunity.png",
  evasion_plot,
  width = 7,
  height = 7,
  bg = "white"
)



omicron_params_plot_za_120d <- sim_omicron_params(
  neut_model,
  draws,
  baseline_immunity = "za",
  days_waning = 120
) %>%
  plot_omicron_params(
    neut_fold_lines = c(3, 5, 10, 15)
  )

ggsave("figures/omicron_params_za_120d.png",
       plot = omicron_params_plot_za_120d,
       width = 7,
       height = 7,
       bg = "white")
