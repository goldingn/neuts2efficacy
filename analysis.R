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

# predict VEs over time for different types of immunity
ve_predictions <- predict_ves(neut_model, draws)
write_csv(
  ve_predictions,
  "outputs/ve_waning_predictions.csv"
)

# plot model predictions waning, with data overlaid
waning_plot <- plot_waning(ve_predictions, neut_model$ve_data_modelling)
ggsave("figures/ve_waning.png",
       plot = waning_plot,
       width = 9,
       height = 6,
       bg = "white")


# plot posterior density of R0 and immune escape


# get posteriors over these and construct a KDE
omicron_params <- sim_omicron_params(neut_model, draws)

# plot in LHSTM style
omicron_kde_lshtm <- omicron_params %>%
  relocate(
    immune_evasion,
    R0_ratio,
    .before = everything()
  ) %>%
  fit_kde()

omicron_params_lshtm_style <- ggplot() +
  geom_polygon(
    aes(x, y),
    data = get_kde_contour(omicron_kde_lshtm, 0.95),
    fill = lighten("seagreen", 0.6)
  ) +
  geom_polygon(
    aes(x, y),
    data = get_kde_contour(omicron_kde_lshtm, 0.5),
    fill = lighten("seagreen", 0.1)
  ) +
  scale_x_continuous(
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

ggsave("figures/omicron_params_lshtm_style.png",
       plot = omicron_params_lshtm_style,
       width = 7,
       height = 7,
       bg = "white")


omicron_kde_bedford <- omicron_params %>%
  relocate(
    R0,
    immune_evasion,
    .before = everything()
  ) %>%
  fit_kde(
    smooth = c(1.8, 0.3)
  )

omicron_params_bedford_style <- ggplot() +
  geom_polygon(
    aes(x, y),
    data = get_kde_contour(omicron_kde_bedford, 0.95),
    fill = lighten("seagreen", 0.6)
  ) +
  geom_polygon(
    aes(x, y),
    data = get_kde_contour(omicron_kde_bedford, 0.5),
    fill = lighten("seagreen", 0.1)
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  geom_vline(
    xintercept = 6,
    linetype = 2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  coord_cartesian(
    xlim = c(1, 12)
  ) +
  ggtitle("Inferred immune escape and relative transmissibility of Omicron vs. Delta") +
  xlab("Omicron transmissibility (R0)") +
  ylab("Omicron immune evasion (reduction in overall VE against transmission, compared to Delta)") +
  theme_minimal()

ggsave("figures/omicron_params_bedford_style.png",
       plot = omicron_params_bedford_style,
       width = 7,
       height = 7,
       bg = "white")


# to do:

# do TP reductions

# add data for variants as neut titre fold of WT, to predict to Delta with C50s

# evaluate prediction of VEs to new variants (alpha, beta) based on variant neut
# titres
