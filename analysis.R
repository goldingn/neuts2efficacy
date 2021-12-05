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

summarise_fit(draws)

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
omicron_kde <- omicron_params %>%
  relocate(
    immune_evasion,
    .after = R0
  ) %>%
  fit_kde()

omicron_params %>%
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
  # scale_y_continuous(
  #   breaks = seq(0, 5, by = 1),
  #   labels = 10 ^ seq(0, 25, by = 5)
  # ) +
  # scale_y_continuous(
  #   labels = scales::percent
  # ) +
  coord_cartesian(
    xlim = c(0, 9),
    # ylim = c(0, 1)
  ) +
  geom_vline(
    xintercept = 6,
    linetype = 2
  ) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  ggtitle("Immune escape and R0 of Omicron vs. Delta") +
  xlab("R0") +
  # ylab("fold reduction in neutralising antibody titre") +
  theme_minimal()

ggsave("figures/contours.png",
       width = 6,
       height = 4,
       bg = "white")

# to do:

# (why) does the population immune fraction dispapear in my model? does it make
# sense to multiply the overall VE for trransmission against this? this
# ignores mixing and the probability that a transmission event occurs?

# if there is no population immunity, the immune escape is irrelevant, so it
# must be explained by R0

# if there is full population immunity, the immune escape is a strong effect, so
# it is explained by both

R = R0 * behavioural_reduction * p_immune * ve_immmunity

# is it linear like this?
(p_immune * reduction_in_transmission_if_immune) + (1 - p_immune) * 1

immunity_effect <- (p_immune * reduction_in_transmission_if_immune + (1 - p_immune))



# express escape as % evasion of immunity, not as fold titre
# VEs for transmission? But bimodal...

# try tighter prior to constrain degree of immune escape

# try running with tighter Reff ratio

# do TP reductions

# add data for variants as neut titre fold of WT, to predict to Delta with C50s

# evaluate prediction of VEs to new variants (alpha, beta) based on variant neut
# titres
