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
neut_model <- build_neut_model(ve_estimates, neut_ratios_vaccine)

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

# to do:

# do TP reductions

# add data for variants as neut titre fold of WT, to predict to Delta with C50s

# evaluate prediction of VEs to new variants (alpha, beta) based on variant neut
# titres
