# calibrate Deb's model against VE waning estimates

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# estimates from Khoury et al
khoury_natmed_estimates <- readRDS("data/SummaryTable_Efficacy_NeutRatio_SD_SEM.RDS")

# get distributions over ratios of neutralising antibody titres to convalescents
# (absolute neut titres are not comparable between studies, so relative to
# convalescent/other vaccine is the only comparable metric)
neut_ratios <- get_khoury_neut_ratios(khoury_natmed_estimates)

# VEs against different outcomes disease for different products and doses
# from national plan parameters doc
ve_estimates <- tibble::tribble(
  ~outcome, ~product, ~dose, ~ve,
  "acquisition", "AZ", 1, 0.46,
  "acquisition", "AZ", 2, 0.67,
  "acquisition", "Pfizer", 1, 0.57,
  "acquisition", "Pfizer", 2, 0.8,
  "transmission", "AZ", 1, 0.05,
  "transmission", "AZ", 2, 0.24,
  "transmission", "Pfizer", 1, 0.17,
  "transmission", "Pfizer", 2, 0.50,
  "symptoms", "AZ", 1, 0.4,
  "symptoms", "AZ", 2, 0.71,
  "symptoms", "Pfizer", 1, 0.58,
  "symptoms", "Pfizer", 2, 0.84,
  "hospitalisation", "AZ", 1, 0.81,
  "hospitalisation", "AZ", 2, 0.93,
  "hospitalisation", "Pfizer", 1, 0.92,
  "hospitalisation", "Pfizer", 2, 0.97,
  "death", "AZ", 1, 0.88,
  "death", "AZ", 2, 0.93,
  "death", "Pfizer", 1, 0.89,
  "death", "Pfizer", 2, 0.95
) %>%
  mutate(
    days = 0,
    .before = ve
  )

# data for 16+ dose 2 waning from Andrews
ve_waning_andrews <- tibble::tribble(
  ~outcome, ~product, ~weeks, ~ve,
  "symptoms", "AZ", "1", 0.627,
  "symptoms", "AZ", "2-9", 0.667,
  "symptoms", "AZ", "10-14", 0.593,
  "symptoms", "AZ", "15-19", 0.526,
  "symptoms", "AZ", "20+", 0.473,
  "symptoms", "Pfizer", "1", 0.924,
  "symptoms", "Pfizer", "2-9", 0.898,
  "symptoms", "Pfizer", "10-14", 0.803,
  "symptoms", "Pfizer", "15-19", 0.734,
  "symptoms", "Pfizer", "20+", 0.697,
  "hospitalisation", "AZ", "1", 0.939,
  "hospitalisation", "AZ", "2-9", 0.952,
  "hospitalisation", "AZ", "10-14", 0.914,
  "hospitalisation", "AZ", "15-19", 0.868,
  "hospitalisation", "AZ", "20+", 0.770,
  "hospitalisation", "Pfizer", "1", 0.997,
  "hospitalisation", "Pfizer", "2-9", 0.984,
  "hospitalisation", "Pfizer", "10-14", 0.965,
  "hospitalisation", "Pfizer", "15-19", 0.944,
  "hospitalisation", "Pfizer", "20+", 0.927,
  "death", "AZ", "2-9", 0.941,
  "death", "AZ", "10-14", 0.924,
  "death", "AZ", "15-19", 0.891,
  "death", "AZ", "20+", 0.787,
  "death", "Pfizer", "2-9", 0.982,
  "death", "Pfizer", "10-14", 0.952,
  "death", "Pfizer", "15-19", 0.939,
  "death", "Pfizer", "20+", 0.904,
) %>%
  mutate(
    # add dose and days info
    dose = 2,
    .after = product
  ) %>%
  mutate(
    # convert 'weeks' to days post optimum (14 days after second dose) text says
    # week 1 is 7-13 days, so presumably week 0 is 0-6 days after the 14-day
    # lag. Compute midpoint in days of other periods accordingly. For 20+ weeks,
    # assume the midpoint of 20-26 weeks. Data includes events up to 3rd
    # September 2021, and most vaccinations starting Jan 4 2021 (supp 1), so
    # subtracting 6 weeks for time to immunity, longest post-vaccination period
    # could be is ~26 months.
    days = case_when(
      weeks == "1" ~ mean(7:13),
      weeks == "2-9" ~ mean((7*2):(7*10 - 1)),
      weeks == "10-14" ~ mean((7*10):(7*15 - 1)),
      weeks == "15-19" ~ mean((7*15):(7*20 - 1)),
      weeks == "20+" ~ mean((7*20):(7*27 - 1)),
    ),
    .after = dose
  ) %>%
  select(
    -weeks
  )


# build a greta model to infer mean neuts, c50s, and the slope of the
# relationship between VEs and neut titres

# record unique values of levels; parameter orders correspond to these
outcomes <- unique(ve_estimates$outcome)
doses <- unique(ve_estimates$dose)
products <- unique(ve_estimates$product)

# set up ve data for modelling
ve_data_modelling <- ve_estimates %>%
  bind_rows(
    ve_waning_andrews
  ) %>%
  mutate(
    # get indices to outcomes, products, doses, to pull together predictions
    outcomes_idx = match(outcome, outcomes),
    doses_idx = match(dose, doses),
    product_idx = match(product, products),
    # logit transform ves to better define likelihood
    ve_logit = qlogis(ve)
  )

# get index to neut ratios, to ensure they are in the correct order
neut_ratios_vaccine <- neut_ratios %>%
  filter(product != "infection")
neut_ratio_vaccine_idx <- match(neut_ratios_vaccine$product, products)

# define greta model

# halflife of the neutralising antibody titre. Khoury found 108 days from
# convalescent this was initialising a bit strangely when defined as N(108, 10),
# so define parameter as standard normmal and transform
neut_halflife_raw <- normal(0, 1)
neut_halflife <- 108 + neut_halflife_raw * 10

# convert to exponential rate parameter
neut_decay <- log(2) / neut_halflife

# log slope of the logit mapping from neuts to VEs
log_k <- normal(0, 1)

# observation error on logit VEs
ve_logit_obs_sd <- normal(0, 1, truncation = c(0, Inf))

# c50s for different outcomes
c50s <- normal(0, 1, dim = 5)

# mean log10 peak neuts for second doses of AZ and Pfizer use the Khoury mean
# (of the log10 ratio to convalescent) and standard error (of the ratio to
# convalescent)
dose_2_mean_log10_neuts <- normal(
  mean = neut_ratios_vaccine$mean_log10_ratio_neut[neut_ratio_vaccine_idx],
  sd = neut_ratios_vaccine$sem_log10_neut[neut_ratio_vaccine_idx]
)

# difference in peak neut titre between first and second doses of each
# (constrain dose 2 to not be lower than dose 1)
dose_2_vs_1_mean_log10_neut_increase <- normal(0, 1, dim = 2, truncation = c(0, Inf))

# mean log10 peak neuts for first doses of AZ and Pfizer respectively
dose_1_mean_log10_neuts <- dose_2_mean_log10_neuts - dose_2_vs_1_mean_log10_neut_increase

# matrix of dose-by-product neuts for lookup; rows are products, columns are
# doses
vaccine_peak_mean_log10_neuts_mat <- cbind(
  dose_1_mean_log10_neuts, dose_2_mean_log10_neuts
)

# pull out vectors of parameters
vaccine_idx <- cbind(
  ve_data_modelling$product_idx,
  ve_data_modelling$doses_idx
)

c50_vec <- c50s[ve_data_modelling$outcomes_idx]
peak_mean_log10_neut_vec <- vaccine_peak_mean_log10_neuts_mat[vaccine_idx]

# mean_log10_neut_vec <- peak_mean_log10_neut_vec
# peak_mean_log10_neut_vec + log10(exp(-neut_decay * ve_data_modelling$days))
mean_log10_neut_vec <- log10_neut_over_time(
  time = ve_data_modelling$days,
  maximum_log10_neut = peak_mean_log10_neut_vec,
  decay = neut_decay
)

# population standard deviation of log10 neut titres
sd_log10_neut_titres <- neut_ratios$sd_log10_ratio_neut[1]

# expected VEs for all combinations
ve_expected <- ve_from_mean_log10_neut(
  mean_log10_neut_vec = mean_log10_neut_vec,
  c50_vec = c50_vec,
  sd_log10_neut = sd_log10_neut_titres,
  log_k = log_k,
  method = "gaussian"
)

# convert to logit scale
ve_expected_logit <- log(ve_expected / (1 - ve_expected))

# likelihood for VEs on logit scale, to better capture variation on observation standard error
distribution(ve_data_modelling$ve_logit) <- normal(ve_expected_logit, ve_logit_obs_sd)

# define and fit model
m <- model(
  log_k,
  neut_halflife,
  ve_logit_obs_sd,
  c50s,
  dose_1_mean_log10_neuts,
  dose_2_mean_log10_neuts
)
draws <- mcmc(m, chains = 10)

# check fit
coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)
coda::effectiveSize(draws)
bayesplot::mcmc_trace(draws)
summary(draws)

# to do:

# add booster and convalescent mean log10 neut variables post-hoc. 5-fold neuts
# of pfizer dose 2 for booster (add log10(5) to log10 neut) and 0 log10 neuts
# for convalescent

# arrange products and doses into a vector for easier lookup

# organise the code to make ve predictions based on outcome, product, dose, time
# since dose etc. (need to add in decay parameter)

# predict to full range of values to plot

# show each as a ribbon (red for AZ, purple for Pfizer, grey for infection) of
#   CIs
# show the number of doses by the darkness of the colours
# add dots for the assumed initial VEs from the parameters doc, and show the
#   Andrews waning data, with bars

# then do TP reductions


ve_expected_sims <- calculate(ve_expected, values = draws, nsim = 1000)[[1]][, , 1]
ve_predictions <- colMeans(ve_expected_sims)
ve_predictions_quants <- apply(ve_expected_sims, 2, quantile, c(0.025, 0.975))


# compare observed and predicted VEs
ve_data_modelling %>%
  mutate(
    predicted_ve = ve_predictions,
    predicted_ve_lower = ve_predictions_quants[1, ],
    predicted_ve_upper = ve_predictions_quants[2, ],
    .after = ve
  ) %>%
  mutate(
    difference = predicted_ve - ve,
    percent_difference = 100 * (1 - predicted_ve / ve)
  ) %>%
  print(
    n = Inf
  ) %>%
  ggplot(
    aes(
      x = ve,
      y = predicted_ve,
      colour = days
    )
  ) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = predicted_ve_lower,
      ymax = predicted_ve_upper,
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1
  ) +
  theme_minimal()

# posterior predictive check on ECDF of observed VEs
bayesplot::ppc_ecdf_overlay(ve_data_modelling$ve, ve_expected_sims)






#
# # create a table of all the outcomes, products, doses, and days post dose and the corresponding VEs
# # later: add in lag to initial immunity
# ve_predictions <- expand_grid(
#   outcome = unique(ve_estimates$outcome),
#   product = unique(ve_estimates$product),
#   dose = 3:1,
#   days = 0:365
# ) %>%
#   # remove the dose 3 AZ entry, since it doesn't exist
#   filter(
#     !(product == "AZ" & dose == 3)
#   ) %>%
#   left_join(
#     optimal_log10_neuts,
#     by = c("product", "dose")
#   ) %>%
#   left_join(
#     c50_estimates,
#     by = "outcome"
#   ) %>%
#   mutate(
#     current_log10_neut = log10_neut_over_time(
#       time = days,
#       maximum_log10_neut = optimal_log10_neut,
#       decay = decay_rate_estimate$neut_decay_rate
#     )
#   ) %>%
#   rowwise() %>%
#   mutate(
#     ve = ve_from_mean_log10_neut(
#       mean_log10_neut = current_log10_neut,
#       sd_log10_neut = sd_log10_neut_titres,
#       log_k = log_k,
#       c50 = c50
#     )
#   ) %>%
#   ungroup() %>%
#   select(
#     -optimal_log10_neut,
#     -c50,
#     -current_log10_neut
#   )
#
# write_csv(
#   ve_predictions,
#   "outputs/ve_waning_predictions.csv"
# )
#
# ve_predictions %>%
#   mutate(
#     dose = factor(
#       dose,
#       levels = 3:1
#     ),
#     outcome = factor(
#       outcome,
#       levels = c(
#         "death",
#         "hospitalisation",
#         "symptoms",
#         "acquisition",
#         "transmission"
#       )
#     )
#   ) %>%
#   ggplot(
#     aes(
#       x = days,
#       y = ve,
#       colour = product,
#       linetype = dose
#     )
#   ) +
#   facet_wrap(~outcome) +
#   ylab("Vaccine efficacy") +
#   xlab("Days since vaccination") +
#   geom_line() +
#   scale_y_continuous(labels = scales::percent) +
#   theme_minimal()
#
# ggsave("figures/ve_waning.png",
#        width = 9,
#        height = 6,
#        bg = "white")

