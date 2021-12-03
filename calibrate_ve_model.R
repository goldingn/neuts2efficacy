# calibrate Deb's model against VE waning estimates

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# get distributions over ratios of neutralising antibody titres to convalescents
# from Khoury et al. (absolute neut titres are not comparable between studies,
# so relative to convalescent/other vaccine is the only comparable metric)
neut_ratios <- get_khoury_neut_ratios(khoury_natmed_estimates)

# VEs against different outcomes for different products and doses, on different
# days since peak vaccination impact
ve_estimates <- get_ve_estimates()

# build a greta model to infer mean neuts, c50s, and the slope of the
# relationship between VEs and neut titres

# record unique values of levels; parameter orders correspond to these
outcomes <- unique(ve_estimates$outcome)
doses <- unique(ve_estimates$dose)
products <- unique(ve_estimates$product)

# set up ve data for modelling
ve_data_modelling <- ve_estimates %>%
  prep_ve_data_for_modelling() %>%
  mutate(
    # get indices to outcomes, products, doses, to pull together predictions
    outcomes_idx = match(outcome, outcomes),
    doses_idx = match(dose, doses),
    product_idx = match(product, products),
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

# a shared additional observation error on logit VEs (since uncertainty in
# estimates doesn't reflect differences in outcome definitions between studies
# etc.)
ve_logit_obs_sd_shared <- normal(0, 1, truncation = c(0, Inf))

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

# expand out the c50 parameters to match the data
c50_vec <- c50s[ve_data_modelling$outcomes_idx]

# get peak mean log10 neut titres for each vaccine
peak_mean_log10_neut_vec <- vaccine_peak_mean_log10_neuts_mat[vaccine_idx]

# compute expected values on given days post peak
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

# get the observation standard deviation for each observation; combining error
# from both the estimation (known for each estimate) with error between studies
# in outcome definitions and sample biases etc
ve_logit_obs_sd <- sqrt(ve_data_modelling$ve_logit_sd ^ 2 + ve_logit_obs_sd_shared ^ 2)

# likelihood for VEs, with observation error given by estimate SD
distribution(ve_data_modelling$ve_logit) <- normal(ve_expected_logit, ve_logit_obs_sd)

# define and fit model
m <- model(
  log_k,
  neut_halflife,
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

# compare observed and predicted VEs
ve_fitted_sims <- calculate(ve_expected, values = draws, nsim = 1000)[[1]][, , 1]
ve_fitted <- colMeans(ve_fitted_sims)
ve_fitted_quants <- apply(ve_fitted_sims, 2, quantile, c(0.025, 0.975))

ve_data_modelling %>%
  mutate(
    predicted_ve = ve_fitted,
    predicted_ve_lower = ve_fitted_quants[1, ],
    predicted_ve_upper = ve_fitted_quants[2, ],
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
      ymax = predicted_ve_upper
    )
  ) +
  geom_errorbar(
    aes(
      xmin = ve_lower,
      xmax = ve_upper
    )
  ) +
  geom_abline(
    intercept = 0,
    slope = 1
  ) +
  theme_minimal()

# posterior predictive check on ECDF of observed VEs
bayesplot::ppc_ecdf_overlay(ve_data_modelling$ve, ve_fitted_sims)

log10_booster_multipler <- log10(5)

# prepare for prediction
log10_neuts_list <- list(
  az_dose_1 = dose_1_mean_log10_neuts[products == "AZ"],
  az_dose_2 = dose_2_mean_log10_neuts[products == "AZ"],
  pfizer_dose_1 = dose_1_mean_log10_neuts[products == "Pfizer"],
  pfizer_dose_2 = dose_2_mean_log10_neuts[products == "Pfizer"],
  booster = dose_2_mean_log10_neuts[products == "Pfizer"] + log10_booster_multipler,
  infection = 0
)

immunities <- names(log10_neuts_list)
peak_mean_log10_neuts_all <- do.call(c, log10_neuts_list)

# create a table of all the outcomes, immunity types, and days post dose, to compute the VEs
ve_prediction_data <- expand_grid(
  outcome = outcomes,
  immunity = immunities,
  days = 0:365
) %>%
  mutate(
    outcome_idx = match(outcome, outcomes),
    immunity_idx = match(immunity, immunities)
  )

mean_log10_neuts_all <- log10_neut_over_time(
  time = ve_prediction_data$days,
  maximum_log10_neut = peak_mean_log10_neuts_all[ve_prediction_data$immunity_idx],
  decay = neut_decay
)

# predict VEs
ve_predict <- ve_from_mean_log10_neut(
  mean_log10_neut_vec = mean_log10_neuts_all,
  c50_vec = c50s[ve_prediction_data$outcome_idx],
  log_k = log_k,
  sd_log10_neut = sd_log10_neut_titres,
  method = "gaussian"
)

ve_predict_sims <- calculate(ve_predict,
                             values = draws,
                             nsim = 1000)[[1]][, , 1]

ve_predictions <- ve_prediction_data %>%
  mutate(
    ve_predict_mean = colMeans(ve_predict_sims),
    ve_predict_sd = apply(ve_predict_sims, 2, stats::sd),
    ve_predict_lower_50 = apply(ve_predict_sims, 2, quantile, 0.25),
    ve_predict_upper_50 = apply(ve_predict_sims, 2, quantile, 0.75),
    ve_predict_lower_90 = apply(ve_predict_sims, 2, quantile, 0.05),
    ve_predict_upper_90 = apply(ve_predict_sims, 2, quantile, 0.95)
  ) %>%
  mutate(
    immunity_type = case_when(
      immunity == "az_dose_2" ~ "AZ vaccine dose 2",
      immunity == "az_dose_1" ~ "AZ vaccine dose 1",
      immunity == "pfizer_dose_2" ~ "Pfizer vaccine dose 2",
      immunity == "pfizer_dose_1" ~ "Pfizer vaccine dose 1",
      immunity == "booster" ~ "mRNA booster",
      immunity == "infection" ~ "Infection"
    )
  )

write_csv(
  ve_predictions,
  "outputs/ve_waning_predictions.csv"
)

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

ve_predictions %>%
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
      levels = c(
        "mRNA booster",
        "Pfizer vaccine dose 2",
        "Pfizer vaccine dose 1",
        "AZ vaccine dose 2",
        "AZ vaccine dose 1",
        "Infection"
      )
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
    colour = grey(0.4)
  ) +
  geom_errorbarh(
    aes(
      xmin = days_earliest,
      xmax = days_latest,
      y = ve
    ),
    height = 0,
    alpha = 0.5,
    data = ve_data_plotting
  ) +
  geom_errorbar(
    aes(
      ymin = ve_lower,
      ymax = ve_upper,
      x = days
    ),
    width = 0,
    alpha = 0.5,
    data = ve_data_plotting
  ) +
  geom_point(
    aes(
      y = ve
    ),
    shape = 21,
    data = ve_data_plotting
  ) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  scale_fill_manual(
    values = c(
      "mRNA booster" = lighten("darkorchid4", 0.1),
      "Pfizer vaccine dose 2" = lighten("darkorchid1", 0.1),
      "Infection" = grey(0.8),
      "AZ vaccine dose 2" = lighten("firebrick1", 0.1),
      "Pfizer vaccine dose 1" = lighten("darkorchid1", 0.8),
      "AZ vaccine dose 1" = lighten("firebrick1", 0.8)
    )
  ) +
  scale_alpha_manual(
    values = c("two doses" = 0.8, "one dose" = 0.1)
  ) +
  ylab("Efficacy against Delta variant") +
  xlab("Days since vaccination") +
  ggtitle("Population average protection against negative outcomes") +
  theme_minimal()

ggsave("figures/ve_waning.png",
       width = 10,
       height = 6,
       bg = "white")


# to do:

# add ve estimate uncertainty to the plot

# add a hierarchical structure over independent log_ks for outcomes

# do TP reductions

# look at implementing waning
