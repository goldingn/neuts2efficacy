# calibrate Deb's model against VE waning estimates
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# standard deviation of population distribution in log10 neutralising antibody
# titres
sd_log10_neut_titres <- 0.4647092
# estimated slope of the sigmoidal relationship between neuts and VEs
log_k <- 1.130661
# estimated log10 neut titre C50 for protection against symptomatic disease
c50_symptoms <- -0.6966127

# list VE against symptomatic disease for different products and doses
# from national plan parameters doc
ve_estimates <- tibble::tribble(
  ~outcome, ~product, ~dose, ~ve,
  "acquisition", "AZ", 1, 0.46,
  "acquisition", "AZ", 2, 0.67,
  "acquisition", "Pfizer", 1, 0.57,
  "acquisition", "Pfizer", 2, 0.8,
  "transmission", "AZ", 1, 0.05,  # 0.02
  "transmission", "AZ", 2, 0.24,  # 0.36
  "transmission", "Pfizer", 1, 0.17, # 0.13
  "transmission", "Pfizer", 2, 0.50, # 0.65
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
)

# compute neutralising antibody titres induced by AZ and Pfizer shortly after
# 1st and 2nd doses:
optimal_log10_neuts <- ve_estimates %>%
  filter(
    outcome == "symptoms"
  ) %>%
  rowwise() %>%
  mutate(
    optimal_log10_neut = find_log10_neut_for_product(
      ve_symptomatic = ve,
      c50_symptoms = c50_symptoms
    )
  ) %>%
  select(
    -outcome,
    -ve
  ) %>%
  ungroup() %>%
  # add a booster assumption - 5-fold the neuts of Pfizer dose 2
  bind_rows(
    filter(.,
           product == "Pfizer",
           dose == 2) %>%
      mutate(
        dose = 3,
        optimal_log10_neut = log10(5 * 10 ^ optimal_log10_neut)
      )
  )

# calibrate c50s for the other parameters from their initial VEs and the modeled neut titres
c50_estimates = ve_estimates %>%
  left_join(
    optimal_log10_neuts,
    by = c("product", "dose")
  ) %>%
  group_by(outcome) %>%
  summarise(
    c50 = find_c50(ve, optimal_log10_neut),
    .groups = "drop"
  )

# compare observed and predicted VEs
evaluate_ves <- ve_estimates %>%
  left_join(
    optimal_log10_neuts,
    by = c("product", "dose")
  ) %>%
  left_join(
    c50_estimates,
    by = "outcome"
  ) %>%
  rowwise() %>%
  mutate(
    predicted_ve = ve_from_mean_log10_neut(
      mean_log10_neut = optimal_log10_neut,
      sd_log10_neut = sd_log10_neut_titres,
      log_k = log_k,
      c50 = c50
    ),
    .after = ve
  ) %>%
  ungroup() %>%
  select(
    -optimal_log10_neut,
    -c50
  ) %>%
  mutate(
    difference = predicted_ve - ve,
    percent_difference = 100 * (1 - predicted_ve / ve)
  )

# calibrate the neut decay estimate for vaccination-derived immunity

# data for 16+ dose 2 waning from Andrews
ve_waning <- tibble::tribble(
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
    ))

# there are some small differences in these initial VEs - likely due to
# differences in definition of outcomes, so include a correction factor when
# fitting
ve_estimates %>%
  filter(
    outcome %in% c("symptoms", "hospitalisation"),
    dose == 2
  ) %>%
  select(
    -dose
  )
ve_waning %>%
  filter(
    weeks == "1"
  ) %>%
  select(
    -weeks,
    -days
  )

# estimate decay rate by least squares
decay_rate_estimate <- ve_waning %>%
  mutate(dose = 2) %>%
  left_join(
    optimal_log10_neuts,
    by = c("product", "dose")
  ) %>%
  left_join(
    c50_estimates,
    by = "outcome"
  ) %>%
  summarise(
    neut_decay_rate = find_decay_rate(.)
  ) %>%
  mutate(
    neut_half_life = log(2) / neut_decay_rate
  )

decay_rate_estimate

# very similar daily rate of decay of neutralising antibody titre from natural
# infection by Khoury et al. (108 day half life)
log(2) / 108


# create a table of all the outcomes, products, doses, and days post dose and the corresponding VEs
# later: add in lag to initial immunity
ve_predictions <- expand_grid(
  outcome = unique(ve_estimates$outcome),
  product = unique(ve_estimates$product),
  dose = 3:1,
  days = 0:365
) %>%
  # remove the dose 3 AZ entry, since it doesn't exist
  filter(
    !(product == "AZ" & dose == 3)
  ) %>%
  left_join(
    optimal_log10_neuts,
    by = c("product", "dose")
  ) %>%
  left_join(
    c50_estimates,
    by = "outcome"
  ) %>%
  mutate(
    current_log10_neut = log10_neut_over_time(
      time = days,
      maximum_log10_neut = optimal_log10_neut,
      decay = decay_rate_estimate$neut_decay_rate
    )
  ) %>%
  rowwise() %>%
  mutate(
    ve = ve_from_mean_log10_neut(
      mean_log10_neut = current_log10_neut,
      sd_log10_neut = sd_log10_neut_titres,
      log_k = log_k,
      c50 = c50
    )
  ) %>%
  ungroup() %>%
  select(
    -optimal_log10_neut,
    -c50,
    -current_log10_neut
  )

write_csv(
  ve_predictions,
  "outputs/ve_waning_predictions.csv"
)

ve_predictions %>%
  mutate(
    dose = factor(
      dose,
      levels = 3:1
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
    )
  ) %>%
  ggplot(
    aes(
      x = days,
      y = ve,
      colour = product,
      linetype = dose
    )
  ) +
  facet_wrap(~outcome) +
  ylab("Vaccine efficacy") +
  xlab("Days since vaccination") +
  geom_line() +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal()

ggsave("figures/ve_waning.png",
       width = 9,
       height = 6,
       bg = "white")

