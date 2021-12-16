# use model estimates to compute vaccine reductions in TP for Aus, and the
# impacts of NPIs

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# get aggregated vaccination scenarios for Australia (very slow first time, then cached)
vaccine_scenarios <- get_vaccine_scenarios()

# calculate the cohorts as at the current date
vaccine_cohorts_now <- get_vaccine_cohorts_at_date(vaccine_scenarios, as.Date("2021-12-08") - 30)

coverage_now <- get_coverage(vaccine_cohorts_now)

write_csv(
  coverage_now,
  "outputs/current_coverage.csv"
)

ves_now <- get_vaccine_efficacies(vaccine_cohorts_now)

write_csv(
  ves_now,
  "outputs/current_ves.csv"
)

vaccine_transmission_effects_now <- get_vaccine_transmission_effects(
  ves = ves_now,
  coverage = coverage_now
)

vaccine_transmission_effects_now %>%
  mutate(
    vaccine_transmission_reduction = 1 - vaccination_effect_multiplier
  ) %>%
  ggplot(
    aes(
      x = scenario,
      y = vaccine_transmission_reduction,
    )
  ) +
  facet_wrap(
    variant ~ omicron_scenario,
    scales = "free_y"
  ) +
  geom_point() +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 0.1
    )
  ) +
  ylab("Reduction in transmission") +
  xlab("Scenario number") +
  theme_minimal()


# pull together NPI assumptions from previous work

# make upside down bar plots
