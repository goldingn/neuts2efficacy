# use model estimates to compute vaccine reductions in TP for Aus, and the
# impacts of NPIs

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# get aggregated vaccination scenarios for Australia (very slow first time, then cached)
vaccine_scenarios <- get_vaccine_scenarios()

# calculate the cohorts as at the current date
target_date <- as.Date("2022-03-31")
vaccine_cohorts_now <- get_vaccine_cohorts_at_date(
  vaccine_scenarios,
  target_date
)

coverage_now <- get_coverage(vaccine_cohorts_now)

# write_csv(
#   coverage_now,
#   "outputs/current_coverage.csv"
# )

ves_now <- get_vaccine_efficacies(vaccine_cohorts_now)

# write_csv(
#   ves_now,
#   "outputs/current_ves.csv"
# )

vaccine_transmission_effects_now <- get_vaccine_transmission_effects(
  ves = ves_now,
  coverage = coverage_now
)

# plot the vaccination effects by scenario
transmission_reduction_plot <- vaccine_transmission_effects_now %>%
  filter(
    scenario %in% c(1, 2, 13, 14)
  ) %>%
  mutate(
    scenario = factor(scenario),
    omicron_scenario = factor(
      omicron_scenario,
      levels = c(
        "optimistic",
        "intermediate",
        "pessimistic"
      ),
      labels = c(
        "Optimistic",
        "Intermediate",
        "Pessimistic"
      )
    )
  ) %>%
  ggplot(
    aes(
      x = scenario,
      y = vaccination_effect,
      colour = variant
    )
  ) +
  facet_grid(
    ~omicron_scenario
  ) +
  geom_point(
    size = 3
  ) +
  scale_y_continuous(
    labels = scales::percent_format(
      accuracy = 1
    ),
    limits = c(0, 1)
  ) +
  ylab("Reduction in transmission") +
  xlab("Scenario number") +
  ggtitle(
    "Reduction in transmission potential due to vaccination",
    paste("on", format(target_date, "%b %d %Y"))
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank()
  )

ggsave(
  "figures/transmission_reduction_scenarios.png",
  transmission_reduction_plot,
  bg = "white",
  width = 8,
  height = 4
)

omicron_R0_ratio <- get_omicron_params_wide() %>%
  select(
    omicron_scenario,
    R0_ratio
  )


vaccine_transmission_effects_now %>%
  left_join(
    omicron_R0_ratio,
    by = "omicron_scenario"
  ) %>%
  full_join(
    tibble(
      ascertainment = c(0, 0.5, 1)
    ),
    by = character()
  ) %>%
  mutate(
    R0_ratio = if_else(
      variant == "Delta",
      1,
      R0_ratio
    ),
    # use previously estimated Delta R0 of 8
    R0 = 8 * R0_ratio,
    # effect of baseline distancing behaviours (reduced to R0 6.32 in previous work)
    baseline_distancing_effect = 1 - 6.32 / 8,
    # partial TTIQ effect for detected cases (estimated at 43%), modified by ascertainment
    ttiq_effect = 0.43 * ascertainment,
    # calculate Rts for plotting
    Rt_baseline = R0 * (1 - baseline_distancing_effect),
    Rt_baseline_ttiq = Rt_baseline * (1 - ttiq_effect),
    Rt_baseline_ttiq_vaccination = Rt_baseline_ttiq * (1 - vaccination_effect),
  ) %>%
  # summarise table as output
  filter(
    scenario %in% c(1, 2, 13, 14),
    variant == "Omicron"
  ) %>%
  select(
    scenario,
    omicron_scenario,
    R0,
    vaccination_effect,
    ttiq_effect,
    TP = Rt_baseline_ttiq_vaccination
  ) %>%
  arrange(
    omicron_scenario,
    scenario
  ) %>%
  print(n = Inf)

# make upside down bar plots
