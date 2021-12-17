# use model estimates to compute vaccine reductions in TP for Aus, and the
# impacts of NPIs

# load packages and functions
source("packages.R")
. <- list.files("R", full.names = TRUE) %>%
  lapply(source)

# get aggregated vaccination scenarios for Australia (very slow first time, then cached)
vaccine_scenarios <- get_vaccine_scenarios()

# calculate the cohorts as at the current date
# target_date <- as.Date("2022-01-31")
# target_date <- as.Date("2022-03-31")
target_date <- Sys.Date()
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
  paste0("figures/transmission_reduction_scenarios_",
         target_date,
         ".png"),
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

low_phsm_effect <- 1 - 4.746546833 / 6.317224934
medium_phsm_effect <- 1 - 3.842466025 / 6.317224934
high_phsm_effect <- 1 - 2.795920425 / 6.317224934

tps <- vaccine_transmission_effects_now %>%
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
    Rt_baseline_ttiq_vaccination_low = Rt_baseline_ttiq_vaccination * (1 - low_phsm_effect),
    Rt_baseline_ttiq_vaccination_medium = Rt_baseline_ttiq_vaccination * (1 - medium_phsm_effect),
    Rt_baseline_ttiq_vaccination_high = Rt_baseline_ttiq_vaccination * (1 - high_phsm_effect)
  )


# summarise table as output
tps %>%
  filter(
    scenario %in% c(1, 2, 13, 14),
    omicron_scenario == "intermediate",
    variant == "Omicron"
  ) %>%
  select(
    scenario,
    omicron_scenario,
    R0,
    vaccination_effect,
    ttiq_effect,
    TP_baseline = Rt_baseline_ttiq_vaccination,
    TP_low = Rt_baseline_ttiq_vaccination_low,
    TP_medium = Rt_baseline_ttiq_vaccination_medium,
    TP_high = Rt_baseline_ttiq_vaccination_high,
  ) %>%
  arrange(
    omicron_scenario,
    scenario
  ) %>%
  print(n = Inf)

# make a washed-out colour, by a vector of amounts between 0 and 1
washout <- function(colour, amount = 0.7) {

  stopifnot(
    all(amount >= 0) & all(amount <= 1)
  )

  n <- 1000
  indices <- pmax(1, ceiling(amount * n))

  palette <- colorRampPalette(c(colour, "white"))
  palette(n)[indices]

}


add_context_hline <- function(
  p,
  at,
  label,
  colour = grey(0),
  linetype = 2,
  text_size = 2.5
) {
  p +
    geom_hline(
      yintercept = at,
      col = colour,
      linetype = linetype
    ) +
    annotate(
      "text",
      label = label,
      x = -0.5,
      y = at,
      hjust = 0,
      vjust = -0.5,
      col = colour,
      size = text_size
    )
}

control_plot_theme <- function() {
  cowplot::theme_cowplot() +
    # turn off the x axis and add some space for annotation on RHS
    theme(
      plot.margin = unit(
        c(1, 1, 0, 1),
        "line"
      ),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y.right = element_text(
        size = 12,
        margin = margin(l = 10),
      )
    )
}

control_base_plot <- function(
  data,
  ylab = "Transmission potential"
) {
  data %>%
    ggplot(
      aes(
        x = scenario,
        middle = as.numeric(NA),
        ymin = as.numeric(NA),
        ymax = as.numeric(NA),
        width = 0.6
      )
    ) +
    scale_x_discrete() +
    xlab("") +
    scale_y_continuous(
      position = "right",
      breaks = c(0.5, 1, 2, 4, 6, 8),
      trans = 'log'
    ) +
    coord_cartesian(clip = "off") +
    ylab(ylab) +
    control_plot_theme()
}

add_single_box <- function(
  p,
  top,
  bottom,
  text_main = "",
  use_reduction_text = FALSE,
  only_scenarios = c(),
  text_size = 3,
  box_colour = grey(0.9),
  text_colour = grey(0.3),
  border_colour = grey(0.6)
) {

  top <- enquo(top)
  bottom <- enquo(bottom)
  p <- p +
    geom_boxplot(
      aes(
        lower = !!bottom,
        upper = !!top
      ),
      stat = "identity",
      fill = box_colour,
      col = border_colour
    )

  if (use_reduction_text) {

    # restriction labels
    p <- p + geom_text(
      aes(
        label = ifelse(
          scenario %in% only_scenarios,
          sprintf(
            "%s\n%s",
            text_main,
            scenario
          ),
          NA
        ),
        y = !!bottom * (!!top / !!bottom) ^ 0.5  # (midpoint on log scale!)
      ),
      size = text_size,
      col = text_colour
    )
  } else {
    p <- p + geom_text(
      aes(
        label = ifelse(
          scenario %in% only_scenarios,
          text_main,
          NA
        ),
        y = !!bottom * (!!top / !!bottom) ^ 0.5  # (midpoint on log scale!)
      ),
      size = text_size,
      col = text_colour
    )
  }

  p

}

add_stacked_box <- function(
  p,
  bottom,  # variable for bottom of the box
  top,  # variable for top of the box
  reference,  # variable against which to calculate % reduction of 'bottom'
  text_main = "",
  text_size = 3,
  only_scenarios = c(),
  box_colour = grey(0.9),
  border_colour = grey(0.6),
  text_colour = grey(0.3),
  use_reduction_text = FALSE
) {

  bottom <- enquo(bottom)
  top <- enquo(top)
  reference <- enquo(reference)

  p <- p +
    geom_boxplot(
      aes(
        lower = !!bottom,
        upper = !!top,
      ),
      stat = "identity",
      fill = box_colour,
      col = border_colour
    )

  if (use_reduction_text) {
    p <- p +
      geom_text(
        aes(
          label = ifelse(
            scenario %in% only_scenarios,
            sprintf(
              "%s\n%i%s",
              text_main,
              round(100 * (1 - !!bottom/!!reference)),
              "% lower TP"
            ),
            NA
          ),
          y = !!bottom * (!!top /!!bottom) ^ 0.5  # (midpoint on log scale!)
        ),
        size = text_size,
        col = text_colour
      )
  } else {
    p <- p +
      geom_text(
        aes(
          label = ifelse(
            scenario %in% only_scenarios,
            text_main,
            NA
          ),
          y = !!bottom * (!!top /!!bottom) ^ 0.5  # (midpoint on log scale!)
        ),
        size = text_size,
        col = text_colour
      )
  }

  p

}

# boxplot-y figure on intervention efects on R effective
library(tidyverse)

colours <- RColorBrewer::brewer.pal(3, "Set2")

baseline_colour <- washout(colours[2], 0.8)
vaccine_colour <- washout(colours[3], 0.5)
vaccine_dark_colour <- washout(colours[3], 0.3)
phsm_colours <- washout(colours[1], c(0.5, 0.25, 0.1))
phsm_dark_colours <- washout(colours[1], c(0.4, 0.15, 0))

border_colour <- grey(0.6)
r0_colour <- grey(0.5)
label_colour <- grey(0.3)

text_size <- 2.5

# set up TP outputs for plotting
tps_plot <- tps %>%
  # summarise table as output
  filter(
    variant == "Omicron",
    scenario %in% c(1, 2, 13, 14)
  )


# make plots with both optimal and partial TTIQ
for (ascertainment_plot in c(0, 0.5, 1)) {
  for (omicron_scenario_plot in c("intermediate", "optimistic", "pessimistic")) {

    tps_this_plot <- tps_plot %>%
      filter(
        ascertainment == ascertainment_plot,
        omicron_scenario == omicron_scenario_plot
      ) %>%
      mutate(
        scenario = factor(
          scenario,
          levels = c(1, 2, 13, 14),
          labels = c(
            "6 months\n60% uptake",
            "6 months\n80% uptake",
            "3 months\n60% uptake",
            "3 months\n80% uptake"
          )
        )
      )

    this_r0 <- tps_this_plot$R0[1]

    tp_plot <- tps_this_plot %>%
      control_base_plot() %>%
      add_context_hline(
        label = "Control",
        at = 1,
        linetype = 2,
        text_size = text_size * 1.3
      ) %>%
      # add the baseline + ttiq effect as a box
      add_single_box(
        top = R0,
        bottom = Rt_baseline_ttiq,
        box_colour = baseline_colour,
        only_scenarios = "6 months\n60% uptake",
        text_main = paste0(
          "baseline\nPHSM\n&\n",
          100 * ascertainment_plot,
          "%\nascertainment"
        )
      ) %>%
      add_single_box(
        top = Rt_baseline_ttiq,
        bottom = Rt_baseline_ttiq_vaccination,
        box_colour = vaccine_colour,
        text_main = "boosters",
        only_scenarios = unique(tps_this_plot$scenario),
        use_reduction_text = TRUE
      ) %>%
      add_stacked_box(
        top = Rt_baseline_ttiq_vaccination,
        bottom = Rt_baseline_ttiq_vaccination_low,
        reference = Rt_baseline_ttiq_vaccination,
        text_main = "low\nPHSM",
        only_scenarios = "6 months\n60% uptake",
        box_colour = phsm_colours[1]
      ) %>%
      add_stacked_box(
        top = Rt_baseline_ttiq_vaccination_low,
        bottom = Rt_baseline_ttiq_vaccination_medium,
        reference = Rt_baseline_ttiq_vaccination,
        text_main = "medium\nPHSM",
        only_scenarios = "6 months\n60% uptake",
        box_colour = phsm_colours[2]
      ) %>%
      add_stacked_box(
        top = Rt_baseline_ttiq_vaccination_medium,
        bottom = Rt_baseline_ttiq_vaccination_high,
        reference = Rt_baseline_ttiq_vaccination,
        text_main = "high\nPHSM",
        only_scenarios = "6 months\n60% uptake",
        box_colour = phsm_colours[3]
      ) %>%
      add_context_hline(
        label = "Omicron R0",
        at = this_r0,
        linetype = 2,
        text_size = text_size * 1.3
      ) %>%
      add_context_hline(
        label = "Delta R0",
        at = 8,
        linetype = 2,
        text_size = text_size * 1.3
      ) +
      ggtitle(
        paste("Transmission potential as at", format(target_date, "%B %d")),
        paste(omicron_scenario_plot, "Omicron scenario")
      )

    ggsave(
      paste0(
        "figures/tp_figures/tp_plot_asc",
        ascertainment_plot,
        "_",
        omicron_scenario_plot,
        "_",
        format(target_date, "%Y-%m-%d"),
        ".png"
      ),
      tp_plot,
      width = 8,
      height = 7,
      bg = "white"
    )

  }
}
