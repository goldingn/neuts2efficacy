#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author Nick Golding
#' @export
make_ngm <- function() {

  lookups <- get_quantium_lookups()

  age_breaks <- lookups$age$age_band %>%
    str_split_fixed("-", 2) %>%
    `[`(, 1) %>%
    str_remove("\\+") %>%
    as.numeric() %>%
    c(Inf)

  # fit polymod model
  setting_models <- fit_setting_contacts(
    contact_data_list = get_polymod_setting_data(),
    population = get_polymod_population()
  )

  australia_pop <- abs_pop_age_lga_2020 %>%
    group_by(age_group) %>%
    summarise(
      population = sum(population)
    ) %>%
    mutate(
      lower.age.limit = readr::parse_number(as.character(age_group))
    ) %>%
    mutate(
      country = "Australia"
    )

  # get Australia-wide setting-specific contact matrices
  # with adjustment to the household size
  australia_setting_contact_matrices <- predict_setting_contacts(
    contact_model = setting_models,
    population = australia_pop,
    per_capita_household_size = get_per_capita_household_size(),
    age_breaks = age_breaks
  )

  # get the setting-specific contact matrices
  setting_transmission_matrices <- get_setting_transmission_matrices(
    age_breaks = age_breaks,
    asymptomatic_relative_infectiousness = 0.5,
  )

  settings <- c("home", "school", "work", "other")

  # setting specific unscaled NGMs
  australia_setting_ngms <- mapply(
    "*",
    australia_setting_contact_matrices[settings],
    setting_transmission_matrices[settings],
    SIMPLIFY = FALSE
  )

  # overall NGM
  australia_ngm <- Reduce("+", australia_setting_ngms)

  australia_ngm

}
