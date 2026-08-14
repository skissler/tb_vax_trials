library(shiny)
library(tidyverse)
library(readxl)
library(zoo)

# Shares the same simulation code as the main analysis (run_analysis.qmd) - no separate copy to
# keep in sync. Path is relative to this app's own working directory (shiny/).
source("../code/utils.R")

# Load data once at startup
# unwpp_slim.csv is a pre-filtered version of the UN WPP Excel file containing
# only the countries used by the app — generate it locally with:
# read_excel("data/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx", skip=16) %>%
#   filter(Year==2023 & Type=="Country/Area") %>%
#   rename(country="Region, subregion, country or area *") %>%
#   filter(country %in% c("South Africa","India","Kenya","Tanzania")) %>%
#   mutate(across(all_of(as.character(0:99)), ~as.numeric(.x)*1000)) %>%
#   select(country, all_of(as.character(0:99))) %>%
#   write_csv("shiny/data/unwpp_slim.csv")
unwpp <- read_csv("data/unwpp_slim.csv", show_col_types=FALSE)

progressor_type_by_age_studies <- read.csv("data/progressor-type-by-age-estimates-from-literature.csv") %>%
  numericize(age) %>%
  impute_probability_of_being_slow()

# Named choice vectors for UI
country_choices <- c("South Africa", "India", "Kenya", "Tanzania", "uniform")

# Age-varying FOI options via make_foi(a0_4, a5_24, a25_34, a35_49, a50plus) - includes the
# literature-fitted, paper-final values per country alongside generic flat-rate options.
foi_lookup <- list(
  "South Africa (6,5,7,7,5%, literature-fitted)" = make_foi(6,5,7,7,5),
  "Kenya (5,2,5,4,2%, literature-fitted)"        = make_foi(5,2,5,4,2),
  "India (2,1,3,1,1%, literature-fitted)"        = make_foi(2,1,3,1,1),
  "Uniform 2%"                                   = define_foi_by_age("uniform_2"),
  "Uniform 5%"                                   = define_foi_by_age("uniform_5")
)
foi_choices <- names(foi_lookup)

pslow_choices <- c(
  "5% (uniform)"                       = "uniform",
  "4-14% age-varying (Vynnycky-1997)"  = "Vynnycky-1997",
  "4-19% age-varying (Vynnycky-high)"  = "Vynnycky-high"
)

muslow_choices <- c(
  "0.001 (base)"       = "base",
  "0.0001 (min)"       = "min",
  "0.00527 (max)"      = "max",
  "0.003 (used in paper's main results)" = "0.003"
)

# Reversion parameters - not exposed as inputs; fixed to the values used throughout the paper's
# main results (see pars_lit in run_analysis.qmd).
default_p_revert  <- 0.075
default_mu_revert <- 0.75
