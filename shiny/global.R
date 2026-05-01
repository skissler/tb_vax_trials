library(shiny)
library(tidyverse)
library(readxl)
library(zoo)

source("utils.R")

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

foi_choices <- c(
  "1% (Low burden)"                             = "uniform1",
  "2% (Low burden)"                             = "uniform2",
  "7% (High burden baseline)"                   = "base",
  "5% up to age 25, then 1% (High burden min)"  = "min",
  "10% (High burden max)"                       = "max"
)

pslow_choices <- c(
  "5% (uniform)"                    = "uniform",
  "4-14% age-varying (Vynnycky-1997)" = "Vynnycky-1997"
)

muslow_choices <- c(
  "0.001 (base)" = "base",
  "0.0001 (min)" = "min",
  "0.00527 (max)" = "max"
)
