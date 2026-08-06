# Exploratory analysis used to choose age-varying FOI structures for South Africa, India, and Kenya.
# For each country, candidate FOI structures from the parameter-fitting output (chunk 4 / 6a in
# run_analysis.qmd) are scored against external infection-prevalence-by-age data (where available)
# or against the literature threshold bounds (Kenya). Not part of the main analysis pipeline - run
# interactively, after run_analysis.qmd chunks 1-3 and the relevant parameter_fit chunk(s).

library(tidyverse)
library(readxl)
source('code/utils.R')

unwpp <- read_excel("data/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx", skip = 16) %>%
  filter(Year == 2023 & Type == "Country/Area") %>% rename(country="Region, subregion, country or area *") %>%
  mutate(across(all_of(as.character(0:99)), ~ as.numeric(.x) * 1000)) %>% select(country, all_of(as.character(0:99))) %>% arrange(country)
Dagnew_reversion_data <- left_join(read.csv("data/Dagnew_2026_IGRA_status_at_day_1_table_s4.csv"), read.csv("data/Dagnew_2026_IGRA_reversion_conversion_table_s6.csv"), by = c("country", "site_no" = "site_number")) %>%
  select(country, site_no, N_D1, n_positive_D1, N_M12, reversion_n) %>%
  group_by(country) %>% summarise(across(c(N_D1, n_positive_D1, N_M12, reversion_n), sum)) %>%
  mutate(prob_of_reversion_once_infected = reversion_n / (n_positive_D1*N_M12/N_D1))
mu_revert <- 0.75

# Model-implied current-infection prevalence by age, accounting for reversion:
# prevalence(a) = sum_{t=0}^{a-1} pmf(t) * P(not yet reverted by age a | infected at t)
compute_prevalence_curve <- function(rho, p_revert, mu_revert, ages=0:80) {
  pmf <- compute_survival_fn(rho)$survival_pmf[1:100]
  sapply(ages, function(a) {
    if (a < 1) return(0)
    t <- 0:(a-1)
    p_not_reverted <- 1 - p_revert*(1 - exp(-mu_revert*(a - t)))
    sum(pmf[t+1] * p_not_reverted)
  })
}

# Parse a foi_structure label ("a0_4_a5_24_a25_34_a35_49_a50plus") back into a FOI vector via make_foi()
parse_foi_structure <- function(foi_structure) {
  parts <- as.numeric(strsplit(foi_structure, "_")[[1]])
  make_foi(parts[1], parts[2], parts[3], parts[4], parts[5])
}


## ---------------------------------------------------------------------------
## South Africa: fit against Wood et al. 2010 IGRA prevalence-by-age curve
## ---------------------------------------------------------------------------

p_revert_sa <- Dagnew_reversion_data %>% filter(country=="South Africa") %>% pull(prob_of_reversion_once_infected)
wood <- read.csv("data/Wood-prevalence-plot-data_plotdigitizer.csv")
names(wood) <- c("age","prevalence")
wood <- wood %>% distinct(age, .keep_all=TRUE) %>% arrange(age)

sa_fit <- read.csv("output/parameter_fit_results.csv") %>% filter(conditions_met==TRUE)
sa_unique_foi <- unique(sa_fit$foi_structure)

sa_rmse <- setNames(sapply(sa_unique_foi, function(fs) {
  rho <- parse_foi_structure(fs)
  curve <- compute_prevalence_curve(rho, p_revert_sa, mu_revert, ages=0:60)
  pred <- approx(x=0:60, y=curve, xout=wood$age)$y
  sqrt(mean((pred - wood$prevalence)^2))
}), sa_unique_foi)

sa_results <- tibble(foi_structure=names(sa_rmse), rmse=sa_rmse) %>% arrange(rmse)
print(head(sa_results, 15))
# NB: the 50+ band is unidentifiable from this data (Wood's curve only extends to ~age 40), and none of
# the passing structures reproduce the post-age-30 decline in Wood's curve (reversion too weak relative
# to residual FOI at older ages within the fixed p_revert/mu_revert used here).
# Chosen: 4-6-4-7-5% (0-4/5-24/25-34/35-49/50+)


## ---------------------------------------------------------------------------
## India: fit against Selvaraju et al. 2023 (National TB Prevalence Survey, Table 2, adjusted TBI)
## ---------------------------------------------------------------------------

p_revert_india <- Dagnew_reversion_data %>% filter(country=="Bangladesh") %>% pull(prob_of_reversion_once_infected)  # Bangladesh used as proxy for India

selvaraju <- tibble(
  age_mid    = c(19.5, 29.5, 39.5, 49.5, 59.5, 70),
  prevalence = c(0.128, 0.212, 0.260, 0.293, 0.285, 0.287)  # adjusted (clustering-adjusted) TBI, Table 2
)

ki_fit <- read.csv("output/kenya_india_parameter_fit_results.csv")
india <- ki_fit %>% filter(country=="India", slow_structure=="0.003", between(ARTI, arti_lo, arti_hi))
india_unique_foi <- unique(india$foi_structure)

india_metrics <- lapply(india_unique_foi, function(fs) {
  rho <- parse_foi_structure(fs)
  curve <- compute_prevalence_curve(rho, p_revert_india, mu_revert, ages=0:80)
  pred <- approx(x=0:80, y=curve, xout=selvaraju$age_mid)$y
  tibble(foi_structure=fs, rmse=sqrt(mean((pred - selvaraju$prevalence)^2)), trend_cor=cor(pred, selvaraju$prevalence))
}) %>% bind_rows()

india <- india %>% left_join(india_metrics, by="foi_structure")

# Cases_lo=160 is incompatible with a good absolute fit (best unconstrained RMSE~0.06 gives only ~69
# cases); since the prevalence survey is more trustworthy than the case-count threshold (which may be
# dominated by high-burden pockets this national model doesn't capture), match the age TREND instead,
# with a soft minimum of >=100 cases for trial feasibility.
for (ps in c("Vynnycky-1997", "uniform")) {
  cat("\n---", ps, "---\n")
  print(india %>% filter(progression_structure==ps, cases>=100, cases<130) %>%
          arrange(desc(trend_cor)) %>% head(3) %>%
          select(foi_structure, trend_cor, rmse, cases, inf_prev_15to24, inf_prev_25to34, ARTI))
}
# Chosen: Vynnycky-1997 -> 1-2-4-2-1%, uniform -> 1-3-4-2-1%


## ---------------------------------------------------------------------------
## Kenya: no age-specific prevalence data available - match median of all threshold bounds instead,
## with a trend-plausibility check against the Wood/Selvaraju fits above
## ---------------------------------------------------------------------------

p_revert_kenya <- Dagnew_reversion_data %>% filter(country=="Kenya") %>% pull(prob_of_reversion_once_infected)

kenya <- ki_fit %>% filter(country=="Kenya", slow_structure=="0.003") %>%
  mutate(
    pos_15to24 = (inf_prev_15to24 - min_15to24) / (max_15to24 - min_15to24),
    pos_25to34 = (inf_prev_25to34 - min_25to34) / (max_25to34 - min_25to34),
    pos_arti   = (ARTI - arti_lo) / (arti_hi - arti_lo),
    pos_cases  = (cases - cases_lo) / (cases_hi - cases_lo),
    in_bounds  = between(inf_prev_15to24, min_15to24, max_15to24) &
                 between(inf_prev_25to34, min_25to34, max_25to34) &
                 between(ARTI, arti_lo, arti_hi) &
                 between(cases, cases_lo, cases_hi),
    mean_dev    = (abs(pos_15to24-0.5) + abs(pos_25to34-0.5) + abs(pos_arti-0.5) + abs(pos_cases-0.5))/4,
    trend_ratio = inf_prev_25to34 / inf_prev_15to24,  # sanity check vs Wood (~1.16) / Selvaraju (~1.66)
    a50plus     = as.numeric(sapply(strsplit(foi_structure, "_"), `[`, 5))
  )

# Restrict to the lowest available 50+ band FOI (2%, the floor of the explored grid) so the resulting
# age curve flattens at older ages rather than climbing at a constant rate throughout - a better
# qualitative match to the deceleration seen in both Wood and Selvaraju's older-age data.
for (ps in c("Vynnycky-1997", "uniform")) {
  cat("\n---", ps, "---\n")
  print(kenya %>% filter(progression_structure==ps, in_bounds, a50plus<=2) %>%
          arrange(mean_dev) %>% head(5) %>%
          select(foi_structure, mean_dev, inf_prev_15to24, inf_prev_25to34, trend_ratio, ARTI, cases))
}
# Chosen: Vynnycky-1997 -> 5-2-5-4-2%, uniform -> 2-3-5-4-2%
