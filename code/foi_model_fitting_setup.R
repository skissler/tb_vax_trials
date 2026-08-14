# Age-varying FOI parameter fitting for South Africa, Kenya, and India, as a 3-stage process:
#   Stage 1 (reps=5):   broad grid search across all candidate FOI structures - cheap, coarse.
#   Stage 2 (reps=25):  re-run near-miss candidates from stage 1 (widened bounds, to avoid excluding
#                        true positives lost to stage-1 noise) - narrows down to a much smaller set.
#   Stage 3 (reps=100): re-run whatever passes stage 2 for a high-precision final check. Candidates
#                        that look like they pass at low reps often don't survive this stage - see
#                        South Africa's Vynnycky-high/mu_s=0.003 case, where a 25-rep "63 passing"
#                        result collapsed to 0 at 100 reps (selection bias from filtering on a noisy
#                        statistic - see the GitHub issue writeup for the full story).
#
# Run stage 1 anywhere. Stages 2 and 3 are the expensive/slow ones and are meant to be run
# separately (e.g. on a machine that can be left running) - each stage writes its own CSV, so the
# process can be picked up in a later session as long as the previous stage's CSV is present.

library(tidyverse)
library(readxl)
library(future.apply)
library(progressr)
source('code/utils.R')

plan(multisession, workers=8)

## ---------------------------------------------------------------------------
## Shared setup
## ---------------------------------------------------------------------------

unwpp <- read_excel("data/WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx", skip = 16) %>%
  filter(Year == 2023 & Type == "Country/Area") %>% rename(country="Region, subregion, country or area *") %>%
  mutate(across(all_of(as.character(0:99)), ~ as.numeric(.x) * 1000)) %>% select(country, all_of(as.character(0:99))) %>% arrange(country)
progressor_type_by_age_studies <- read.csv(file="data/progressor-type-by-age-estimates-from-literature.csv") %>%
  numericize(age) %>% impute_probability_of_being_slow()
Dagnew_reversion_data <- left_join(read.csv("data/Dagnew_2026_IGRA_status_at_day_1_table_s4.csv"), read.csv("data/Dagnew_2026_IGRA_reversion_conversion_table_s6.csv"), by = c("country", "site_no" = "site_number")) %>%
  select(country, site_no, N_D1, n_positive_D1, N_M12, reversion_n) %>%
  group_by(country) %>% summarise(across(c(N_D1, n_positive_D1, N_M12, reversion_n), sum)) %>%
  mutate(prob_of_reversion_once_infected = reversion_n / (n_positive_D1*N_M12/N_D1))
Dagnew_positivity_data <- read.csv("data/Dagnew_2026_IGRA_positivity_by_age_table_s5.csv") %>%
  select(country, site_no, age_group, pct_positive, ci_lo_positive, ci_hi_positive) %>%
  mutate(across(c(pct_positive, ci_lo_positive, ci_hi_positive), ~ .x / 100)) %>%
  group_by(country, age_group) %>% summarise(min_positivity = min(ci_lo_positive), max_positivity = max(ci_hi_positive)) %>%
  ungroup()

# Per-country literature bounds. India uses Selvaraju et al 2023 for 15-24/25-34 (overrides Dagnew);
# ARTI/cases bounds as established earlier in the project (see run_analysis.qmd chunk 6_scenarios_setup).
country_bounds <- list(
  "South Africa" = list(
    inf_prev_15to24 = Dagnew_positivity_data %>% filter(country=="South Africa", age_group=="15-24") %>% select(min=min_positivity, max=max_positivity),
    inf_prev_25to34 = Dagnew_positivity_data %>% filter(country=="South Africa", age_group=="25-34") %>% select(min=min_positivity, max=max_positivity),
    arti_lo=0.04, arti_hi=0.06, cases_lo=242, cases_hi=571),
  "Kenya" = list(
    inf_prev_15to24 = Dagnew_positivity_data %>% filter(country=="Kenya", age_group=="15-24") %>% select(min=min_positivity, max=max_positivity),
    inf_prev_25to34 = Dagnew_positivity_data %>% filter(country=="Kenya", age_group=="25-34") %>% select(min=min_positivity, max=max_positivity),
    arti_lo=0.016, arti_hi=0.036, cases_lo=109, cases_hi=339),
  "India" = list(
    inf_prev_15to24 = tibble(min=0.098, max=0.158),  # Selvaraju et al 2023
    inf_prev_25to34 = tibble(min=0.171, max=0.252),  # Selvaraju et al 2023
    arti_lo=0.00, arti_hi=0.02, cases_lo=160, cases_hi=218)
)

# Reversion proxy: Bangladesh used for India (see earlier discussion in run_analysis.qmd)
p_revert_by_country <- list(
  "South Africa" = Dagnew_reversion_data %>% filter(country=="South Africa") %>% pull(prob_of_reversion_once_infected),
  "Kenya"        = Dagnew_reversion_data %>% filter(country=="Kenya")        %>% pull(prob_of_reversion_once_infected),
  "India"        = Dagnew_reversion_data %>% filter(country=="Bangladesh")   %>% pull(prob_of_reversion_once_infected))

agedist_by_country <- setNames(
  lapply(names(country_bounds), function(c) extract_age_distribution(unwpp, c)),
  names(country_bounds))

build_base_pars <- function(country, progression_structure, slow_structure) {
  list(minage=15, maxage=49, rho=NULL,
       p_slow=extract_probability_of_being_slow(progressor_type_by_age_studies, progression_structure),
       mu_slow=define_mu_slow(slow_structure),
       mu_fast=1.5, p_revert=p_revert_by_country[[country]], mu_revert=0.75,
       sigma=100, trial_length=3)
}

parse_foi_structure <- function(foi_structure) {
  parts <- as.numeric(strsplit(foi_structure, "_")[[1]])
  make_foi(parts[1], parts[2], parts[3], parts[4], parts[5])
}

# Runs sim_stoch_over_sigma for each row of foi_grid (a data.frame with a0_4/a5_24/a25_34/a35_49/a50plus
# columns, as from make_foi()'s grid mode) or a plain named list of FOI vectors, at the given rep count.
# Returns one summary row per FOI structure (median inf_prev/ARTI/cases across reps, at sig=80).
run_fit_grid <- function(country, pars_base, foi_structures, agedist, reps, sig=80, casetarget=50) {
  grid <- expand.grid(foi_structure=names(foi_structures), stringsAsFactors=FALSE)
  results <- future_lapply(seq_len(nrow(grid)), function(i) {
    pars_i <- pars_base
    pars_i$rho <- foi_structures[[grid$foi_structure[i]]]
    stochastic_df <- sim_stoch_over_sigma(pars=pars_i, sigmavec=sig, casetarget=casetarget, agedist=agedist, reps=reps)
    from_model <- summarise_case_incidence_and_ARTI(stochastic_df, pars=pars_i, reps=reps, sig=sig)
    tibble(country=country, foi_structure=grid$foi_structure[i],
           inf_prev_15to24=median(stochastic_df$inf_prev_15to24),
           inf_prev_25to34=median(stochastic_df$inf_prev_25to34),
           ARTI=median(from_model$ARTI), cases=median(from_model$cases))
  }, future.seed=TRUE)
  bind_rows(results)
}

# Flags conditions_met using country_bounds; call after run_fit_grid()
apply_bounds <- function(df, country) {
  b <- country_bounds[[country]]
  df %>% mutate(
    conditions_met = between(inf_prev_15to24, b$inf_prev_15to24$min, b$inf_prev_15to24$max) &
      between(inf_prev_25to34, b$inf_prev_25to34$min, b$inf_prev_25to34$max) &
      between(ARTI, b$arti_lo, b$arti_hi) &
      between(cases, b$cases_lo, b$cases_hi))
}

# Loads any previously-saved higher-rep results for a (country, progression_structure, slow_structure)
# combination, so a later run can skip re-testing candidates already confirmed one way or the other.
# Looks for files matching output/foi_fitting_history/{country}_{progression_structure}_{slow_structure}_history_*.csv
# (see foi_overnight_run.R for how these get created) - handles the differing column-naming
# conventions across ad-hoc analysis scripts (some use "_med" suffixes, some use plain names; some
# use "conditions_met", others "conditions_met_25reps"/"conditions_met_100reps"). Returns one row per
# foi_structure, keeping only the highest reps_used available for each. Returns a 0-row tibble (not
# an error) if nothing's been saved yet for this combination.
load_test_history <- function(country, progression_structure, slow_structure) {
  prefix <- sprintf("output/foi_fitting_history/%s_%s_%s_history_", tolower(gsub(" ", "", country)), progression_structure, slow_structure)
  files <- Sys.glob(paste0(prefix, "*.csv"))
  if (length(files) == 0) return(tibble(foi_structure=character(), inf_prev_15to24=numeric(),
    inf_prev_25to34=numeric(), ARTI=numeric(), cases=numeric(), conditions_met=logical(), reps_used=integer(),
    country=character(), progression_structure=character(), slow_structure=character()))

  reps_from_filename <- function(f) as.integer(sub(".*reps(\\d+).*", "\\1", f))

  standardise_one <- function(f) {
    df <- read.csv(f)
    names(df) <- sub("_med$", "", names(df))                 # drop "_med" suffix if present
    names(df)[grepl("^conditions_met", names(df))] <- "conditions_met"  # unify conditions_met* -> conditions_met
    reps <- reps_from_filename(f)
    # .env$ forces these to come from the function's own arguments, not from same-named columns some
    # of these files already carry (e.g. history_full_reps100.csv has its own slow_structure column,
    # written/read back as numeric - without .env$, dplyr's data-masking picks that up instead of the
    # character argument here, causing a bind_rows() type mismatch across files below).
    df %>% mutate(reps_used = reps, country = .env$country, progression_structure = .env$progression_structure,
                  slow_structure = .env$slow_structure) %>%
      select(foi_structure, inf_prev_15to24, inf_prev_25to34, ARTI, cases, conditions_met, reps_used,
             country, progression_structure, slow_structure)
  }

  bind_rows(lapply(files, standardise_one)) %>%
    arrange(foi_structure, desc(reps_used)) %>%
    distinct(foi_structure, .keep_all=TRUE)
}

# Widened near-miss window used to select stage-1 (or stage-2) survivors to re-test at higher reps -
# generous margins to avoid excluding true positives lost to noise at the lower rep count (see the
# South Africa case: candidates well outside a *strict* pass at 5 reps turned out to pass at 25 reps).
near_miss_candidates <- function(df, country, cases_margin=c(-60,60), inf_prev_25to34_margin=c(-0.10,0.10)) {
  b <- country_bounds[[country]]
  df %>% filter(
    cases >= b$cases_lo + cases_margin[1], cases <= b$cases_hi + cases_margin[2],
    inf_prev_25to34 >= b$inf_prev_25to34$min + inf_prev_25to34_margin[1],
    inf_prev_25to34 <= b$inf_prev_25to34$max + inf_prev_25to34_margin[2],
    between(inf_prev_15to24, b$inf_prev_15to24$min, b$inf_prev_15to24$max),
    between(ARTI, b$arti_lo, b$arti_hi))
}


## ---------------------------------------------------------------------------
## Stage 1: broad grid search (reps=5) - safe to run anywhere, including this session
## ---------------------------------------------------------------------------

# FOI grids per country. South Africa searches the widest range; Kenya/India ranges are narrowed
# from an earlier uniform-FOI sweep (see run_analysis.qmd chunk 6a for that sweep).
foi_grid_by_country <- list(
  "South Africa" = bind_rows(
    tibble(a0_4=1:10, a5_24=1:10, a25_34=1:10, a35_49=1:10, a50plus=1:10),
    expand.grid(a0_4=3:7, a5_24=3:7, a25_34=3:7, a35_49=3:7, a50plus=3:7)) %>% distinct(),
  "Kenya" = bind_rows(
    tibble(a0_4=1:10, a5_24=1:10, a25_34=1:10, a35_49=1:10, a50plus=1:10),
    expand.grid(a0_4=2:5, a5_24=2:5, a25_34=2:5, a35_49=2:5, a50plus=2:5)) %>% distinct(),
  "India" = bind_rows(
    tibble(a0_4=1:10, a5_24=1:10, a25_34=1:10, a35_49=1:10, a50plus=1:10),
    expand.grid(a0_4=1:4, a5_24=1:4, a25_34=1:4, a35_49=1:4, a50plus=1:4)) %>% distinct()
)

# Progression/slow structures now narrowed to just Vynnycky-high / mu_s in {"0.003", "max"} for South
# Africa, per the finding that lower mu_s (or Vynnycky-1997, which has a lower fast-progression
# fraction) only reduces case counts further - the opposite of what's needed. Adjust per country as
# appropriate before running.
run_stage1 <- function(country, progression_structure, slow_structure, reps=5) {
  pars_base <- build_base_pars(country, progression_structure, slow_structure)
  foi_structures <- make_foi(foi_grid_by_country[[country]])
  agedist <- agedist_by_country[[country]]

  results <- run_fit_grid(country, pars_base, foi_structures, agedist, reps=reps) %>%
    mutate(progression_structure=progression_structure, slow_structure=slow_structure) %>%
    apply_bounds(country)

  outfile <- sprintf("output/foi_fitting_history/%s_%s_%s_stage1_reps%d.csv",
                      tolower(gsub(" ", "", country)), progression_structure, slow_structure, reps)
  write.csv(results, file=outfile, row.names=FALSE)
  cat(country, progression_structure, slow_structure, "- stage 1 (reps=", reps, "): ",
      sum(results$conditions_met), "of", nrow(results), "passing. Saved to", outfile, "\n")
  results
}

# Loads an already-completed stage-1 (reps=5) result straight from the existing master CSVs, instead
# of re-running run_stage1() from scratch. South Africa's grid search (chunk 4 in run_analysis.qmd)
# already covers all 4 slow_structures x 2 progression_structures; Kenya/India's (chunk 6a) covers
# 2 progression_structures x 3 slow_structures. Check the printed n before assuming a combination is
# covered - if it's 0 rows, that combination genuinely hasn't been run and needs run_stage1() instead.
load_existing_stage1 <- function(country, progression_structure, slow_structure) {
  master_file <- if (country == "South Africa") "output/parameter_fit_results.csv" else "output/kenya_india_parameter_fit_results.csv"
  df <- read.csv(master_file) %>%
    filter(progression_structure == !!progression_structure, slow_structure == !!slow_structure)
  if (!"country" %in% names(df)) df$country <- country  # South Africa's master file has no country column (SA-only)
  df <- df %>% filter(country == !!country)
  if (!"reps_used" %in% names(df)) df$reps_used <- 5
  cat(country, progression_structure, slow_structure, "- loaded", nrow(df), "existing stage-1 rows from", master_file, "\n")
  df %>% select(country, foi_structure, inf_prev_15to24, inf_prev_25to34, ARTI, cases, conditions_met, reps_used) %>%
    mutate(progression_structure=progression_structure, slow_structure=slow_structure)
}

# Example (run when ready):
# sa_max_stage1   <- load_existing_stage1("South Africa", "Vynnycky-high", "max")     # already have this - no need to rerun
# sa_0003_stage1  <- load_existing_stage1("South Africa", "Vynnycky-high", "0.003")   # already have this too


## ---------------------------------------------------------------------------
## Stage 2: near-miss re-run at reps=25 - slower, intended to be run separately
## ---------------------------------------------------------------------------

run_stage2 <- function(country, progression_structure, slow_structure, stage1_results, reps=25, ...) {
  pars_base <- build_base_pars(country, progression_structure, slow_structure)
  agedist <- agedist_by_country[[country]]

  candidates <- near_miss_candidates(stage1_results, country, ...)
  cat(country, progression_structure, slow_structure, "- stage 2: re-running", nrow(candidates), "near-miss candidates at reps=", reps, "\n")
  foi_structures <- setNames(lapply(candidates$foi_structure, parse_foi_structure), candidates$foi_structure)

  results <- run_fit_grid(country, pars_base, foi_structures, agedist, reps=reps) %>%
    mutate(progression_structure=progression_structure, slow_structure=slow_structure) %>%
    apply_bounds(country)

  outfile <- sprintf("output/foi_fitting_history/%s_%s_%s_stage2_reps%d.csv",
                      tolower(gsub(" ", "", country)), progression_structure, slow_structure, reps)
  write.csv(results, file=outfile, row.names=FALSE)
  cat(country, progression_structure, slow_structure, "- stage 2 (reps=", reps, "): ",
      sum(results$conditions_met), "of", nrow(results), "passing. Saved to", outfile, "\n")
  results
}

# Example (run separately - this is the expensive stage):
# sa_max_stage1  <- read.csv("output/foi_fitting_history/southafrica_Vynnycky-high_max_stage1_reps5.csv")
# sa_max_stage2  <- run_stage2("South Africa", "Vynnycky-high", "max", sa_max_stage1, reps=25)


## ---------------------------------------------------------------------------
## Stage 3: high-precision confirmation at reps=100 - slowest, run separately.
## Re-runs stage 2's passing (or, better, all near-miss) candidates - passing candidates from stage 2
## are NOT guaranteed to hold up here, per the South Africa Vynnycky-high/0.003 experience.
## ---------------------------------------------------------------------------

run_stage3 <- function(country, progression_structure, slow_structure, stage2_results, reps=100, only_passing=FALSE) {
  pars_base <- build_base_pars(country, progression_structure, slow_structure)
  agedist <- agedist_by_country[[country]]

  candidates <- if (only_passing) stage2_results %>% filter(conditions_met) else stage2_results
  cat(country, progression_structure, slow_structure, "- stage 3: re-running", nrow(candidates), "candidates at reps=", reps, "\n")
  foi_structures <- setNames(lapply(candidates$foi_structure, parse_foi_structure), candidates$foi_structure)

  results <- run_fit_grid(country, pars_base, foi_structures, agedist, reps=reps) %>%
    mutate(progression_structure=progression_structure, slow_structure=slow_structure) %>%
    apply_bounds(country)

  outfile <- sprintf("output/foi_fitting_history/%s_%s_%s_stage3_reps%d.csv",
                      tolower(gsub(" ", "", country)), progression_structure, slow_structure, reps)
  write.csv(results, file=outfile, row.names=FALSE)
  cat(country, progression_structure, slow_structure, "- stage 3 (reps=", reps, "): ",
      sum(results$conditions_met), "of", nrow(results), "passing. Saved to", outfile, "\n")
  results
}

# Example (run separately - the slowest stage):
# sa_max_stage2  <- read.csv("output/foi_fitting_history/southafrica_Vynnycky-high_max_stage2_reps25.csv")
# sa_max_stage3  <- run_stage3("South Africa", "Vynnycky-high", "max", sa_max_stage2, reps=100)


## ---------------------------------------------------------------------------
## Merge a completed stage's results back into the master parameter_fit_results.csv, tagging each
## updated row with which stage/rep count it reflects. Run after stage 2 or 3 completes.
## ---------------------------------------------------------------------------

# Per-row version of apply_bounds() - looks up each row's own country in country_bounds instead of
# applying a single country's bounds to the whole data frame. Needed for master files that hold more
# than one country (e.g. kenya_india_parameter_fit_results.csv) - using the single-country apply_bounds()
# on those would silently recompute conditions_met for every row using the wrong country's bounds.
apply_bounds_per_row <- function(df) {
  df %>% rowwise() %>% mutate(conditions_met = {
    b <- country_bounds[[country]]
    between(inf_prev_15to24, b$inf_prev_15to24$min, b$inf_prev_15to24$max) &
      between(inf_prev_25to34, b$inf_prev_25to34$min, b$inf_prev_25to34$max) &
      between(ARTI, b$arti_lo, b$arti_hi) &
      between(cases, b$cases_lo, b$cases_hi)
  }) %>% ungroup()
}

merge_into_master <- function(new_results, reps_used, master_file="output/parameter_fit_results.csv") {
  master <- read.csv(master_file)
  if (!"reps_used" %in% names(master)) master$reps_used <- 5
  if (!"country" %in% names(master)) master$country <- new_results$country[1]  # SA-only master file has no country column

  target <- master$country == new_results$country[1] &
    master$progression_structure == new_results$progression_structure[1] &
    master$slow_structure == new_results$slow_structure[1] &
    master$foi_structure %in% new_results$foi_structure
  match_idx <- match(master$foi_structure[target], new_results$foi_structure)

  idx <- which(target)
  master$inf_prev_15to24[idx] <- new_results$inf_prev_15to24[match_idx]
  master$inf_prev_25to34[idx] <- new_results$inf_prev_25to34[match_idx]
  master$ARTI[idx]            <- new_results$ARTI[match_idx]
  master$cases[idx]           <- new_results$cases[match_idx]
  master$reps_used[idx]       <- reps_used

  # Recompute conditions_met per row using each row's own country's bounds - safe for both the
  # single-country South Africa file and the multi-country Kenya/India file.
  master <- apply_bounds_per_row(master)

  write.csv(master, file=master_file, row.names=FALSE)
  cat("Merged", sum(target), "rows (reps_used=", reps_used, ") into", master_file, "\n")
  master
}

# Example:
# merge_into_master(sa_max_stage3, reps_used=100)
