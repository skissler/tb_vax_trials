# Standalone overnight run: finds robust FOI choices for South Africa, Kenya, and India.
# Run this directly, e.g. from a terminal:
#   "C:/Program Files/R/R-4.5.0/bin/Rscript.exe" code/foi_overnight_run.R
# Not tied to any chat session - safe to close everything else, just keep this process and the
# machine running (disable sleep/hibernate; a laptop lid closing will suspend R too).
#
# Progress is logged to output/overnight_run_log.txt with timestamps - check that file if you want
# to see how far it got before it's finished. Final answers are written to
# output/final_foi_recommendations.csv (machine-readable) and .txt (human-readable summary).
#
# What this does and why:
#   - South Africa (Vynnycky-high / mu_s=0.003 AND max) - for each, reuses any saved history (see the
#     output/southafrica_Vynnycky-high_*_history_*.csv files) to skip candidates already confirmed at
#     25 or 100 reps, and only fills in the gaps: near-miss candidates never tested at 100 reps yet.
#     For 0.003 this means most of the 1115-candidate near-miss pool (already at 25 reps) just needs
#     the ~9 untested-at-100-reps candidates filled in; for max, nothing's been tested beyond 5 reps
#     yet, so it runs the full stage1(reuse)->stage2(25 reps)->stage3(100 reps) pipeline from scratch.
#   - Kenya and India - their currently-chosen FOI structures were selected using different metrics
#     (median-closeness for Kenya, trend-correlation for India), not a simple pass/fail bounds check,
#     so they're less exposed to the sharp-threshold selection bias that hit South Africa's 0.003 case.
#     Still, neither has been checked beyond 5 reps, so this verifies the current choice plus a small
#     neighbourhood (+/-1 on each FOI band) at 25 then 100 reps, and reports whether the choice holds
#     or a neighbour turns out better under scrutiny.

start_time <- Sys.time()
log_file <- "output/overnight_run_log.txt"
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", paste0(..., collapse=""))
  cat(msg, "\n")
  cat(msg, "\n", file=log_file, append=TRUE)
}
unlink(log_file)
log_msg("Overnight run started")

source('code/foi_model_fitting.R')
log_msg("Sourced foi_model_fitting.R")

final_recommendations <- list()

## ---------------------------------------------------------------------------
## South Africa: Vynnycky-high / mu_s in {0.003, max}
## ---------------------------------------------------------------------------

wood <- read.csv("data/Wood-prevalence-plot-data_plotdigitizer.csv")
names(wood) <- c("age","prevalence")
wood <- wood %>% distinct(age, .keep_all=TRUE) %>% arrange(age)
compute_prevalence_curve <- function(rho, p_revert, mu_revert, ages=0:60) {
  pmf <- compute_survival_fn(rho)$survival_pmf[1:100]
  sapply(ages, function(a) {
    if (a < 1) return(0)
    t <- 0:(a-1)
    p_not_reverted <- 1 - p_revert*(1 - exp(-mu_revert*(a - t)))
    sum(pmf[t+1] * p_not_reverted)
  })
}

sa_pools <- list()

for (ms in c("0.003", "max")) {
  tryCatch({
    log_msg("=== South Africa (Vynnycky-high, ", ms, ") ===")
    history <- load_test_history("South Africa", "Vynnycky-high", ms)

    if (nrow(history) > 0 && any(history$reps_used >= 25)) {
      log_msg("Found existing history: ", nrow(history), " candidates, reps up to ", max(history$reps_used))
      stage_low <- history %>% filter(reps_used >= 25)
      near_miss <- near_miss_candidates(stage_low, "South Africa")
      already_100 <- history %>% filter(reps_used >= 100) %>% pull(foi_structure)
      to_test <- near_miss %>% filter(!foi_structure %in% already_100)
      log_msg("Near-miss pool: ", nrow(near_miss), "; already confirmed at reps=100: ", length(already_100),
              "; remaining to test at reps=100: ", nrow(to_test))
    } else {
      log_msg("No existing higher-rep history - starting from stage 1 (reps=5)")
      stage1 <- load_existing_stage1("South Africa", "Vynnycky-high", ms)
      log_msg("Stage 1: ", sum(stage1$conditions_met), " of ", nrow(stage1), " pass at reps=5")
      near_miss_5 <- near_miss_candidates(stage1, "South Africa")
      log_msg("Stage 2: running ", nrow(near_miss_5), " near-miss candidates at reps=25")
      stage2 <- run_stage2("South Africa", "Vynnycky-high", ms, stage1, reps=25)
      log_msg("Stage 2 done: ", sum(stage2$conditions_met), " of ", nrow(stage2), " pass at reps=25")
      to_test <- near_miss_candidates(stage2, "South Africa")
      already_100 <- character(0)
      log_msg("Stage 3: all ", nrow(to_test), " near-miss candidates need reps=100 (none tested yet)")
    }

    if (nrow(to_test) > 0) {
      stage3_new <- run_stage3("South Africa", "Vynnycky-high", ms, to_test, reps=100, only_passing=FALSE) %>%
        mutate(reps_used=100)
      log_msg("Newly tested at reps=100: ", sum(stage3_new$conditions_met), " of ", nrow(stage3_new), " pass")
    } else {
      log_msg("Nothing new to test - all near-miss candidates already confirmed at reps=100")
      stage3_new <- tibble()
    }

    already_100_data <- if (nrow(history) > 0) history %>% filter(reps_used >= 100) %>%
      mutate(progression_structure="Vynnycky-high", slow_structure=ms) else tibble()
    stage3_all <- bind_rows(already_100_data, stage3_new) %>% distinct(foi_structure, .keep_all=TRUE)

    log_msg("Full reps=100 picture for mu_s=", ms, ": ", sum(stage3_all$conditions_met), " of ", nrow(stage3_all), " pass")
    outfile <- sprintf("output/southafrica_Vynnycky-high_%s_history_full_reps100.csv", ms)
    write.csv(stage3_all, file=outfile, row.names=FALSE)

    sa_pools[[ms]] <- stage3_all
  }, error = function(e) log_msg("ERROR in South Africa (", ms, ") block: ", conditionMessage(e)))
}

tryCatch({
  sa_all <- bind_rows(sa_pools)
  sa_passing <- sa_all %>% filter(conditions_met)
  sa_status <- "ROBUST PASS"
  if (nrow(sa_passing) == 0) {
    log_msg("Nothing passes at reps=100 in either mu_s=0.003 or max - falling back to closest-miss candidate")
    sa_passing <- sa_all %>%
      mutate(cases_gap = pmax(242 - cases, 0), inf_gap = pmax(inf_prev_25to34 - 0.751, 0)) %>%
      arrange(cases_gap + inf_gap * 10)
    sa_status <- "NO ROBUST PASS - closest miss reported"
  }

  # Fit to Wood et al 2010 curve among whatever's in the final pool, to pick the best of them
  sa_passing <- sa_passing %>% rowwise() %>%
    mutate(rmse = {
      rho <- parse_foi_structure(foi_structure)
      curve <- compute_prevalence_curve(rho, p_revert_by_country[["South Africa"]], 0.75, ages=0:60)
      pred <- approx(x=0:60, y=curve, xout=wood$age)$y
      sqrt(mean((pred - wood$prevalence)^2))
    }) %>% ungroup() %>% arrange(rmse)

  sa_best <- sa_passing %>% slice(1)
  log_msg("South Africa best: ", sa_best$foi_structure, " (mu_s=", sa_best$slow_structure, ", ", sa_status,
          "), RMSE=", round(sa_best$rmse,4), ", cases=", round(sa_best$cases,1), ", ARTI=", round(sa_best$ARTI,4),
          ", inf_prev_15to24=", round(sa_best$inf_prev_15to24,3), ", inf_prev_25to34=", round(sa_best$inf_prev_25to34,3))

  final_recommendations[["South Africa"]] <- sa_best %>%
    mutate(country="South Africa", progression_structure="Vynnycky-high", status=sa_status) %>%
    select(country, progression_structure, slow_structure, foi_structure, status,
           inf_prev_15to24, inf_prev_25to34, ARTI, cases, any_of("rmse"))

}, error = function(e) log_msg("ERROR in South Africa final selection: ", conditionMessage(e)))


## ---------------------------------------------------------------------------
## Kenya and India: verify already-chosen FOI + small neighbourhood
## ---------------------------------------------------------------------------

chosen_foi <- list(
  "Kenya" = list("Vynnycky-1997" = c(5,2,5,4,2), "uniform" = c(2,3,5,4,2)),
  "India" = list("Vynnycky-1997" = c(1,2,4,2,1), "uniform" = c(1,3,4,2,1))
)

# Small neighbourhood: +/-1 on each band around the chosen point (3^5 = 243 combos, minus any with
# non-positive values), capped to keep runtime bounded on top of South Africa's search above.
build_neighbourhood <- function(centre) {
  bands <- lapply(centre, function(v) unique(pmax(v + c(-1,0,1), 1)))
  grid <- expand.grid(bands)
  names(grid) <- c("a0_4","a5_24","a25_34","a35_49","a50plus")
  grid
}

for (country in c("Kenya", "India")) {
  for (ps in names(chosen_foi[[country]])) {
    tryCatch({
      log_msg("=== ", country, " (", ps, ", mu_s=0.003) ===")
      centre <- chosen_foi[[country]][[ps]]
      centre_label <- paste(centre, collapse="_")
      log_msg("Chosen FOI: ", centre_label, " - verifying plus a +/-1 neighbourhood")

      pars_base <- build_base_pars(country, ps, "0.003")
      agedist <- agedist_by_country[[country]]
      neigh_grid <- build_neighbourhood(centre)
      foi_structures <- make_foi(neigh_grid)

      log_msg("Stage 2: running ", length(foi_structures), " neighbourhood candidates at reps=25")
      stage2 <- run_fit_grid(country, pars_base, foi_structures, agedist, reps=25) %>%
        mutate(progression_structure=ps, slow_structure="0.003") %>%
        apply_bounds(country)
      outfile2 <- sprintf("output/%s_%s_0.003_neighbourhood_stage2_reps25.csv", tolower(country), ps)
      write.csv(stage2, file=outfile2, row.names=FALSE)
      log_msg("Stage 2 done: ", sum(stage2$conditions_met), " of ", nrow(stage2), " pass at reps=25. Saved to ", outfile2)

      # Confirm ALL neighbourhood candidates at reps=100 (small set, so no need to pre-filter)
      log_msg("Stage 3: running all ", nrow(stage2), " neighbourhood candidates at reps=100")
      stage3 <- run_fit_grid(country, pars_base, foi_structures, agedist, reps=100) %>%
        mutate(progression_structure=ps, slow_structure="0.003") %>%
        apply_bounds(country)
      outfile3 <- sprintf("output/%s_%s_0.003_neighbourhood_stage3_reps100.csv", tolower(country), ps)
      write.csv(stage3, file=outfile3, row.names=FALSE)
      log_msg("Stage 3 done: ", sum(stage3$conditions_met), " of ", nrow(stage3), " pass at reps=100. Saved to ", outfile3)

      centre_row <- stage3 %>% filter(foi_structure == centre_label)
      centre_holds <- nrow(centre_row) > 0 && centre_row$conditions_met[1]

      cases_mid <- (country_bounds[[country]]$cases_lo + country_bounds[[country]]$cases_hi) / 2
      passing <- stage3 %>% filter(conditions_met) %>% mutate(cases_dev = abs(cases - cases_mid)) %>% arrange(cases_dev)

      if (centre_holds) {
        status <- "ORIGINAL CHOICE CONFIRMED at reps=100"
        best <- centre_row
      } else if (nrow(passing) > 0) {
        status <- "ORIGINAL CHOICE DID NOT HOLD - neighbour substituted"
        best <- passing %>% slice(1)
      } else {
        # Small neighbourhood came up empty entirely - widen to a full near-miss search across the
        # existing grid, same pattern as South Africa's block above, instead of just reporting a
        # nearby FOI that doesn't actually pass anything.
        log_msg(country, " (", ps, ") neighbourhood fully failed - widening to full-grid near-miss search")
        wide_stage1 <- load_existing_stage1(country, ps, "0.003")
        wide_near_miss <- near_miss_candidates(wide_stage1, country)
        log_msg("Widened stage 2: re-running ", nrow(wide_near_miss), " near-miss candidates at reps=25")
        wide_stage2 <- run_stage2(country, ps, "0.003", wide_stage1, reps=25)
        outfile2w <- sprintf("output/%s_%s_0.003_widened_stage2_reps25.csv", tolower(country), ps)
        write.csv(wide_stage2, file=outfile2w, row.names=FALSE)
        log_msg("Widened stage 2 done: ", sum(wide_stage2$conditions_met), " of ", nrow(wide_stage2), " pass at reps=25")

        wide_survivors <- wide_stage2 %>% filter(conditions_met)
        if (nrow(wide_survivors) == 0) wide_survivors <- near_miss_candidates(wide_stage2, country)
        log_msg("Widened stage 3: re-running ", nrow(wide_survivors), " candidates at reps=100")
        wide_stage3 <- run_stage3(country, ps, "0.003", wide_survivors, reps=100, only_passing=FALSE)
        outfile3w <- sprintf("output/%s_%s_0.003_widened_stage3_reps100.csv", tolower(country), ps)
        write.csv(wide_stage3, file=outfile3w, row.names=FALSE)
        log_msg("Widened stage 3 done: ", sum(wide_stage3$conditions_met), " of ", nrow(wide_stage3), " pass at reps=100")

        wide_passing <- wide_stage3 %>% filter(conditions_met) %>% mutate(cases_dev = abs(cases - cases_mid)) %>% arrange(cases_dev)
        if (nrow(wide_passing) > 0) {
          status <- "NEIGHBOURHOOD FAILED - found robust pass via widened full-grid search"
          best <- wide_passing %>% slice(1)
        } else {
          status <- "NO ROBUST PASS even after widened search - closest miss reported"
          best <- wide_stage3 %>% mutate(cases_dev = abs(cases - cases_mid)) %>% arrange(cases_dev) %>% slice(1)
        }
      }
      log_msg(country, " (", ps, ") result: ", status, " -> ", best$foi_structure,
              ", cases=", round(best$cases,1), ", ARTI=", round(best$ARTI,4),
              ", inf_prev_15to24=", round(best$inf_prev_15to24,3), ", inf_prev_25to34=", round(best$inf_prev_25to34,3))

      final_recommendations[[paste(country, ps)]] <- best %>%
        mutate(country=country, status=status) %>%
        select(country, progression_structure, slow_structure, foi_structure, status,
               inf_prev_15to24, inf_prev_25to34, ARTI, cases)

    }, error = function(e) log_msg("ERROR in ", country, " (", ps, ") block: ", conditionMessage(e)))
  }
}


## ---------------------------------------------------------------------------
## Final summary
## ---------------------------------------------------------------------------

tryCatch({
  final_df <- bind_rows(final_recommendations)
  write.csv(final_df, file="output/final_foi_recommendations.csv", row.names=FALSE)

  summary_lines <- c(
    paste0("FOI fitting overnight run - completed ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("Total runtime: ", round(difftime(Sys.time(), start_time, units="hours"), 2), " hours"),
    "",
    "FINAL RECOMMENDATIONS:",
    ""
  )
  for (i in seq_len(nrow(final_df))) {
    r <- final_df[i,]
    summary_lines <- c(summary_lines, sprintf(
      "%-14s %-16s foi=%-14s [%s]\n  cases=%.1f  ARTI=%.4f  inf_prev_15to24=%.3f  inf_prev_25to34=%.3f\n",
      r$country, r$progression_structure, r$foi_structure, r$status,
      r$cases, r$ARTI, r$inf_prev_15to24, r$inf_prev_25to34))
  }
  writeLines(summary_lines, "output/final_foi_recommendations.txt")
  log_msg("Final recommendations written to output/final_foi_recommendations.csv and .txt")
  log_msg("=== OVERNIGHT RUN COMPLETE ===")
}, error = function(e) log_msg("ERROR writing final summary: ", conditionMessage(e)))
