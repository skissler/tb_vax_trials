# FOI fitting history

These CSVs are intermediate results from the age-varying force-of-infection (FOI) fitting process
for South Africa, Kenya, and India — not final pipeline outputs. They record a 5 -> 25 -> 100 rep
re-verification process, run via `code/foi_overnight_run.R` (using the library in
`code/foi_model_fitting_setup.R`).

They're kept for transparency: re-testing at higher rep counts revealed that a candidate FOI
combination which appeared to pass all literature-calibration bounds at 5 reps did not hold up at
100 reps (a "winner's curse" effect from selecting a candidate out of a large pool using a single
noisy statistic). South Africa's original best-fit combination is the clearest example of this,
which is why it required a widened re-search (see the `southafrica_*_history_*` files) rather than
a single stage1 -> stage2 -> stage3 pass.

## File naming

- `<country>_<progression_structure>_<slow_structure>_stage2_reps25.csv` / `stage3_reps100.csv` —
  standard 3-stage pipeline output for a given (country, p_slow, mu_slow) combination, at the
  rep count given in the filename.
- `<country>_<progression_structure>_<slow_structure>_neighbourhood_stage*.csv` — for Kenya/India,
  a small +/-1-per-band neighbourhood search around the already-chosen FOI, run to check it holds
  up under closer scrutiny.
- `<country>_<progression_structure>_<slow_structure>_widened_stage*.csv` — run when a
  neighbourhood search came up empty; widens back out to a full near-miss search across the whole
  candidate grid.
- `southafrica_Vynnycky-high_<slow_structure>_history_*.csv` — South Africa's re-verification
  history specifically, spanning several rounds of increasingly-targeted re-testing as the
  winner's-curse issue was diagnosed (see `code/foi_model_fitting_setup.R`'s `load_test_history()`,
  which reads all files matching this pattern back in and de-duplicates to the highest rep count
  available per candidate).

The final, chosen FOI values for each country are in `output/final_foi_chosen.csv` and
`output/final_foi_recommendations.txt` (one level up), not here.
