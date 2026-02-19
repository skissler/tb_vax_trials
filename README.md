# MINA_TB README


This code base `MINA_TB` aims to estimate the required size of
Tuberculosis vaccine trials under different assumptions. Note: this
README file renders to README.md.

## Conventions

### Folders

Code and data are stored in the respective `code/` / `data/` folders.
Model outputs are saved in `outputs/`, usually as .csv files, and
figures in `figures/`. See `notes/` folder for details on how the model
code has developed over time.

### Code

`code/trialsim_rote.R` was Stephen’s main code for running the analysis.
This has been copied to `run_analysis.qmd` in the parent folder for ease
of editing (02/19/2026).

All helper functions are stored in `code/utils.R`. If a piece of code is
re-used more than twice, it should ideally be written up as a function
and added to `utils.R`. Naming conventions for functions are as follows:

`sim_` = run the code

`plot_` = make a plot

`_theory` = theoretical/analytical approach, using probability
expressions derived in `documentation/MINA-TB.pdf`

`_trials` = stochastic simulation-based approach, also described in
`documentation/MINA-TB.pdf`

The probability expressions derived in `documentation/MINA-TB.pdf` are
also saved as the following functions:
`p_inf_and_type_and_asymp_given_age`, `p_asymp_given_age`, and
`p_inf_and_type_given_asymp_and_age`.

Lastly, naming of parameter sets uses the following conventions:

`_lit` = parameter values are taken from the literature

`_nejm` = parameter values are backcalculated from the M72 NEJM trial
(Van Der Meeren et al., 2018)

## Example workflow

The file `code/trialsim_rote.R` has been preserved as an example
workflow. Key steps are as follows:

1.  import functions

``` r
source('code/utils.R')
```

    Warning: package 'forcats' was built under R version 4.5.1

    ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
    ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ✔ purrr     1.0.4     
    ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ✖ dplyr::filter() masks stats::filter()
    ✖ dplyr::lag()    masks stats::lag()
    ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

2.  define parameter values

``` r
pars_nejm_igra <- list(
    minage=18,
    maxage=49,
    rho=0.0275,
    p_slow=0.5,
    mu_slow=0.0001,
    mu_fast=1.5,
    sigma=100,
    ve=0.55,
    trial_length=3)
pars_nejm_tasa <- pars_nejm_igra
pars_nejm_tasa$sigma <- 2

pars_lit_igra <- list(
    minage=18,
    maxage=49,
    rho=0.0025,
    p_slow=0.95,
    mu_slow=0.0001,
    mu_fast=1.5,
    sigma=100,
    ve=0.55,
    trial_length=3)
pars_lit_tasa <- pars_lit_igra 
pars_lit_tasa$sigma <- 2
```

3.  define an age distribution e.g. uniform

``` r
agedist <- rep(1, 80)
names(agedist) <- 1:80
agedist <- agedist/sum(agedist)
```

4.  simulate trials under NEJM conditions (theoretical and stochastic
    approaches)

``` r
#' Simulate the trials approach with NEJM params, over the duration of positivity sigma OR pre-load saved csv file
# trial_df_nejm <- sim_trials_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, reps=25)
trial_df_nejm <- read_csv("output/trial_df_nejm.csv")
```

    Rows: 2000 Columns: 7
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    dbl (7): sigma, rep, n_tested, n_recruited, n_fast, n_slow, n_overall

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#' Simulate the analytical approach with NEJM params, over the duration of positivity sigma
theoretical_df_nejm <- sim_theory_over_sigma(pars=pars_nejm_igra, sigmavec=seq(from=.5, to=80, by=0.01), agedist=agedist)
theoretical_df_nejm %>% filter(sigma %in% c(0.5, 1, 1.5, 2, 5))
```

    # A tibble: 5 × 6
      sigma n_tested n_recruited n_fast n_slow n_overall
      <dbl>    <dbl>       <dbl>  <dbl>  <dbl>     <dbl>
    1   0.5   17720.        121.     50   71.1    25003.
    2   1     11982.        147.     50   96.9    16906.
    3   1.5   10373.        177.     50  127.     14636.
    4   2      9747.        210.     50  160.     13752.
    5   5      9240.        445.     50  395.     13037.

``` r
#' Plot all results
fig_trial_theory_nejm <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, 
                                           cols=c("n_tested","n_recruited","n_fast","n_overall")) + 
    labs(title="M72/AS01E trial emulation")
fig_trial_theory_nejm
```

![](README_files/figure-commonmark/eg_workflow_4-1.png)

``` r
ggsave(fig_trial_theory_nejm, file="figures/trial_theory_nejm.pdf", width=5, height=5/1.6)

fig_trial_theory_log_nejm <- fig_trial_theory_nejm + scale_y_continuous(trans="log10")
fig_trial_theory_log_nejm
```

![](README_files/figure-commonmark/eg_workflow_4-2.png)

``` r
ggsave(fig_trial_theory_log_nejm, file="figures/trial_theory_log_nejm.pdf", width=5, height=5/1.6)

fig_screenslope_nejm <- plot_screenslope(theoretical_df_nejm)
fig_screenslope_nejm
```

![](README_files/figure-commonmark/eg_workflow_4-3.png)

``` r
fig_trial_theory_nejm_r21 <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, cols=c("n_tested","n_recruited")) + 
    labs(title="M72/AS01E trial emulation")
fig_trial_theory_nejm_r21
```

![](README_files/figure-commonmark/eg_workflow_4-4.png)

``` r
ggsave(fig_trial_theory_nejm_r21, file="figures/trial_theory_nejm_r21.pdf", width=5, height=5/1.6)

fig_trial_theory_nejm_r21_log <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, 
                                                   cols=c("n_tested","n_recruited")) + 
    labs(title="M72/AS01E trial emulation") + 
    scale_y_continuous(trans="log10")
fig_trial_theory_nejm_r21_log
```

![](README_files/figure-commonmark/eg_workflow_4-5.png)

``` r
ggsave(fig_trial_theory_nejm_r21_log, file="figures/trial_theory_nejm_r21_log.pdf", width=5, height=5/1.6)
```

5.  simulate trials under literature conditions (theoretical and
    stochastic approaches)

``` r
#' Simulate the trials approach with lit params, over the duration of positivity sigma OR pre-load saved csv file
# trial_df_lit <- sim_trials_over_sigma(pars=pars_lit_igra, sigmavec=1:80, reps=25)
trial_df_lit <- read_csv("output/trial_df_lit.csv")
```

    Rows: 2000 Columns: 7
    ── Column specification ────────────────────────────────────────────────────────
    Delimiter: ","
    dbl (7): sigma, rep, n_tested, n_recruited, n_fast, n_slow, n_overall

    ℹ Use `spec()` to retrieve the full column specification for this data.
    ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
#' Simulate the analytical approach with lit params, over the duration of positivity sigma
theoretical_df_lit <- sim_theory_over_sigma(pars=pars_lit_igra, sigmavec=seq(from=.5, to=80, by=0.01), agedist=agedist)

#' Plot all results
fig_trial_theory_lit <- plot_trial_theory(trial_df_lit, theoretical_df_lit, cols=c("n_tested","n_recruited","n_fast","n_overall")) + 
    labs(title="Literature parameters")
fig_trial_theory_lit
```

![](README_files/figure-commonmark/eg_workflow_5-1.png)

``` r
ggsave(fig_trial_theory_lit, file="figures/trial_theory_lit.pdf", width=5, height=5/1.6)

fig_trial_theory_log_lit <- fig_trial_theory_lit + scale_y_continuous(trans="log10")
fig_trial_theory_log_lit
```

![](README_files/figure-commonmark/eg_workflow_5-2.png)

``` r
ggsave(fig_trial_theory_log_lit, file="figures/trial_theory_log_lit.pdf", width=5, height=5/1.6)

fig_screenslope_lit <- plot_screenslope(theoretical_df_lit)
fig_screenslope_lit
```

![](README_files/figure-commonmark/eg_workflow_5-3.png)

Possible extensions to the above basic workflow are to include custom
age distribution in `sim_trial` and custom fast targets for
`sim_trials_over_sigma` and `sim_theory_over_sigma`.
