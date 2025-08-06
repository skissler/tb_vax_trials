# ==============================================================================
# Import, set parameter values 
# ==============================================================================

library(tidyverse) 
source('code/utils.R')

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

# Define an age distribution (uniform for now): 
agedist <- rep(1, 80)
names(agedist) <- 1:80
agedist <- agedist/sum(agedist)

# ==============================================================================
# Simulate trials under NEJM conditions
# ==============================================================================

trial_df <- sim_trials_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, reps=25)
# trial_df <- read_csv("output/trial_df_nejm.csv")
theoretical_df <- sim_theory_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, agedist=agedist)
fig_trial_theory <- plot_trial_theory(trial_df, theoretical_df) 

# To add: 
# custom age distribution to sim_trial
# custom fast target for sim_trials_over_sigma and sim_theory_over_sigma























# ==============================================================================
# Simulate trials under literature conditions
# ==============================================================================