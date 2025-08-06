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

# ==============================================================================
# Simulate trials under NEJM conditions
# ==============================================================================

# Define an age distribution (uniform for now): 
agedist <- rep(1, 80)
names(agedist) <- 1:80
agedist <- agedist/sum(agedist)

trial_df <- sim_trials_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, reps=25)
theoretical_df <- sim_theory_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, agedist=agedist)

plot_trial_theory <- function(trial_df, theoretical_df){
	
}


trial_df_toplot <- trial_df %>% 
	select(sigma, n_tested, n_recruited, n_fast) %>% 
	pivot_longer(-sigma) %>% 
	mutate(name=case_when(
		name=="n_tested"~"Tested",
		name=="n_recruited"~"Recruited",
		name=="n_fast"~"Fast",
		name=="n_slow"~"Slow",
		name=="n_overall"~"Overall"
		))

fig_trial <- ggplot(trial_df_toplot, aes(x=sigma, y=value, col=name)) + 
	geom_point(size=0.5, alpha=0.6) + 
	geom_line(stat="smooth", method="loess", linewidth=1) + 
	theme_classic() 

# trial_df <- read_csv("output/trial_df_nejm.csv")





# To add: 
# custom age distribution to sim_trial
# custom fast target for sim_trials_over_sigma and sim_theory_over_sigma























# ==============================================================================
# Simulate trials under literature conditions
# ==============================================================================