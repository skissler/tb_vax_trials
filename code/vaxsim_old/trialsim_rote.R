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

# trial_df_nejm <- sim_trials_over_sigma(pars=pars_nejm_igra, sigmavec=1:80, reps=25)
trial_df_nejm <- read_csv("output/trial_df_nejm.csv")
theoretical_df_nejm <- sim_theory_over_sigma(pars=pars_nejm_igra, sigmavec=seq(from=.5, to=80, by=0.01), agedist=agedist)
fig_trial_theory_nejm <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, cols=c("n_tested","n_recruited","n_fast","n_overall")) + 
	labs(title="M72/AS01E trial emulation")
ggsave(fig_trial_theory_nejm, file="figures/trial_theory_nejm.pdf", width=5, height=5/1.6)

fig_trial_theory_log_nejm <- fig_trial_theory_nejm + scale_y_continuous(trans="log10")
ggsave(fig_trial_theory_log_nejm, file="figures/trial_theory_log_nejm.pdf", width=5, height=5/1.6)

fig_screenslope_nejm <- plot_screenslope(theoretical_df_nejm)

theoretical_df_nejm %>% 
	filter(sigma %in% c(0.5, 1, 1.5, 2, 5))



fig_trial_theory_nejm_r21 <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, cols=c("n_tested","n_recruited")) + 
	labs(title="M72/AS01E trial emulation")
ggsave(fig_trial_theory_nejm_r21, file="figures/trial_theory_nejm_r21.pdf", width=5, height=5/1.6)

fig_trial_theory_nejm_r21_log <- plot_trial_theory(trial_df_nejm, theoretical_df_nejm, cols=c("n_tested","n_recruited")) + 
	labs(title="M72/AS01E trial emulation") + 
	scale_y_continuous(trans="log10")
ggsave(fig_trial_theory_nejm_r21_log, file="figures/trial_theory_nejm_r21_log.pdf", width=5, height=5/1.6)

# To add: 
# custom age distribution to sim_trial
# custom fast target for sim_trials_over_sigma and sim_theory_over_sigma

# ==============================================================================
# Simulate trials under literature conditions
# ==============================================================================

# trial_df_lit <- sim_trials_over_sigma(pars=pars_lit_igra, sigmavec=1:80, reps=25)
trial_df_lit <- read_csv("output/trial_df_lit.csv")
theoretical_df_lit <- sim_theory_over_sigma(pars=pars_lit_igra, sigmavec=seq(from=.5, to=80, by=0.01), agedist=agedist)
fig_trial_theory_lit <- plot_trial_theory(trial_df_lit, theoretical_df_lit, cols=c("n_tested","n_recruited","n_fast","n_overall")) + 
	labs(title="Literature parameters")
ggsave(fig_trial_theory_lit, file="figures/trial_theory_lit.pdf", width=5, height=5/1.6)

fig_trial_theory_log_lit <- fig_trial_theory_lit + scale_y_continuous(trans="log10")
ggsave(fig_trial_theory_log_lit, file="figures/trial_theory_log_lit.pdf", width=5, height=5/1.6)

fig_screenslope_lit <- plot_screenslope(theoretical_df_lit)
