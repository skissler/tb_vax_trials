# ==============================================================================
# Import 
# ==============================================================================

library(tidyverse) 
source('code/utils.R')

# ==============================================================================
# Initial simulation 
# ==============================================================================

# Let's look at test windows of 2 to 30 years: ---------------------------------
sigmavals <- seq(from=2, to=30, by=0.1)

# Some literature-derived parameters: ------------------------------------------
# pfastvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp_and_age(ptype="fast", age=30, sigma=x, p_slow=0.95, incidence=0.0025, prograte_slow=0.0001, prograte_fast=1.5)}))
# pslowvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp_and_age(ptype="slow", age=30, sigma=x, p_slow=0.95, incidence=0.0025, prograte_slow=0.0001, prograte_fast=1.5)}))

# Some parameters to match the NEJM trial: -------------------------------------
pfastvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp_and_age(ptype="fast", age=30, sigma=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
pslowvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp_and_age(ptype="slow", age=30, sigma=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))


# Some quick counts, to estimate tests needed to get 40 fast progressors: ------
ntestvec <- 40/pfastvec
nfastvec <- ntestvec*pfastvec
nslowvec <- ntestvec*pslowvec
nposvec <- nfastvec + nslowvec
nnegvec <- ntestvec - nposvec

# Aggregate into a data frame for plotting: 
trial_df <- tibble(
	sigma=sigmavals,
	ntests=ntestvec,
	nfast=nfastvec,
	nslow=nslowvec,
	npos=nposvec,
	nneg=nnegvec)

fig_trial <- trial_df %>% 
	select(sigma, `Tested` = ntests, `Recruited (positive test)` = npos, `Fast progressors` = nfast, `Slow progressors` = nslow) %>% 
	pivot_longer(-sigma) %>% 
	mutate(name=factor(name, levels=c("Tested", "Recruited (positive test)", "Slow progressors", "Fast progressors"))) %>% 
	ggplot(aes(x = sigma, y = value, col = name, linetype = name)) + 
		geom_line(linewidth = 0.75, alpha=0.8) + 
		scale_linetype_manual(values = c(
			"Tested" = "solid", 
			"Recruited (positive test)" = "solid", 
			"Fast progressors" = "dashed", 
			"Slow progressors" = "dashed"
		)) + 
		scale_color_manual(values=c(
			"Tested" = "black", 
			"Recruited (positive test)" = "purple", 
			"Fast progressors" = "red", 
			"Slow progressors" = "blue"
		)) + 
		labs(
			x="Duration of test positivity\n(years post-infection)",
			y="Number of people") + 
		theme_classic() + 
		theme(legend.title=element_blank())

# ==============================================================================
# Simulate with uncertainty and age distribution
# ==============================================================================

# Define an age distribution (uniform for now): 
agedist <- rep(1, 80)
names(agedist) <- 1:80
agedist <- agedist/sum(agedist)

# Restrict to eligible age groups: 
eligible <- agedist[names(agedist) %in% 18:49]
eligible <- eligible/sum(eligible)

# Sample a person's age: 
# age <- as.numeric(sample(names(eligible), size=1, prob=eligible))

# Calculate the probabillity they're asymptomatic: 
# p_asymp_given_age(age=age, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)

# Calculate a vector for drawing asymptomatic people of age a: 
p_asymp_given_age_vec <- unlist(lapply(as.numeric(names(eligible)), function(x){
	p_asymp_given_age(age=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)
	}))
p_age_given_asymp_vec <- p_asymp_given_age_vec * eligible / sum(p_asymp_given_age_vec * eligible)
names(p_asymp_given_age_vec) <- names(eligible)

# Calculate the probability vectors of testing positive and being fast/slow given asymptomatic and age a: 
p_inf_and_slow_given_asymp_and_age_vec <- unlist(lapply(
	as.numeric(names(eligible)),
	function(x){p_inf_and_type_given_asymp_and_age(ptype="slow", age=x, sigma=age_i, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
names(p_inf_and_slow_given_asymp_and_age_vec) <- names(eligible)

p_inf_and_fast_given_asymp_and_age_vec <- unlist(lapply(
	as.numeric(names(eligible)),
	function(x){p_inf_and_type_given_asymp_and_age(ptype="fast", age=x, sigma=age_i, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
names(p_inf_and_fast_given_asymp_and_age_vec) <- names(eligible)



trial_df <- tibble()
n_recruited <- 0
n_fast <- 0
n_slow <- 0
n_tested <- 0
while(n_fast < 50){
	# LOOPING NOW HERE
	# If a person is asymptomatic, we can draw their age using this vector. 
	age_i <- as.numeric(sample(names(eligible), size=1, prob=p_age_given_asymp_vec))

	# We can quickly calculate the outcome of their test: 
	draw_i <- runif(1) 

	pfast_i <- p_inf_and_fast_given_asymp_and_age_vec[as.character(age_i)]
	pslow_i <- p_inf_and_slow_given_asymp_and_age_vec[as.character(age_i)]

	if(draw_i < pfast_i) {
		type_i <- "fast"
		n_fast <- n_fast + 1

		recruited_i <- TRUE
		n_recruited <- n_recruited + 1
	} else if(draw_i < (pfast_i + pslow_i)) {
		type_i <- "slow"
		n_slow <- n_slow + 1

		recruited_i <- TRUE
		n_recruited <- n_recruited + 1
	} else {
		type_i <- NA_character_
		recruited_i <- FALSE
	}

	trial_df <- bind_rows(trial_df, tibble(
		age_i=age_i,
		type_i=type_i,
		recruited_i=recruited_i
		))

	n_tested <- n_tested + 1
}




# Calculate the probability they test positive and are fast/slow, given asymp
# pslow_i <- p_inf_and_type_given_asymp_and_age(ptype="slow", age=age_i, sigma=age_i, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)

# pfast_i <- p_inf_and_type_given_asymp_and_age(ptype="fast", age=age_i, sigma=age_i, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)











# Scratch ======================================================================
# agevec <- 1:100
# p_asymp_given_age_vec <- unlist(lapply(agevec, function(x){p_asymp_given_age(age=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
# fig_p_asymp <- tibble(age=agevec, p_asymp=p_asymp_given_age_vec) |> 
# 	ggplot(aes(x=age, y=p_asymp)) + 
# 		geom_line() + 
# 		scale_y_continuous(limits=c(0,1)) + 
# 		theme_classic() 












































