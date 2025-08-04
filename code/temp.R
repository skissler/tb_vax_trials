library(tidyverse) 

source('code/utils.R')

pars_nejm_igra <- list(
	minage=18,
	maxage=50,
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
	maxage=50,
	rho=0.0025,
	p_slow=0.95,
	mu_slow=0.0001,
	mu_fast=1.5,
	sigma=100,
	ve=0.55,
	trial_length=3)
pars_lit_tasa <- pars_lit_igra 
pars_lit_tasa$sigma <- 2

sim_trial <- function(pars, fasttarget=50){
	with(as.list(pars), {

	# Initialize tracking variables  
	n_tested <- 0
	n_recruited <- 0
	n_fast <- 0
	n_overall <- 0
	recruited_df <- tibble()

	while(n_fast < fasttarget){

		# Grab their age 
		age <- runif(1, min=minage, max=maxage)

		# Grab their progression status 
		progressor_type <- sample(c("slow","fast"), 
			size=1, prob=c(p_slow, 1-p_slow))

		# Simulate their time of infection
		tinf <- rexp(1, rho)

		# Simulate their time of symptoms 
		tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))

		# ELIGIBILITY
		if(tsymp > age){
			# They are asymptomatic, so let's test them: 
			n_tested <- n_tested + 1
			if((tinf >= age-sigma) & (tinf <= age)){
				# They were infected in past sigma years, so test positive. Recruit!
				n_recruited <- n_recruited + 1
				if(progressor_type=="fast"){n_fast <- n_fast + 1}
				recruited_df <- bind_rows(recruited_df, tibble(
					id=n_recruited,
					age=age,
					progressor_type=progressor_type,
					tinf=tinf
					))
			}
		}

		n_overall <- n_overall + 1

	}

	out <- list(n_tested=n_tested, recruited_df=recruited_df)
	return(out)

	})
}



sim_trial_fast <- function(pars, fasttarget=50){
	with(as.list(pars), {

	# Initialize tracking variables  
	n_tested <- 0
	n_recruited <- 0
	n_fast <- 0
	n_overall <- 0

	capacity <- 1e6
	recruited_list <- vector("list", capacity)

	while(n_fast < fasttarget){

		# Grab their age 
		age <- runif(1, min=minage, max=maxage)

		# Grab their progression status 
		progressor_type <- sample(c("slow","fast"), 
			size=1, prob=c(p_slow, 1-p_slow))

		# Simulate their time of infection
		tinf <- rexp(1, rho)

		# Simulate their time of symptoms 
		tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))

		# ELIGIBILITY
		if(tsymp > age){
			# They are asymptomatic, so let's test them: 
			n_tested <- n_tested + 1
			if((tinf >= age-sigma) & (tinf <= age)){
				# They were infected in past sigma years, so test positive. Recruit!
				n_recruited <- n_recruited + 1
				if(progressor_type=="fast"){n_fast <- n_fast + 1}

				# Expand list if needed
				if (n_recruited > capacity) {
					capacity <- capacity * 2
					length(recruited_list) <- capacity  
				}

				recruited_list[[n_recruited]] <- list(
					id=n_recruited,
					age=age,
					progressor_type=progressor_type,
					tinf=tinf
					)

			}
		}

		n_overall <- n_overall + 1

	}

	recruited_df <- bind_rows(recruited_list[1:(n_recruited)])

	out <- list(n_tested=n_tested, n_recruited=n_recruited, n_fast=n_fast, n_slow=n_recruited-n_fast, n_overall=n_overall, recruited_df=recruited_df)
	return(out)

	})
}

# ==============================================================================
# Some testing 
# ==============================================================================


sigmavec <- 1:80
reps <- 25

trial_list <- vector("list", length(sigmavec)*reps)

counter <- 1
for(sigma in sigmavec){
	
	these_pars <- pars_nejm_igra
	these_pars$sigma <- sigma 
	
	for(rep in 1:reps){
		
		trial_output <- sim_trial_fast(these_pars)

		trial_list[[counter]] <- list(
			sigma=sigma,
			rep=rep,
			n_tested=trial_output$n_tested,
			n_recruited=trial_output$n_recruited,
			n_fast=trial_output$n_fast,
			n_slow=trial_output$n_slow,
			n_overall=trial_output$n_overall)
		counter <- counter + 1
	}
	print(sigma)
}

trial_df <- bind_rows(trial_list)

trial_means_df <- trial_df %>% 
	group_by(sigma) %>% 
	summarise(
		n_tested=mean(n_tested), 
		n_recruited=mean(n_recruited),
		n_fast=mean(n_fast),
		n_slow=mean(n_slow),
		n_overall=mean(n_overall))

trial_df_toplot <- trial_df %>% 
	select(sigma, n_tested, n_recruited, n_fast) %>% 
	pivot_longer(-sigma)

trial_means_df_toplot <- trial_means_df %>% 	
	select(sigma, n_tested, n_recruited, n_fast) %>% 
	pivot_longer(-sigma)

fig_trial <- ggplot(trial_df_toplot, aes(x=sigma, y=value, col=name)) + 
	geom_point(size=0.5, alpha=0.6) + 
	geom_line(stat="smooth", method="loess", linewidth=1) + 
	theme_classic() 

fig_trial_means <- ggplot() + 
	geom_point(data=trial_df_toplot, aes(x=sigma, y=value, col=name), size=0.5, alpha=0.6) + 
	# geom_point(data=trial_means_df_toplot, aes(x=sigma, y=value, col=name), size=2, alpha=0.6) + 
	geom_line(data=trial_means_df_toplot, aes(x=sigma, y=value, col=name), alpha=0.6, linewidth=2) + 
	theme_classic() 

# ==============================================================================
# Try including a theoretical line: 

# Define an age distribution (uniform for now): 
agedist <- rep(1, 80)
names(agedist) <- 1:80
agedist <- agedist/sum(agedist)

# Restrict to eligible age groups: 
eligible <- agedist[names(agedist) %in% 18:49]
eligible <- eligible/sum(eligible)

theoretical_df <- vector("list", length(sigmavec))
counter <- 1
for(sigma in sigmavec){
	# Calculate a vector for drawing asymptomatic people of age a: 
	p_asymp_given_age_vec <- unlist(lapply(as.numeric(names(eligible)), function(x){
		p_asymp_given_age(age=x, p_slow=0.5, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)
		}))
	p_age_given_asymp_vec <- p_asymp_given_age_vec * eligible / sum(p_asymp_given_age_vec * eligible)
	names(p_asymp_given_age_vec) <- names(eligible)

	# Calculate the probability vectors of testing positive and being fast/slow given asymptomatic and age a: 
	p_inf_and_slow_given_asymp_and_age_vec <- unlist(lapply(
		as.numeric(names(eligible)),
		function(x){p_inf_and_type_given_asymp_and_age(ptype="slow", age=x, sigma=sigma, p_slow=0.5, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
	names(p_inf_and_slow_given_asymp_and_age_vec) <- names(eligible)

	p_inf_and_fast_given_asymp_and_age_vec <- unlist(lapply(
		as.numeric(names(eligible)),
		function(x){p_inf_and_type_given_asymp_and_age(ptype="fast", age=x, sigma=sigma, p_slow=0.5, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
	names(p_inf_and_fast_given_asymp_and_age_vec) <- names(eligible)

	p_fast <- sum(p_inf_and_fast_given_asymp_and_age_vec * p_age_given_asymp_vec)
	p_slow <- sum(p_inf_and_slow_given_asymp_and_age_vec * p_age_given_asymp_vec)
	p_pos  <- p_fast + p_slow

	tests_to_50_fast    <- 50 / p_fast
	recruits_to_50_fast <- 50 / (p_fast / p_pos)

	theoretical_df[[counter]] <- list(
		sigma=sigma, 
		n_tested=tests_to_50_fast, 
		n_recruited=recruits_to_50_fast,
		n_fast=50)
	counter <- counter + 1

}

theoretical_df <- bind_rows(theoretical_df)

theoretical_df_toplot <- theoretical_df %>% 
	pivot_longer(-sigma)

fig_trial_theory <- ggplot() + 
	geom_point(data=trial_df_toplot, aes(x=sigma, y=value, col=name), size=0.5, alpha=0.6) + 
	geom_line(data=theoretical_df_toplot, aes(x=sigma, y=value, col=name), linewidth=1, alpha=0.6) + 
	theme_classic() 

