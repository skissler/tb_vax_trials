library(tidyverse) 

p_inf_and_type_and_asymp_given_age <- function(
	ptype, age, sigma, p_slow, incidence, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow 
	xi_f <- 1 - p_slow
	rho <- incidence 
	mu_s <- prograte_slow 
	mu_f <- prograte_fast 

	# Make sure sigma isn't greater than age: 	
	sigma <- min(sigma, age) 
  
  	# Check that the ptype is valid: 
	if (!(ptype %in% c("slow", "fast"))) stop("Invalid ptype")

	# Define the xi and mu for this type: 
	if (ptype == "slow") {
		xi <- xi_s
		mu <- mu_s
	} else {
		xi <- xi_f
		mu <- mu_f
	}

	# Calculate P(infected in [a, a-sigma], type, asymp at age a)
	if (abs(rho - mu) < 1e-8) {
		# Special case: rho == mu
		out <- xi * rho * sigma * exp(-rho * a)
	} else {
		out <- xi * rho / (rho - mu) * (
		  exp(-rho * (a - sigma)) * exp(-mu * sigma) - exp(-rho * a)
	)
	}

	return(out)

}

p_asymp_given_age <- function(
	age, p_slow, incidence, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow 
	xi_f <- 1 - p_slow
	rho <- incidence 
	mu_s <- prograte_slow 
	mu_f <- prograte_fast 

	# fast progressors
	if (abs(rho - mu_f) < 1e-8) {
		term_f <- xi_f * rho * a * exp(-rho * a)
	} else {
		term_f <- xi_f * rho / (mu_f - rho) * (exp(-rho * a) - exp(-mu_f * a))
	}

	# slow progressors
	if (abs(rho - mu_s) < 1e-8) {
		term_s <- xi_s * rho * a * exp(-rho * a)
	} else {
		term_s <- xi_s * rho / (mu_s - rho) * (exp(-rho * a) - exp(-mu_s * a))
	}

	# full probability
	out <- term_f + term_s + exp(-rho * a)

	return(out)

}

p_inf_and_type_given_asymp_and_age <- function(
	ptype, age, sigma, p_slow, incidence, prograte_slow, prograte_fast){
	  
	if (!(ptype %in% c("slow", "fast"))) stop("Invalid ptype")

	# Calculate the joint probability, i.e. 
	# P(infected in [a, a-sigma], type, asymp at age a)
	num <- p_inf_and_type_and_asymp_given_age(ptype=ptype, age=age, sigma=sigma, p_slow=p_slow, incidence=incidence, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	# Calculate the marginal probability of being asymptomatic at age a: 
	den <- p_asymp_given_age(age=age, p_slow=p_slow, incidence=incidence, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	out <- num / den
	return(out)
}

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

sim_trials_over_sigma <- function(pars, sigmavec, reps=25){

	trial_list <- vector("list", length(sigmavec)*reps)

	counter <- 1
	for(sigma in sigmavec){
		
		these_pars <- pars
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

	return(trial_df)
}

sim_theory_over_sigma <- function(pars, sigmavec, agedist){
	
	# Restrict to eligible age groups: 
	eligible <- agedist[names(agedist) %in% pars$minage:pars$maxage]
	eligible <- eligible/sum(eligible)

	theoretical_df <- vector("list", length(sigmavec))
	counter <- 1
	for(sigma in sigmavec){
		# Calculate a vector for drawing asymptomatic people of age a: 
		p_asymp_given_age_vec <- unlist(lapply(as.numeric(names(eligible)), function(x){
			p_asymp_given_age(age=x, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)
			}))
		p_age_given_asymp_vec <- p_asymp_given_age_vec * eligible / sum(p_asymp_given_age_vec * eligible)
		names(p_asymp_given_age_vec) <- names(eligible)

		# Calculate the probability vectors of testing positive and being fast/slow given asymptomatic and age a: 
		p_inf_and_slow_given_asymp_and_age_vec <- unlist(lapply(
			as.numeric(names(eligible)),
			function(x){p_inf_and_type_given_asymp_and_age(ptype="slow", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)}))
		names(p_inf_and_slow_given_asymp_and_age_vec) <- names(eligible)

		p_inf_and_fast_given_asymp_and_age_vec <- unlist(lapply(
			as.numeric(names(eligible)),
			function(x){p_inf_and_type_given_asymp_and_age(ptype="fast", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)}))
		names(p_inf_and_fast_given_asymp_and_age_vec) <- names(eligible)

		p_asymp <- sum(p_asymp_given_age_vec * eligible)

		p_fast <- sum(p_inf_and_fast_given_asymp_and_age_vec * p_age_given_asymp_vec)
		p_slow <- sum(p_inf_and_slow_given_asymp_and_age_vec * p_age_given_asymp_vec)
		p_pos  <- p_fast + p_slow

		tests_to_50_fast    <- 50 / p_fast
		recruits_to_50_fast <- 50 / (p_fast / p_pos)
		overall_to_50_fast <- 50 / (p_asymp * p_fast)

		theoretical_df[[counter]] <- list(
			sigma=sigma, 
			n_tested=tests_to_50_fast, 
			n_recruited=recruits_to_50_fast,
			n_fast=50,
			n_slow=recruits_to_50_fast-50,
			n_overall=overall_to_50_fast)
		counter <- counter + 1

	}

	theoretical_df <- bind_rows(theoretical_df)
	return(theoretical_df)

}

plot_trial_theory <- function(trial_df, theoretical_df, cols=c("n_tested","n_recruited","n_fast")){
	trial_df_toplot <- trial_df %>% 
		select(sigma, all_of(cols)) %>% 
		pivot_longer(-sigma) %>% 
		mutate(name=case_when(
			name=="n_tested"~"Screened",
			name=="n_recruited"~"Recruited",
			name=="n_fast"~"Fast",
			name=="n_slow"~"Slow",
			name=="n_overall"~"Overall"
			))

	theoretical_df_toplot <- theoretical_df %>% 
		select(sigma, all_of(cols)) %>% 
		pivot_longer(-sigma) %>% 
		mutate(name=case_when(
			name=="n_tested"~"Screened",
			name=="n_recruited"~"Recruited",
			name=="n_fast"~"Fast",
			name=="n_slow"~"Slow",
			name=="n_overall"~"Overall"
			))

	fig_trial_theory <- ggplot() + 
		geom_point(data=trial_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Overall","Screened","Recruited","Slow","Fast"))), size=0.5, alpha=0.2) + 
		geom_line(data=theoretical_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Overall","Screened","Recruited","Slow","Fast"))), linewidth=1, alpha=1) + 
		scale_color_manual(values=c("Overall"="green","Screened"="black","Recruited"="blue","Slow"="magenta","Fast"="red")) + 
		geom_vline(aes(xintercept=2), col="black", linetype="dashed", alpha=0.5) + 
		geom_vline(aes(xintercept=80), col="black", linetype="dashed", alpha=0.5) + 
		theme_classic() + 
		theme(legend.title=element_blank()) + 
		labs(x="Test span (years)", y="People")

	return(fig_trial_theory)
}

plot_screenslope <- function(theoretical_df){
	fig_screenslope <- theoretical_df %>% 
	select(sigma, n_tested) %>% 
	mutate(sigma_diff = sigma - lag(sigma)) %>% 
	mutate(tested_diff=n_tested - lag(n_tested)) %>% 
	mutate(slope=tested_diff/sigma_diff) %>% 
	filter(!is.na(slope)) %>% 
	ggplot(aes(x=sigma, y=slope)) + 
		geom_line(linewidth=1) + 
		theme_classic() + 
		labs(x="Test span (years)", y="Slope of screening line")
	return(fig_screenslope)
}