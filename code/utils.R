library(tidyverse) 

# Probability expression functions
p_inf_and_type_and_asymp_given_age <- function(  # Used to calculate no. needed to sample (page 13)
	prog_type, age, sigma, p_slow, incidence, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow 
	xi_f <- 1 - p_slow
	rho <- incidence[age]  # age-specific incidence
	mu_s <- prograte_slow 
	mu_f <- prograte_fast 

	# Make sure sigma isn't greater than age: 	
	sigma <- min(sigma, age) 
  
  # Check that the prog_type is valid: 
	if (!(prog_type %in% c("slow", "fast"))) stop("Invalid prog_type")

	# Define the xi and mu for this type: 
	if (prog_type == "slow") {
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

p_asymp_given_age <- function(   # Used in calculating no. needed to screen (page 14)
	age, p_slow, incidence, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow 
	xi_f <- 1 - p_slow
	rho <- incidence[age]  # age-specific incidence
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

	# full probability (fasts + slows + never infecteds)
	out <- term_f + term_s + exp(-rho * a)

	return(out)
}

p_inf_and_type_given_asymp_and_age <- function(  # Used in calculating no. needed to screen and to enroll (page 14/~17)
	prog_type, age, sigma, p_slow, incidence, prograte_slow, prograte_fast){
	  
	if (!(prog_type %in% c("slow", "fast"))) stop("Invalid prog_type")

	# Calculate the joint probability, i.e. 
	# P(infected in [a, a-sigma], type, asymp at age a)
	num <- p_inf_and_type_and_asymp_given_age(prog_type=prog_type, age=age, sigma=sigma, p_slow=p_slow, incidence=incidence, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	# Calculate the marginal probability of being asymptomatic at age a: 
	den <- p_asymp_given_age(age=age, p_slow=p_slow, incidence=incidence, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	out <- num / den
	
	return(out)
}

# Age distribution functions
extract_age_distribution <- function(method="uniform", life_exp=80){
  if (method=="uniform") {
    agedist <- rep(1, life_exp)
  } else if (method=="test") {
    agedist <- c(1:life_exp)
  }
  
  names(agedist) <- 1:length(agedist)
  agedist <- agedist/sum(agedist)
  return(agedist)
}

extract_incidence_by_age <- function(method="uniform", life_exp=80){
  if (method=="uniform") {
    inc_by_age <- rep(1, life_exp)
  } else if (method=="test") {
    global_TB_report <- read.csv(file = "data/TB_burden_age_sex_2026-02-20.csv")  ## NEED TO FILTER
    inc_by_age <- c(1:life_exp)
  }
  
  names(inc_by_age) <- 1:life_exp
  inc_by_age <- inc_by_age/sum(inc_by_age)
  return(inc_by_age)
}

# Simulate functions
sim_analytic_over_sigma <- function(pars, sigmavec, agedist){
  # Restrict to eligible age groups: 
  eligible <- agedist[names(agedist) %in% pars$minage:pars$maxage]
  eligible <- eligible/sum(eligible)
  
  analytical_df <- vector("list", length(sigmavec))
  counter <- 1
  for(sigma in sigmavec){
    # Calculate a vector for drawing asymptomatic people of age a: 
    p_asymp_given_age_vec <- unlist(lapply(as.numeric(names(eligible)), function(x){  # probability for each age x
      p_asymp_given_age(age=x, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)
    }))
    p_age_given_asymp_vec <- p_asymp_given_age_vec * eligible / sum(p_asymp_given_age_vec * eligible)
    names(p_asymp_given_age_vec) <- names(eligible)
    
    # Calculate the probability vectors of testing positive and being fast/slow given asymptomatic and age a: 
    p_inf_and_slow_given_asymp_and_age_vec <- unlist(lapply(
      as.numeric(names(eligible)),
      function(x){p_inf_and_type_given_asymp_and_age(prog_type="slow", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)}))
    names(p_inf_and_slow_given_asymp_and_age_vec) <- names(eligible)
    
    p_inf_and_fast_given_asymp_and_age_vec <- unlist(lapply(
      as.numeric(names(eligible)),
      function(x){p_inf_and_type_given_asymp_and_age(prog_type="fast", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast)}))
    names(p_inf_and_fast_given_asymp_and_age_vec) <- names(eligible)
    
    p_asymp <- sum(p_asymp_given_age_vec * eligible)
    
    p_fast <- sum(p_inf_and_fast_given_asymp_and_age_vec * p_age_given_asymp_vec)
    p_slow <- sum(p_inf_and_slow_given_asymp_and_age_vec * p_age_given_asymp_vec)
    p_pos  <- p_fast + p_slow
    
    tests_to_50_fast    <- 50 / p_fast
    recruits_to_50_fast <- 50 / (p_fast / p_pos)
    overall_to_50_fast <- 50 / (p_asymp * p_fast)
    
    analytical_df[[counter]] <- list(
      sigma=sigma, 
      n_tested=tests_to_50_fast, 
      n_recruited=recruits_to_50_fast,
      n_fast=50,
      n_slow=recruits_to_50_fast-50,
      n_overall=overall_to_50_fast)
    counter <- counter + 1
  }
  
  analytical_df <- bind_rows(analytical_df)
  return(analytical_df)
}

sim_stoch <- function(pars, fasttarget=50, agedist){
	with(as.list(pars), {

	# Initialize tracking variables  
	n_tested <- 0
	n_recruited <- 0
	n_fast <- 0
	n_overall <- 0

	capacity <- 1e6
	recruited_list <- vector("list", capacity)  # specifying list size in advance for speed

	while(n_fast < fasttarget){
	  # Grab their age from the age distribution
	  age <- sample(minage:maxage, size=1, prob=agedist[minage:maxage])  # R normalises the subset automatically

		# Grab their progression status 
		progressor_type <- sample(c("slow","fast"), 
			size=1, prob=c(p_slow, 1-p_slow))

		# Simulate their time of infection, using age-specific incidence rho[age]
		tinf <- rexp(1, rho[age])

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

	recruited_df <- bind_rows(recruited_list[1:(n_recruited)])  # bind_rows outside of the while loop is much faster

	out <- list(n_tested=n_tested, n_recruited=n_recruited, n_fast=n_fast, n_slow=n_recruited-n_fast, n_overall=n_overall, recruited_df=recruited_df)
	
	return(out)
	})
}

sim_stoch_over_sigma <- function(pars, sigmavec, fasttarget=50, agedist, reps=25){
	stoch_list <- vector("list", length(sigmavec)*reps)

	counter <- 1
	for(sigma in sigmavec){
		
		these_pars <- pars
		these_pars$sigma <- sigma 
		
		for(rep in 1:reps){
			
			stoch_output <- sim_stoch(these_pars, fasttarget=50, agedist)

			stoch_list[[counter]] <- list(
				sigma=sigma,
				rep=rep,
				n_tested=stoch_output$n_tested,
				n_recruited=stoch_output$n_recruited,
				n_fast=stoch_output$n_fast,
				n_slow=stoch_output$n_slow,
				n_overall=stoch_output$n_overall)
			counter <- counter + 1
		}
		print(sigma)
	}

	stochastic_df <- bind_rows(stoch_list)

	return(stochastic_df)
}

# Plot functions
plot_stochastic_analytic <- function(stochastic_df, analytical_df, cols=c("n_tested","n_recruited","n_fast")){
	stochastic_df_toplot <- stochastic_df %>% 
		select(sigma, all_of(cols)) %>% 
		pivot_longer(-sigma) %>% 
		mutate(name=case_when(
			name=="n_tested"~"Screened",
			name=="n_recruited"~"Recruited",
			name=="n_fast"~"Fast",
			name=="n_slow"~"Slow",
			name=="n_overall"~"Overall"
			))

	analytical_df_toplot <- analytical_df %>% 
		select(sigma, all_of(cols)) %>% 
		pivot_longer(-sigma) %>% 
		mutate(name=case_when(
			name=="n_tested"~"Screened",
			name=="n_recruited"~"Recruited",
			name=="n_fast"~"Fast",
			name=="n_slow"~"Slow",
			name=="n_overall"~"Overall"
			))

	fig_stochastic_analytic <- ggplot() + 
		geom_point(data=stochastic_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Overall","Screened","Recruited","Slow","Fast"))), size=0.5, alpha=0.2) + 
		geom_line(data=analytical_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Overall","Screened","Recruited","Slow","Fast"))), linewidth=1, alpha=1) + 
		scale_color_manual(values=c("Overall"="green","Screened"="black","Recruited"="blue","Slow"="magenta","Fast"="red")) + 
		geom_vline(aes(xintercept=2), col="black", linetype="dashed", alpha=0.5) + 
		geom_vline(aes(xintercept=80), col="black", linetype="dashed", alpha=0.5) + 
		theme_classic() + 
		theme(legend.title=element_blank()) + 
		labs(x="Test span (years)", y="People")

	return(fig_stochastic_analytic)
}

plot_screenslope <- function(analytical_df){
	fig_screenslope <- analytical_df %>% 
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