# Probability expression functions
p_inf_and_type_and_asymp_given_age <- function(  # Used to calculate no. needed to sample (page 13)
	prog_type, age, sigma, p_slow, foi, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow[as.character(age)]
	xi_f <- 1 - p_slow[as.character(age)]
	rho <- foi[as.character(age)]  # age-specific foi - only meaningful if uniform
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
		out <- xi * rho * sigma * exp(-rho * a)  # from l'hopital's rule differentiating e.g. wrt mu
	} else {
		out <- xi * rho / (rho - mu) * (
		  exp(-rho * (a - sigma)) * exp(-mu * sigma) - exp(-rho * a)
	)
	}

	return(out)
}

p_asymp_given_age <- function(   # Used in calculating no. needed to screen (page 14)
	age, p_slow, foi, prograte_slow, prograte_fast){

	# Rename variables 
	a <- age
	xi_s <- p_slow[as.character(age)]
	xi_f <- 1 - p_slow[as.character(age)]
	rho <- foi[as.character(age)]  # age-specific foi - only meaningful if uniform
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

p_inf_and_type_given_asymp_and_age <- function(  # Used in calculating no. needed to screen and to enroll (page 14/~17)  # Unused in code
	prog_type, age, sigma, p_slow, foi, prograte_slow, prograte_fast){
	  
	if (!(prog_type %in% c("slow", "fast"))) stop("Invalid prog_type")

	# Calculate the joint probability, i.e. 
	# P(infected in [a, a-sigma], type, asymp at age a)
	num <- p_inf_and_type_and_asymp_given_age(prog_type=prog_type, age=age, sigma=sigma, p_slow=p_slow, foi=foi, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	# Calculate the marginal probability of being asymptomatic at age a: 
	den <- p_asymp_given_age(age=age, p_slow=p_slow, foi=foi, prograte_slow=prograte_slow, prograte_fast=prograte_fast)

	out <- num / den
	
	return(out)
}

# Age structure functions
extract_age_distribution <- function(df, method="uniform"){
  if (method=="uniform") {
    agedist <- rep(1, 100)
    names(agedist) <- 0:99
  } else if (method %in% df$country) {
    agedist <- df %>% filter(country==method) %>% select(-country) %>% unlist()
  } else {
    print("method must be \"uniform\" or any UN-recognised country")
  }
  
  agedist <- agedist/sum(agedist)
  return(agedist)
}

extract_incidence_by_age <- function(df, method="uniform"){
  if (method=="uniform") {
    inc_by_age <- rep(0.0025, 100)
  } else if (method %in% df$country) {
    inc_by_age <- df %>% filter(country==method) %>% pull(inc)
  } else {
    print("method must be \"uniform\" or any UN-recognised country")
  }
  
  names(inc_by_age) <- 0:99
  return(inc_by_age)
}

get_pop <- function(x, country, unwpp) {
  if (country=="occupied Palestinian territory, including east Jerusalem") country <- "State of Palestine"
  unwpp_row <- unwpp[unwpp$country == country, ]
  ages <- as.character(pull_first(x):pull_last(x))
  as.numeric(rowSums(unwpp_row[, ages]))
}

impute_ARTI_by_age <- function(df) {
  df %>%
    group_by(study, country) %>%
    tidyr::complete(age = 0:99) %>%
    mutate(
      imputed = is.na(ARTI),
      ARTI = zoo::na.approx(ARTI, x = age, rule = 2)
    ) %>%
    ungroup()
}

extract_ARTI_by_age <- function(df, method="uniform", lo=F) {
  if (method=="South Africa" & lo==F) {  # specify which South Africa data to use (Wood or Ncayiyana)
    df <- df %>% filter(study=="Wood-2010") 
  } else if (method=="South Africa" & lo==T) {
    df <- df %>% filter(study=="Ncayiyana-2016")
  }
  
  if (method=="uniform") {
    ARTI_by_age <- rep(0.04, 100)
  } else if (method %in% unique(df$country)) {
    ARTI_by_age <- df %>% filter(country==method) %>% pull(ARTI)
  } else {
    print("method must be \"uniform\" or any country from `unique(ARTI_by_age_studies$country)`")
  }
  
  names(ARTI_by_age) <- 0:99
  return(ARTI_by_age)
}

impute_probability_of_being_slow <- function(df) {
  df %>%
    group_by(study) %>%
    tidyr::complete(age = 0:99) %>%
    mutate(
      p_slow = zoo::na.approx(p_slow, x = age, rule = 2),  # straight line interpolation
      p_fast = zoo::na.approx(p_fast, x = age, rule = 2),  # straight line interpolation
    ) %>%
    ungroup()
}

extract_probability_of_being_slow <- function(df, method="uniform"){
  if (method=="uniform") {
    p_slow_by_age <- rep(0.95, 100)
  } else if (method %in% df$study) {
    p_slow_by_age <- df %>% filter(study==method) %>% pull(p_slow)
  } else if (method=="Vynnycky-high") {
    p_fast_temp <- c(rep(0.04,11), 0.055, 0.07, 0.085, 0.10, 0.115, 0.13, 0.145, 0.16, 0.175, rep(0.19,80))  # 4% for ages 0-10, then linearly increasing, then 19% for age 20+
    p_slow_by_age <- 1 - p_fast_temp
  } else {
    print("method must be \"uniform\" or \"Vynnycky-high\" or any study from `unique(progressor_type_by_age_studies$study)`")
  }
  
  names(p_slow_by_age) <- 0:99
  return(p_slow_by_age)
}

make_foi <- function(a0_4, a5_24=NULL, a25_34=NULL, a35_49=NULL, a50plus=NULL) {
  # Helper: build age-varying FOI vector(s) using standard age bands (0-4, 5-24, 25-34, 35-49, 50+).
  # Values supplied as percentages (e.g. a0_4=2 means 2% per year for ages 0-4).
  # Single vector mode:  make_foi(2, 3, 4, 3, 2)
  # Grid mode:           make_foi(expand.grid(a0_4=1:5, a5_24=1:5, a25_34=1:5, a35_49=1:5, a50plus=1:5))
  #                      returns a named list of FOI vectors, names like "2_3_4_3_2"
  if (is.data.frame(a0_4)) {
    df <- a0_4
    vecs <- lapply(1:nrow(df), function(i) {
      foi <- c(rep(df$a0_4[i]/100, 5), rep(df$a5_24[i]/100, 20), rep(df$a25_34[i]/100, 10),
               rep(df$a35_49[i]/100, 15), rep(df$a50plus[i]/100, 50))
      setNames(foi, 0:99)
    })
    setNames(vecs, paste(df$a0_4, df$a5_24, df$a25_34, df$a35_49, df$a50plus, sep="_"))
  } else {
    foi <- c(rep(a0_4/100, 5), rep(a5_24/100, 20), rep(a25_34/100, 10), rep(a35_49/100, 15), rep(a50plus/100, 50))
    setNames(foi, 0:99)
  }
}

define_foi_by_age <- function(method){
  # raw vector pass-through: if method is already a length-100 numeric vector, use directly
  if (is.numeric(method) && length(method) == 100) {
    foi <- method
  # generic uniform_n: any string matching "uniform" followed by a number (e.g. "uniform_3", "uniform3.5")
  # value is interpreted as % per year (e.g. "uniform_3" -> 3% -> 0.03)
  } else if (grepl("^uniform_?[0-9]+(\\.[0-9]+)?$", method)) {
    n <- as.numeric(sub("^uniform_?", "", method))
    foi <- rep(n/100, 100)
  } else {
    stop(paste0("define_foi_by_age: unrecognised method '", method, "'"))
  }
  
  names(foi) <- 0:99
  return(foi)
}

define_mu_slow <- function(method) {
  if (method=="base") {
    mu_slow=0.001
  } else if (method=="min") {
    mu_slow=0.0001
  } else if (method=="max") {
    mu_slow=0.00527
  } else if (method=="0.003") {
    mu_slow=0.003
  } else if (method=="Vynnycky") {
    print("Functionality not yet added for Vynnycky (i.e. age-varying mu_slow)")
  }
}

# Utility functions
compute_survival_fn <- function(rho) {  # no reversion!
  # Compute survival function and PMF from age-specific FOI vector (length 100, names 0:99)
  # NB: survival_fn[a]  = P(not infected by age a)
  # NB: survival_pmf[a] = P(first infected during year a-1 to a)
  # NB: survival_pmf[101] = P(never infected during lifetime)
  survival_fn  <- exp(-cumsum(rho[1:100]))
  survival_pmf <- c(1, survival_fn[-100]) - survival_fn
  survival_pmf <- c(survival_pmf, survival_fn[100])
  list(survival_fn=survival_fn, survival_pmf=survival_pmf)
}

# Simulate functions
sim_analytic_over_sigma <- function(pars, sigmavec, casetarget=50, agedist){
  # Restrict to eligible age groups
  eligible <- agedist[names(agedist) %in% pars$minage:pars$maxage]
  eligible <- eligible/sum(eligible)
  ages <- as.numeric(names(eligible))

  # Pre-compute sigma-independent quantities outside the for loop
  # Vector for drawing asymptomatic people of age a
  p_asymp_given_age_vec <- vapply(ages, function(x)
    p_asymp_given_age(age=x, p_slow=pars$p_slow, foi=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),  # probability for each age x
    numeric(1))
  names(p_asymp_given_age_vec) <- names(eligible)
  # Integrate over all a
  p_asymp <- sum(p_asymp_given_age_vec * eligible)
  # Fraction of fast inf and slow inf individuals expected to progress to disease during trial
  fast_endpt_fraction <- 1 - exp(-pars$mu_fast * pars$trial_length)
  slow_endpt_fraction <- 1 - exp(-pars$mu_slow * pars$trial_length)
  print(fast_endpt_fraction)
  print(slow_endpt_fraction)

  analytical_df <- vector("list", length(sigmavec))
  counter <- 1
  for(sigma in sigmavec){
    # Calculate numerator P(infected in [a, a-sigma], type, asymp at age a)
    num_slow <- vapply(ages, function(x)
      p_inf_and_type_and_asymp_given_age(prog_type="slow", age=x, sigma=sigma, p_slow=pars$p_slow, foi=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),
      numeric(1))

    num_fast <- vapply(ages, function(x)
      p_inf_and_type_and_asymp_given_age(prog_type="fast", age=x, sigma=sigma, p_slow=pars$p_slow, foi=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),
      numeric(1))
    
    # Calculate probabilities needed for no. overall, no. tested, and no. enrolled
    p_inf_and_fast_and_asymp <- sum(num_fast * eligible)  # P1
    p_inf_and_slow_and_asymp <- sum(num_slow * eligible)  # P1'
    p_inf_and_fast_given_asymp <- sum(num_fast * eligible) / p_asymp  # P2
    p_inf_and_slow_given_asymp <- sum(num_slow * eligible) / p_asymp  # P2'
    
    overall_to_casetarget <- casetarget / (fast_endpt_fraction*p_inf_and_fast_and_asymp + slow_endpt_fraction*p_inf_and_slow_and_asymp)
    tests_to_casetarget    <- casetarget / (fast_endpt_fraction*p_inf_and_fast_given_asymp + slow_endpt_fraction*p_inf_and_slow_given_asymp)
    enrolls_to_casetarget <- casetarget * (p_inf_and_fast_given_asymp + p_inf_and_slow_given_asymp) / (fast_endpt_fraction*p_inf_and_fast_given_asymp + slow_endpt_fraction*p_inf_and_slow_given_asymp)
    
    analytical_df[[counter]] <- list(
      sigma=sigma,
      n_tested=tests_to_casetarget,
      n_enrolled=enrolls_to_casetarget,
      n_trial_cases=casetarget,
      n_overall=overall_to_casetarget)
    counter <- counter + 1
  }

  analytical_df <- bind_rows(analytical_df)
  return(analytical_df)
}

sim_stoch <- function(pars, casetarget=50, agedist, households=F){
	with(as.list(pars), {

	# Initialize tracking variables  
	n_overall <- 0
	n_tested <- 0
	n_enrolled <- 0
	n_fast <- 0
	n_trial_cases <- 0
	n_overall_15to24 <- 0  # only used for estimating inf prev
	n_overall_25to34 <- 0  # only used for estimating inf prev
	n_infected_15to24 <- 0  # only used for estimating inf prev
	n_infected_25to34 <- 0  # only used for estimating inf prev
	still_infected <- 1

	# Pre-compute survival function pmf - NB there is also a separate compute_survival_fn which does the same thing
	survival_fn <- exp(-cumsum(rho[1:100]))  # survival fn with foi by single-year ages
	survival_pmf <- c(1, survival_fn[-100]) - survival_fn  # prob infected during year a
	survival_pmf <- c(survival_pmf, survival_fn[100])  # adding tail for not infected during lifetime
	
	# Pre-compute some of the sampling vectors
	eligible_ages  <- minage:maxage
	eligible_probs <- agedist[as.character(eligible_ages)]
	prog_types     <- c("slow","fast")
	
	# # Pre-compute fraction of fast inf and slow inf individuals expected to progress to disease during trial - currently unused
	# fast_endpt_fraction <- (1 - exp(-(mu_fast + mu_revert2) * trial_length))*mu_fast/(mu_fast + mu_revert2)
	# slow_endpt_fraction <- (1 - exp(-(mu_slow + mu_revert2) * trial_length))*mu_slow/(mu_slow + mu_revert2)
	
	while(n_trial_cases < casetarget){
	  # Grab their age from the age distribution
	  age <- sample(eligible_ages, size=1, prob=eligible_probs)  # R normalises the subset automatically

	  # Simulate their time to infection, using survival_pmf probabilities
	  tinf <- sample(0:100, size=1, prob=survival_pmf)
	  if (tinf < 100) {
	    tinf <- tinf + runif(1)  # continuous within that year
	  } else {
	    tinf <- 100 + rexp(1, rho[100])  # not infected during lifetime; assume constant foi for ages 101+?
	  }
	  
	  # HOUSEHOLDS step
	  if (households) {
	    if (tinf > age) {  # not already infected
	      household_infection <- rbinom(1, size=1, prob=0.3)  # true or false coin flip
	      if (household_infection == TRUE) tinf <- age  # new tinf is current age
	    } else if (tinf <= age) {
	      # do nothing
	    }
	  }
	  
		# Grab their (tinf-dependent) progression status
	  if (tinf < 100) {
	    prog_probs <- c(p_slow[floor(tinf)+1], 1-p_slow[floor(tinf)+1])  # p_slow is a vector of age-varying probability of being slow
	  } else {
	    prog_probs <- c(p_slow[100], 1-p_slow[100])  # not infected during lifetime
	  }
	  progressor_type <- sample(prog_types, size=1, prob=prog_probs)
	  
	  # Grab their reversion status
	  reversion_status <- rbinom(1, size=1, prob=p_revert)  # true or false coin flip
	  
		# Simulate their time to symptoms and time to reversion
		tsymp <- tinf + rexp(1, (if(progressor_type=="slow") {mu_slow} else {mu_fast}))  # fixed rate of progression to disease
		trev <- (if(reversion_status==T) {tinf + rexp(1, mu_revert)} else {Inf})  # Inf = doesn't revert
		# But ensure only the first process of reversion/symptoms occurs
		tsymp <- if (tsymp > trev) Inf else tsymp
		trev <- if (trev > tsymp) Inf else trev
		still_infected <- (if(reversion_status==T & trev <= age) {F} else {T})
		
		# ELIGIBILITY
		if(tsymp > age){
			# They are asymptomatic, so let's test them: 
			n_tested <- n_tested + 1
			if((tinf >= age-sigma) && (tinf <= age)){
			  # They were infected in past sigma years but may have reverted. Check reversion.
				if(still_infected == T){
				  # They were infected in past sigma years and haven't reverted yet, so test positive. Enrol!
				  n_enrolled <- n_enrolled + 1
				  if(progressor_type=="fast") {n_fast <- n_fast + 1}
				  # exp_trial_cases <- fast_endpt_fraction*n_fast + slow_endpt_fraction*(n_enrolled - n_fast)  # expected no of cases during trial - exponential race method for interest
				  if(tsymp > age & tsymp <= age + trial_length) n_trial_cases <- n_trial_cases + 1  # record actual no of trial cases that would occur
				}
			}
		}

		# Also record those from n_overall who are positive, and by age (only used in ARTI/ inf prev calculations)
		if (age >= 15 & age < 25) {
		  n_overall_15to24 <- n_overall_15to24 + 1
		  if (tinf <= age & still_infected == T) n_infected_15to24 <- n_infected_15to24 + 1
		} else if (age >= 25 & age < 35) {
		  n_overall_25to34 <- n_overall_25to34 + 1
		  if (tinf <= age & still_infected == T) n_infected_25to34 <- n_infected_25to34 + 1
		}
		
		n_overall <- n_overall + 1
	}

	inf_prev_15to24 <- n_infected_15to24 / n_overall_15to24
	inf_prev_25to34 <- n_infected_25to34 / n_overall_25to34
	
	out <- list(n_tested=n_tested, n_enrolled=n_enrolled, n_fast=n_fast, n_slow=n_enrolled-n_fast, n_trial_cases=n_trial_cases, n_overall=n_overall, inf_prev_15to24=inf_prev_15to24, inf_prev_25to34=inf_prev_25to34)
	
	return(out)
	})
}
sim_stoch <- compiler::cmpfun(sim_stoch)  # R bytecode compiler (for speed)

sim_stoch_notest <- function(pars, casetarget=50, agedist, households=F){
  with(as.list(pars), {
    
    # Initialize tracking variables
    n_overall <- 0
    n_enrolled <- 0
    n_fast <- 0
    n_slow <- 0
    n_infected <- 0  # only used for calculating ARTI
    n_trial_cases <- 0
    n_overall_15to24 <- 0  # only used for estimating inf prev
    n_overall_25to34 <- 0  # only used for estimating inf prev
    n_infected_15to24 <- 0  # only used for estimating inf prev
    n_infected_25to34 <- 0  # only used for estimating inf prev
    still_infected <- 1
    
    # Pre-compute survival function pmf - NB there is also a separate compute_survival_fn which does the same thing
    survival_fn <- exp(-cumsum(rho[1:100]))  # survival fn with foi by single-year ages
    survival_pmf <- c(1, survival_fn[-100]) - survival_fn  # prob infected during year a
    survival_pmf <- c(survival_pmf, survival_fn[100])  # adding tail for not infected during lifetime
    
    # Pre-compute some of the sampling vectors
    eligible_ages  <- minage:maxage
    eligible_probs <- agedist[as.character(eligible_ages)]
    prog_types     <- c("slow","fast")
    
    while(n_trial_cases < casetarget){
      # Grab their age from the age distribution
      age <- sample(eligible_ages, size=1, prob=eligible_probs)  # R normalises the subset automatically
      
      # Simulate their time to infection, using survival_pmf probabilities
      tinf <- sample(0:100, size=1, prob=survival_pmf)
      if (tinf < 100) {
        tinf <- tinf + runif(1)  # continuous within that year
      } else {
        tinf <- 100 + rexp(1, rho[100])  # not infected during lifetime; assume constant foi for ages 101+?
      }
      
      # HOUSEHOLDS step
      if (households) {
        if (tinf > age) {  # not already infected
          household_infection <- rbinom(1, size=1, prob=0.3)  # true or false coin flip
          if (household_infection == TRUE) tinf <- age  # new tinf is current age
        } else if (tinf <= age) {
          # do nothing
        }
      }
      
      # Grab their (tinf-dependent) progression status 
      if (tinf < 100) {
        prog_probs <- c(p_slow[floor(tinf)+1], 1-p_slow[floor(tinf)+1])  # p_slow is a vector of age-varying probability of being slow
      } else {
        prog_probs <- c(p_slow[100], 1-p_slow[100])  # not infected during lifetime
      }
      progressor_type <- sample(prog_types, size=1, prob=prog_probs)
      
      # Grab their reversion status
      reversion_status <- rbinom(1, size=1, prob=p_revert)  # true or false coin flip
      
      # Simulate their time to symptoms and time to reversion
      tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))  # fixed rate of progression to disease
      trev <- (if(reversion_status==T) {tinf + rexp(1, mu_revert)} else {Inf})  # Inf = doesn't revert
      # But ensure only the first process of reversion/symptoms occurs
      tsymp <- if (tsymp > trev) Inf else tsymp
      trev <- if (trev > tsymp) Inf else trev
      still_infected <- (if(reversion_status==T & trev <= age) {F} else {T})
      
      # ELIGIBILITY
      if(tsymp > age){
        # They are asymptomatic, so let's enrol them: 
        n_enrolled <- n_enrolled + 1
        if (tinf <= age) {
          # They were infected in past sigma years but may have reverted. Check reversion.
          if (still_infected == T) {
            # They were infected in past sigma years and haven't reverted yet.
            n_infected <- n_infected + 1
            if(tsymp > age & tsymp <= age + trial_length) n_trial_cases <- n_trial_cases + 1  # record actual no of trial cases that would occur
          }
        }
      }

      n_overall <- n_overall + 1
      
      # Also record those from n_overall who would have been positive if tested, by age (only used in ARTI/ inf prev calculations)
      if (age >= 15 & age < 25) {
        n_overall_15to24 <- n_overall_15to24 + 1
        if (tinf <= age & still_infected == T) n_infected_15to24 <- n_infected_15to24 + 1
      } else if (age >= 25 & age < 35) {
        n_overall_25to34 <- n_overall_25to34 + 1
        if (tinf <= age & still_infected == T) n_infected_25to34 <- n_infected_25to34 + 1
      }
    }
    
    inf_prev_15to24 <- n_infected_15to24 / n_overall_15to24
    inf_prev_25to34 <- n_infected_25to34 / n_overall_25to34
    
    out <- list(n_enrolled=n_enrolled, n_trial_cases=n_trial_cases, n_overall=n_overall, n_infected=n_infected, inf_prev_15to24=inf_prev_15to24, inf_prev_25to34=inf_prev_25to34)
    
    return(out)
  })
}
sim_stoch_notest <- compiler::cmpfun(sim_stoch_notest)  # R bytecode compiler (for speed)

sim_stoch_trial <- function(pars, enrol_target, agedist) {
  # Like sim_stoch but stops when enrol_target enrolls are reached rather than casetarget cases.
  with(as.list(pars), {

    n_overall   <- 0
    n_tested    <- 0
    n_enrolled <- 0
    n_fast      <- 0
    n_trial_cases <- 0
    still_infected <- 1

    survival_fn  <- exp(-cumsum(rho[1:100]))
    survival_pmf <- c(1, survival_fn[-100]) - survival_fn
    survival_pmf <- c(survival_pmf, survival_fn[100])

    eligible_ages  <- minage:maxage
    eligible_probs <- agedist[as.character(eligible_ages)]
    prog_types     <- c("slow", "fast")

    while (n_enrolled < enrol_target) {
      # Grab their age from the age distribution
      age <- sample(eligible_ages, size=1, prob=eligible_probs)
      
      # Simulate their time to infection, using survival_pmf probabilities
      tinf <- sample(0:100, size=1, prob=survival_pmf)
      if (tinf < 100) {
        tinf <- tinf + runif(1)
      } else {
        tinf <- 100 + rexp(1, rho[100])
      }

      # Grab their (tinf-dependent) progression status
      prog_probs <- if (tinf < 100) c(p_slow[floor(tinf)+1], 1-p_slow[floor(tinf)+1]) else c(p_slow[100], 1-p_slow[100])
      progressor_type <- sample(prog_types, size=1, prob=prog_probs)
      
      # Grab their reversion status
      reversion_status <- rbinom(1, size=1, prob=p_revert)  # true or false coin flip
      
      # Simulate their time to symptoms and time to reversion
      tsymp <- tinf + rexp(1, if (progressor_type=="slow") mu_slow else mu_fast)
      trev <- (if(reversion_status==T) {tinf + rexp(1, mu_revert)} else {Inf})  # Inf = doesn't revert
      # But ensure only the first process of reversion/symptoms occurs
      tsymp <- if (tsymp > trev) Inf else tsymp
      trev <- if (trev > tsymp) Inf else trev
      still_infected <- (if(reversion_status==T & trev <= age) {F} else {T})
      
      # ELIGIBILITY
      if (tsymp > age) {
        # They are asymptomatic, so let's test them:
        n_tested <- n_tested + 1
        if ((tinf >= age - sigma) && (tinf <= age)) {
          # They were infected in past sigma years but may have reverted. Check reversion.
          if (still_infected == T) {
            # They were infected in past sigma years and haven't reverted yet, so test positive. Enrol!
            n_enrolled <- n_enrolled + 1
            if (progressor_type == "fast") n_fast <- n_fast + 1
            if(tsymp > age & tsymp <= age + trial_length) n_trial_cases <- n_trial_cases + 1  # record actual no of trial cases that would occur
          }
        }
      }
      n_overall <- n_overall + 1
    }

    out <- list(n_tested=n_tested, n_enrolled=n_enrolled, n_fast=n_fast, n_slow=n_enrolled-n_fast,
                n_trial_cases=n_trial_cases, n_overall=n_overall)
    return(out)
  })
}
sim_stoch_trial <- compiler::cmpfun(sim_stoch_trial)  # R bytecode compiler (for speed)

sim_stoch_over_sigma <- function(pars, sigmavec, casetarget=50, agedist, reps=25){
	grid <- expand.grid(sigma=sigmavec, rep=1:reps)  # expand.grid is faster than nested for loops
  p <- progressr::progressor(steps=nrow(grid))  # set up progress bar
  
	stoch_list <- future.apply::future_lapply(1:nrow(grid), function(i) {  # put everything into a lapply to use multiple cores
	  these_pars <- pars
	  these_pars$sigma <- grid$sigma[i]
	  stoch_output <- sim_stoch(these_pars, casetarget=casetarget, agedist=agedist)
	  p()  # report progress
	  list(sigma=grid$sigma[i], rep=grid$rep[i],
	    n_tested=stoch_output$n_tested, n_enrolled=stoch_output$n_enrolled,
	    n_fast=stoch_output$n_fast, exp_trial_cases=stoch_output$exp_trial_cases, 
	    n_trial_cases=stoch_output$n_trial_cases, n_overall=stoch_output$n_overall,
	    inf_prev_15to24=stoch_output$inf_prev_15to24, inf_prev_25to34=stoch_output$inf_prev_25to34)
	}, future.seed=T)  # future.seed does something important (each core does its own indep random number generation)
	
	bind_rows(stoch_list)
}

# Plot functions
plot_stochastic_analytic <- function(stochastic_df, analytical_df, cols=c("n_tested","n_enrolled","n_overall","n_trial_cases"), ylimit=F){
	stochastic_df_toplot <- stochastic_df %>% 
		select(sigma, all_of(cols)) %>% 
	  mutate(n_tested_1 = n_tested / n_trial_cases,  # scaling for no tested/enrolled/overall per 1 case
	         n_enrolled_1 = n_enrolled / n_trial_cases,
	         n_overall_1 = n_overall / n_trial_cases) %>%
	  select(sigma, n_tested_1, n_enrolled_1, n_overall_1) %>%
		pivot_longer(-sigma) %>%
		mutate(name=case_when(
			name=="n_tested_1"~"Tested",
			name=="n_enrolled_1"~"Enrolled",
			name=="n_overall_1"~"Contacted"
			))

	analytical_df_toplot <- analytical_df %>% 
		select(sigma, all_of(cols)) %>% 
	  mutate(n_tested_1 = n_tested / n_trial_cases,  # scaling for no tested/enrolled/overall per 1 case
	         n_enrolled_1 = n_enrolled / n_trial_cases,
	         n_overall_1 = n_overall / n_trial_cases) %>%
	  select(sigma, n_tested_1, n_enrolled_1, n_overall_1) %>%
		pivot_longer(-sigma) %>% 
		mutate(name=case_when(
			name=="n_tested_1"~"Tested",
			name=="n_enrolled_1"~"Enrolled",
			name=="n_overall_1"~"Contacted"
			))

	fig_stochastic_analytic <- ggplot() + 
		geom_point(data=stochastic_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Contacted","Tested","Enrolled","Cases"))), size=0.5, alpha=0.1) + 
		geom_line(data=analytical_df_toplot, aes(x=sigma, y=value, col=factor(name, levels=c("Contacted","Tested","Enrolled","Cases"))), linewidth=0.7, alpha=0.5) + 
		scale_color_manual(values=c("Contacted"="green","Tested"="black","Enrolled"="blue","Cases"="red")) + 
		geom_vline(aes(xintercept=2), col="black", linetype="dashed", alpha=0.5) + 
		geom_vline(aes(xintercept=80), col="black", linetype="dashed", alpha=0.5) + 
		theme_classic() + 
		theme(legend.title=element_blank()) + 
		labs(x="Test span (years)", y="Number (for 1 case)")

	if (ylimit) fig_stochastic_analytic <- fig_stochastic_analytic + ylim(0,ylimit)
	return(fig_stochastic_analytic)
}

# Character string manipulation functions
numericize <- function(df, col){  # from Measles code
  df %>% rowwise() %>%
    mutate({{col}} := list(pull_first({{col}}):pull_last({{col}}))) %>%
    unnest_longer({{col}})
}

pull_first <- function(x){  # for character string x
  if (x=="all") return(0)
  as.numeric(stringr::str_extract(x, "^[[:digit:]\\.]+"))
}

pull_last <- function(x, life_exp=99){  # for character string x
  if (x == "all") return(life_exp)
  if (stringr::str_detect(x, "(\\+|plus)$")) return(life_exp)
  as.numeric(stringr::str_extract(x, "[[:digit:]\\.]+$"))
}

pull_mean <- function(x, life_exp=99){
  (pull_first(x) + pull_last(x, life_exp)) / 2
}

# Functions outside the model, at time of recruitment
calculate_case_incidence_at_time_of_recruitment_from_model <- function(df, pars, my_rep=1, sig=80){  # case inc = # of new cases / (pop * trial_length)
  if ("rep" %in% names(df)) df <- df %>% filter(rep == my_rep)
  if ("sigma" %in% names(df)) df <- df %>% filter(sigma == sig)
  model_run <- df %>% unlist()
  # case incidence (new cases)
  incidence <- unname(model_run["n_trial_cases"] / (model_run["n_overall"]*pars$trial_length))
  return(incidence*100000)
}

estimate_ARTI_from_model <- function(df, pars, my_rep=1, sig=80){
  if ("rep" %in% names(df)) df <- df %>% filter(rep == my_rep)
  if ("sigma" %in% names(df)) df <- df %>% filter(sigma == sig)
  model_run <- df %>% unlist()
  # infection prevalence all ages
  if ("n_tested" %in% names(df)) {
    infection_prevalence <- unname((model_run["n_overall"] - model_run["n_tested"] + model_run["n_enrolled"]) / model_run["n_overall"])  # total infection prevalence at start of trial (for all trial ages)
  } else {
    infection_prevalence <- unname((model_run["n_overall"] - model_run["n_enrolled"] + model_run["n_infected"]) / model_run["n_overall"])  # total infection prevalence at start of trial (for all trial ages)
  }
  meanage <- ((pars$maxage - pars$minage)/2) + pars$minage
  ARTI <- 1 - ((1 - infection_prevalence)^(1/meanage))
  # infection prevalence for key age groups
  inf_prev_15to24 <- unname(model_run["inf_prev_15to24"])
  inf_prev_25to34 <- unname(model_run["inf_prev_25to34"])
  return(c(ARTI = ARTI, inf_prev_15to24 = inf_prev_15to24, inf_prev_25to34 = inf_prev_25to34))
}