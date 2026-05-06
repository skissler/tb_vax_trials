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
  } else {
    print("method must be \"uniform\" or any study from `unique(progressor_type_by_age_studies$study)`")
  }
  
  names(p_slow_by_age) <- 0:99
  return(p_slow_by_age)
}

define_foi_by_age <- function(method){
  if (method=="uniform1") {
    foi <- rep(0.01,100)
  } else if (method=="uniform2") {
    foi <- rep(0.02,100)
  } else if (method=="uniform4") {
    foi <- rep(0.04,100)
  } else if (method=="uniform6") {
    foi <- rep(0.06,100)
  } else if (method=="uniform10") {
    foi <- rep(0.1,100)
  } else if (method=="5888") {
    foi <- c(rep(0.05,20), rep(0.08,20), rep(0.08,40), rep(0.08,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="5102030") {
    foi <- c(rep(0.05,20), rep(0.10,20), rep(0.20,40), rep(0.30,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="5204050") {
    foi <- c(rep(0.05,20), rep(0.20,20), rep(0.40,40), rep(0.50,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="5255070") {
    foi <- c(rep(0.05,20), rep(0.25,20), rep(0.50,40), rep(0.70,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="33510") {
    foi <- c(rep(0.03,20), rep(0.03,20), rep(0.05,40), rep(0.10,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="151015") {
    foi <- c(rep(0.01,20), rep(0.05,20), rep(0.10,40), rep(0.15,20))  # age 0-19, 20-39, 40-79, 80+
  } else if (method=="base") {
    foi <- c(rep(0.07,25), rep(0.07,5), rep(0.07,20), rep(0.07,50))  # age 0-24, 25-29, 30-49, 50+
  } else if (method=="min") {
    foi <- c(rep(0.05,25), rep(0.01,5), rep(0.01,20), rep(0.01,50))  # age 0-24, 25-29, 30-49, 50+
  } else if (method=="max") {
    foi <- c(rep(0.10,25), rep(0.10,5), rep(0.10,20), rep(0.10,50))  # age 0-24, 25-29, 30-49, 50+
  } else if (method=="incr") {
    foi <- c(rep(0.07,25), rep(0.08,5), rep(0.09,20), rep(0.10,50))  # age 0-24, 25-29, 30-49, 50+
  } else if (method=="decr") {
    foi <- c(rep(0.07,25), rep(0.05,5), rep(0.03,20), rep(0.01,50))  # age 0-24, 25-29, 30-49, 50+
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
  } else if (method=="Vynnycky") {
    print("Functionality not yet added for Vynnycky (i.e. age-varying mu_slow)")
  }
}

# Utility functions
compute_survival_fn <- function(rho) {
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
sim_analytic_over_sigma <- function(pars, sigmavec, fasttarget=50, agedist){
  # Restrict to eligible age groups
  eligible <- agedist[names(agedist) %in% pars$minage:pars$maxage]
  eligible <- eligible/sum(eligible)
  ages <- as.numeric(names(eligible))

  # Pre-compute sigma-independent quantities outside the for loop
  # Vector for drawing asymtomatic people of age a
  p_asymp_given_age_vec <- vapply(ages, function(x)
    p_asymp_given_age(age=x, p_slow=pars$p_slow, foi=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),  # probability for each age x
    numeric(1))
  names(p_asymp_given_age_vec) <- names(eligible)
  # Integrate over all a
  p_asymp <- sum(p_asymp_given_age_vec * eligible)

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
    
    # Calculate probabilities needed for no. overall, no. tested, and no. recruited
    p_inf_and_fast_and_asymp <- sum(num_fast * eligible)
    p_inf_and_fast_given_asymp <- sum(num_fast * eligible) / p_asymp  # used to be called p_fast but it is actually P_{pos&fast|asymp} from p15
    p_inf_and_slow_given_asymp <- sum(num_slow * eligible) / p_asymp  # used to be called p_slow or p_slow_prop
    p_fast_given_inf_and_asymp <- p_inf_and_fast_given_asymp / (p_inf_and_fast_given_asymp + p_inf_and_slow_given_asymp)
    
    recruits_to_fasttarget <- fasttarget / p_fast_given_inf_and_asymp
    tests_to_fasttarget    <- fasttarget / p_inf_and_fast_given_asymp
    overall_to_fasttarget <- fasttarget / p_inf_and_fast_and_asymp
    
    analytical_df[[counter]] <- list(
      sigma=sigma,
      n_tested=tests_to_fasttarget,
      n_recruited=recruits_to_fasttarget,
      n_fast=fasttarget,
      n_slow=recruits_to_fasttarget - fasttarget,
      n_overall=overall_to_fasttarget)
    counter <- counter + 1
  }

  analytical_df <- bind_rows(analytical_df)
  return(analytical_df)
}

sim_stoch <- function(pars, fasttarget=50, agedist, households=F){
	with(as.list(pars), {

	# Initialize tracking variables  
	n_tested <- 0
	n_recruited <- 0
	n_fast <- 0
	n_overall <- 0

	capacity <- 1e6
	recruited_list <- vector("list", capacity)  # specifying list size in advance for speed

	# Pre-compute survival function pmf - NB there is also a separate compute_survival_fn which does the same thing
	survival_fn <- exp(-cumsum(rho[1:100]))  # survival fn with foi by single-year ages
	survival_pmf <- c(1, survival_fn[-100]) - survival_fn  # prob infected during year a
	survival_pmf <- c(survival_pmf, survival_fn[100])  # adding tail for not infected during lifetime
	
	# Pre-compute some of the sampling vectors
	eligible_ages  <- minage:maxage
	eligible_probs <- agedist[as.character(eligible_ages)]
	prog_types     <- c("slow","fast")
	
	while(n_fast < fasttarget){
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
	  
		# Simulate their time to symptoms
		tsymp <- tinf + rexp(1, (if(progressor_type=="slow") {mu_slow} else {mu_fast}))  # fixed rate of progression to disease
		
		# ELIGIBILITY
		if(tsymp > age){
			# They are asymptomatic, so let's test them: 
			n_tested <- n_tested + 1
			if((tinf >= age-sigma) && (tinf <= age)){
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
					tinf=tinf  # only recording tinf values for those who were recruited
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
sim_stoch <- compiler::cmpfun(sim_stoch)  # R bytecode compiler (for speed)

sim_stoch_notest <- function(pars, fasttarget=50, agedist, households=F){
  with(as.list(pars), {
    
    # Initialize tracking variables
    n_recruited <- 0
    n_fast <- 0
    n_slow <- 0
    n_overall <- 0
    n_infected <- 0  # only used for calculating ARTI
    
    capacity <- 1e6
    
    # Pre-compute survival function pmf - NB there is also a separate compute_survival_fn which does the same thing
    survival_fn <- exp(-cumsum(rho[1:100]))  # survival fn with foi by single-year ages
    survival_pmf <- c(1, survival_fn[-100]) - survival_fn  # prob infected during year a
    survival_pmf <- c(survival_pmf, survival_fn[100])  # adding tail for not infected during lifetime
    
    # Pre-compute some of the sampling vectors
    eligible_ages  <- minage:maxage
    eligible_probs <- agedist[as.character(eligible_ages)]
    prog_types     <- c("slow","fast")
    
    while(n_fast < fasttarget){
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
      
      # Simulate their time to symptoms
      tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))  # fixed rate of progression to disease
      
      # ELIGIBILITY
      if(tsymp > age){
        # They are asymptomatic, so let's recruit them: 
        n_recruited <- n_recruited + 1
        # If they are fast (and actually infected), add 1 to n_fast:
        if (tinf <= age) {
          (if (progressor_type=="fast") {n_fast <- n_fast + 1} else {n_slow <- n_slow + 1})
        }
      }

      n_overall <- n_overall + 1
      
      # Also record those from n_overall who would have been positive if tested (only used in ARTI calculations)
      if (tinf <= age) n_infected <- n_infected + 1
    }
    
    out <- list(n_recruited=n_recruited, n_fast=n_fast, n_slow=n_slow, n_overall=n_overall, n_infected=n_infected)
    
    return(out)
  })
}
sim_stoch_notest <- compiler::cmpfun(sim_stoch_notest)  # R bytecode compiler (for speed)

sim_trial <- function(pars, enrol_target, agedist) {
  # Like sim_stoch but stops when enrol_target recruits are reached rather than fasttarget fast progressors.
  with(as.list(pars), {

    n_tested    <- 0
    n_recruited <- 0
    n_fast      <- 0
    n_overall   <- 0

    capacity       <- enrol_target * 20
    recruited_list <- vector("list", capacity)

    survival_fn  <- exp(-cumsum(rho[1:100]))
    survival_pmf <- c(1, survival_fn[-100]) - survival_fn
    survival_pmf <- c(survival_pmf, survival_fn[100])

    eligible_ages  <- minage:maxage
    eligible_probs <- agedist[as.character(eligible_ages)]
    prog_types     <- c("slow", "fast")

    while (n_recruited < enrol_target) {
      age  <- sample(eligible_ages, size=1, prob=eligible_probs)
      tinf <- sample(0:100, size=1, prob=survival_pmf)
      if (tinf < 100) {
        tinf <- tinf + runif(1)
      } else {
        tinf <- 100 + rexp(1, rho[100])
      }

      prog_probs      <- if (tinf < 100) c(p_slow[floor(tinf)+1], 1-p_slow[floor(tinf)+1]) else c(p_slow[100], 1-p_slow[100])
      progressor_type <- sample(prog_types, size=1, prob=prog_probs)
      tsymp           <- tinf + rexp(1, if (progressor_type=="slow") mu_slow else mu_fast)

      if (tsymp > age) {
        n_tested <- n_tested + 1
        if ((tinf >= age - sigma) && (tinf <= age)) {
          n_recruited <- n_recruited + 1
          if (progressor_type == "fast") n_fast <- n_fast + 1
          if (n_recruited > capacity) { capacity <- capacity * 2; length(recruited_list) <- capacity }
          recruited_list[[n_recruited]] <- list(
            id=n_recruited, age=age, progressor_type=progressor_type, tinf=tinf
          )
        }
      }
      n_overall <- n_overall + 1
    }

    recruited_df <- bind_rows(recruited_list[1:n_recruited])

    out <- list(n_tested=n_tested, n_recruited=n_recruited, n_fast=n_fast, n_slow=n_recruited-n_fast,
                n_overall=n_overall, recruited_df=recruited_df)
    return(out)
  })
}
sim_trial <- compiler::cmpfun(sim_trial)  # R bytecode compiler (for speed)

sim_stoch_over_sigma <- function(pars, sigmavec, fasttarget=50, agedist, reps=25){
	grid <- expand.grid(sigma=sigmavec, rep=1:reps)  # expand.grid is faster than nested for loops
  p <- progressr::progressor(steps=nrow(grid))  # set up progress bar
  
	stoch_list <- future.apply::future_lapply(1:nrow(grid), function(i) {  # put everything into a lapply to use multiple cores
	  these_pars <- pars
	  these_pars$sigma <- grid$sigma[i]
	  stoch_output <- sim_stoch(these_pars, fasttarget=fasttarget, agedist=agedist)
	  p()  # report progress
	  list(sigma=grid$sigma[i], rep=grid$rep[i],
	    n_tested=stoch_output$n_tested, n_recruited=stoch_output$n_recruited,
	    n_fast=stoch_output$n_fast, n_slow=stoch_output$n_slow, n_overall=stoch_output$n_overall)
	}, future.seed=T)  # future.seed does something important (each core does its own indep random number generation)
	
	bind_rows(stoch_list)
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

# Functions outside the model
estimate_case_incidence_from_model <- function(df, pars, my_rep=1, sig=80){  # case inc = # of new cases / (pop * trial_length)
  if ("rep" %in% names(df)) df <- df %>% filter(rep == my_rep)
  if ("sigma" %in% names(df)) df <- df %>% filter(sigma == sig)
  model_run <- df %>% unlist()
  incidence <- unname((model_run["n_fast"] + model_run["n_slow"]*pars$mu_slow*pars$trial_length) / (model_run["n_overall"]*pars$trial_length))
  return(incidence*100000)
}

estimate_inf_prev_from_model <- function(foi, ages = c(20, 30, 40, 60)) {
  # P(ever infected by age a) = 1 - S(a), derived directly from the FOI
  sv <- compute_survival_fn(foi)
  inf_prev <- 1 - sv$survival_fn  # index i = age i (survival_fn[1] = P(not infected by age 1))
  setNames(inf_prev[ages], paste0("inf_prev_age", ages))
}

estimate_ARTI_from_model <- function(df, pars, my_rep=1, sig=80){
  if ("rep" %in% names(df)) df <- df %>% filter(rep == my_rep)
  if ("sigma" %in% names(df)) df <- df %>% filter(sigma == sig)
  model_run <- df %>% unlist()
  if ("n_tested" %in% names(df)) {
    infection_prevalence <- unname((model_run["n_overall"] - model_run["n_tested"] + model_run["n_recruited"]) / model_run["n_overall"])  # total infection prevalence at start of trial (for eligible ages)
  } else {
    infection_prevalence <- unname(model_run["n_infected"] / model_run["n_overall"])  # total infection prevalence at start of trial (for eligible ages)
  }
  meanage <- ((pars$maxage - pars$minage)/2) + pars$minage
  ARTI <- 1 - ((1 - infection_prevalence)^(1/meanage))
  return(c(ARTI = ARTI))
  }