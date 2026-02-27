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
		out <- xi * rho * sigma * exp(-rho * a)  # from l'hopital's rule differentiating e.g. wrt mu
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

# Simulate functions
sim_analytic_over_sigma <- function(pars, sigmavec, agedist){
  # Restrict to eligible age groups
  eligible <- agedist[names(agedist) %in% pars$minage:pars$maxage]
  eligible <- eligible/sum(eligible)
  ages <- as.numeric(names(eligible))

  # Pre-compute sigma-independent quantities outside the for loop
  # Vector for drawing asymtomatic people of age a
  p_asymp_given_age_vec <- vapply(ages, function(x)
    p_asymp_given_age(age=x, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),  # probability for each age x
    numeric(1))
  names(p_asymp_given_age_vec) <- names(eligible)
  # Integrate over all a
  p_asymp <- sum(p_asymp_given_age_vec * eligible)
  # Defn of conditional probability
  p_age_given_asymp_vec <- p_asymp_given_age_vec * eligible / p_asymp

  analytical_df <- vector("list", length(sigmavec))
  counter <- 1
  for(sigma in sigmavec){
    # Calculate numerator P(infected in [a, a-sigma], type, asymp at age a)
    num_slow <- vapply(ages, function(x)
      p_inf_and_type_and_asymp_given_age(prog_type="slow", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),
      numeric(1))

    num_fast <- vapply(ages, function(x)
      p_inf_and_type_and_asymp_given_age(prog_type="fast", age=x, sigma=sigma, p_slow=pars$p_slow, incidence=pars$rho, prograte_slow=pars$mu_slow, prograte_fast=pars$mu_fast),
      numeric(1))

    # Divide by pre-computed denominator to get conditional probabilities
    p_inf_and_slow_given_asymp_and_age_vec <- num_slow / p_asymp_given_age_vec
    p_inf_and_fast_given_asymp_and_age_vec <- num_fast / p_asymp_given_age_vec

    p_fast     <- sum(p_inf_and_fast_given_asymp_and_age_vec * p_age_given_asymp_vec)
    p_slow_prop <- sum(p_inf_and_slow_given_asymp_and_age_vec * p_age_given_asymp_vec)
    p_pos      <- p_fast + p_slow_prop

    tests_to_50_fast    <- 50 / p_fast
    recruits_to_50_fast <- 50 / (p_fast / p_pos)
    overall_to_50_fast  <- 50 / (p_asymp * p_fast)

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

	# Pre-compute survival function pmf
	survival_fn <- exp(-cumsum(rho[1:100]))  # survival fn with case incidence by single-year ages
	survival_pmf <- c(1, survival_fn[-100]) - survival_fn  # prob infected during year a
	survival_pmf <- c(survival_pmf, survival_fn[100])  # adding tail for not infected during lifetime
	
	# Pre-compute some of the sampling vectors
	eligible_ages  <- minage:maxage
	eligible_probs <- agedist[minage:maxage]
	prog_types     <- c("slow","fast")
	prog_probs     <- c(p_slow, 1-p_slow)
	
	while(n_fast < fasttarget){
	  # Grab their age from the age distribution
	  age <- sample(eligible_ages, size=1, prob=eligible_probs)  # R normalises the subset automatically

		# Grab their progression status 
		progressor_type <- sample(prog_types, size=1, prob=prog_probs)

		# Simulate their time to infection, using age-specific case incidence vector rho
		tinf <- sample(0:100, size=1, prob=survival_pmf)
		if (tinf < 100) {
		  tinf <- tinf + runif(1)  # continuous within that year
		} else {
		  tinf <- rexp(1, 0.0025)  # not infected during lifetime; assume constant incidence for ages 101+?
		}
	
		# Simulate their time to infection, working backwards
		tsymp <- tinf + rexp(1, (if(progressor_type=="slow"){mu_slow} else {mu_fast}))  #EDIT HERE FOR VARYING RATE OF PROGRESSION BY AGE

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
