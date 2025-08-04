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