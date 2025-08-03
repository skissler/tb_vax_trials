# Imagine re-creating this trial: 
# https://www.nejm.org/doi/full/10.1056/NEJMoa1803484
# Ages 18-50
# Latent infection as confirmed by IGRA-positive requirement for entry 
# Outcome is bacteriologically confirmed active pulmonary TB disease
# Mean of 2.3 years of follow-up (aim is 3 years)
# 1786 participants per arm 
# Outcomes: 10 vaccinated and 22 placebo participants progressed 

library(tidyverse)
library(survival)

# UNCOMMENT FOR LOOPING MANY TRIALS ###################################
# n_fast_vec <- c()
# for(indexA in 1:100){
########################################################################

minage <- 18
maxage <- 50
rho <- 0.0275 # 0.0275 # .003  # 0.025
p_slow <- 0.5 # 0.95  # 0.75
mu_slow <- 0.0001
mu_fast <- 1.5 # 1.5 # 0.5 
sigma <- 2 # 100
ve <- 0.55
trial_length <- 3

# Initialize tracking variables  
n_tested <- 0
n_recruited <- 0
n_fast <- 0
n_overall <- 0

recruited_df <- tibble()

# for(grabs in 1:100000){
while(n_recruited < 3500){

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
# UNCOMMENT FOR LOOPING MANY TRIALS ###################################
# print(n_tested)
# n_fast_vec <- c(n_fast_vec, n_fast)
# }
########################################################################

rexpV <- Vectorize(rexp)
recruited_df <- recruited_df %>% 
	mutate(vaxstatus = sample(
		c(
		rep(0, ceiling(n_recruited/2)), 
		rep(1, ceiling(n_recruited/2)))[1:n_recruited]
		)) %>% 
	mutate(vaxworked = sample(c(TRUE,FALSE), 
		replace=TRUE, size=nrow(recruited_df), prob=c(ve,1-ve))) %>% 
	mutate(prograte = case_when(progressor_type=="slow"~mu_slow, TRUE~mu_fast)) %>% 
	mutate(prograte = case_when(vaxstatus==1 & vaxworked~mu_slow, TRUE~prograte)) %>% 
	# mutate(prograte = case_when(vaxstatus==1~prograte*(1-ve), TRUE~prograte)) %>% 
	mutate(tsymp_trial = rexpV(1, rate=prograte)) %>% 
	mutate(event=case_when(tsymp_trial < trial_length ~ 1, TRUE~0)) %>% 
	mutate(tsymp_trial=case_when(tsymp_trial < trial_length~tsymp_trial, TRUE~trial_length))


# Fit the Cox proportional hazards model
cox_model <- coxph(Surv(tsymp_trial, event) ~ vaxstatus, data = recruited_df)

# View the model summary
summary(cox_model)

nsummary <- recruited_df %>% filter(event==1) %>% group_by(vaxstatus) %>% summarise(n())
