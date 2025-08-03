library(tidyverse) 

p_inf_and_type_given_asymp <- function(
	ptype, age, sigma, p_slow, incidence, prograte_slow, prograte_fast){
	
	# Rename variables 
	a <- age
	xi_s <- p_slow 
	xi_f <- 1 - p_slow
	rho <- incidence 
	mu_s <- prograte_slow 
	mu_f <- prograte_fast 
  
	if (!(ptype %in% c("slow", "fast"))) stop("Invalid ptype")

	if (ptype == "slow") {
		xi <- xi_s
		mu <- mu_s
	} else {
		xi <- xi_f
		mu <- mu_f
	}

	# Numerator
	if (abs(rho - mu) < 1e-8) {
		# Special case: rho == mu
		num <- xi * rho * sigma * exp(-rho * a)
	} else {
		num <- xi * rho / (rho - mu) * (
		  exp(-rho * (a - sigma)) * exp(-mu * sigma) - exp(-rho * a)
	)
	}

	# Denominator
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

	den <- term_f + term_s + exp(-rho * a)

	out <- num / den
	return(out)
}


sigmavals <- seq(from=2, to=30, by=0.1)
# pfastvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp(ptype="fast", age=30, sigma=x, p_slow=0.95, incidence=0.0025, prograte_slow=0.0001, prograte_fast=1.5)}))
# pslowvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp(ptype="slow", age=30, sigma=x, p_slow=0.95, incidence=0.0025, prograte_slow=0.0001, prograte_fast=1.5)}))

pfastvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp(ptype="fast", age=30, sigma=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))
pslowvec <- unlist(lapply(sigmavals, function(x){p_inf_and_type_given_asymp(ptype="slow", age=30, sigma=x, p_slow=0.55, incidence=0.0275, prograte_slow=0.0001, prograte_fast=1.5)}))

ntestvec <- 40/pfastvec
nfastvec <- ntestvec*pfastvec
nslowvec <- ntestvec*pslowvec
nposvec <- nfastvec + nslowvec
nnegvec <- ntestvec - nposvec

test_df <- tibble(
	sigma=sigmavals,
	ntests=ntestvec,
	nfast=nfastvec,
	nslow=nslowvec,
	npos=nposvec,
	nneg=nnegvec)

test_df %>% 
	ggplot(aes(x=sigma, y=ntests)) + 
		geom_line(linewidth=1) + 
		theme_classic()

test_df %>% 
	select(sigma, ntests, npos) %>% 
	pivot_longer(-sigma) %>% 
	ggplot(aes(x=sigma, y=value, col=name)) + 
		geom_line(linewidth=1) + 
		theme_classic()

test_df %>% 
	select(sigma, `Tested`=ntests, `Recruited (positive test)`=npos, `Fast progressors`=nfast, `Slow progressors`=nslow) %>% 
	pivot_longer(-sigma) %>% 
	mutate(dashed=case_when(name %in% c("Fast progressors","Slow progressors") ~ "Y", TRUE~"N")) %>% 
	ggplot(aes(x=sigma, y=value, col=name, lty=dashed)) + 
		geom_line(linewidth=1) + 
		theme_classic()

fig_comp <- test_df %>% 
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

