library(tidyverse) 

# Some quick visualizations: 

a <- 30 
mu_s <- 0.0001 
mu_f <- 1.5 
sigma <- 2
rho <- .001

pfunc <- function(a, mu, sigma, rho){

	num <- 1/(a*mu)*(1 - exp(-a*mu))*(1 - exp(-rho*a))
	den <- num + exp(-rho*a)
	p <- num/den
	if(sigma < a){
		adj <- (1 - exp(-mu*sigma)) / (1 - exp(-mu*a))
	} else {
		adj <- 1
	}

	return(p*adj)

}

pfunc(a=a, mu=mu_f, sigma=30, rho=rho)


sigmavals <- 1:50
temp <- unlist(lapply(sigmavals, function(x){pfunc(a=a, mu=mu_s, sigma=x, rho=100)}))

plot(sigmavals, temp)










# First question: how many people do we need to sample to get N participants? 











# Randomly pull 2N people from a population with uniform age distribution: 
N <- 100 
ages <- 80*runif(2*N)

