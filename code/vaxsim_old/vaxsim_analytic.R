# Functions for analytic expressions
# A_i: denominator component for type i
A_i <- function(mu, a, rho) {
  rho * exp(-mu * a) * (exp((mu - rho) * a) - 1) / (mu - rho)
}

# N_i: numerator component for type i
N_i <- function(mu, p_i, a, sigma, rho) {
  p_i * rho * exp(-rho * a) * (exp((rho - mu) * sigma) - 1) / (rho - mu)
}

# P(infected in last sigma AND type i | asymptomatic at age a)
# P_infected_sigma_and_type <- function(mu, p_i, a, sigma, rho, mu_s, mu_f, p_s) {
#   N_i_val <- N_i(mu, p_i, a, sigma, rho)
#   D <- p_s * A_i(mu_s, a, rho) + (1 - p_s) * A_i(mu_f, a, rho)
#   N_i_val / D
# }

P_infected_sigma_and_type <- function(ptype, a, sigma, rho, mu_s, mu_f, p_s) {
  if(!ptype %in% c("fast","slow")){stop("Invalid ptype")}
  if(ptype=="fast"){mu <- mu_f} else {mu <- mu_s}
  if(ptype=="fast"){p_i <- 1-p_s} else {p_i <- p_s}
  N_i_val <- N_i(mu, p_i, a, sigma, rho)
  D <- p_s * A_i(mu_s, a, rho) + (1 - p_s) * A_i(mu_f, a, rho)
  N_i_val / D
}

# P(infected in last sigma | asymptomatic at age a)
P_infected_sigma <- function(a, sigma, rho, mu_s, mu_f, p_s) {
  p_f <- 1 - p_s
  N_s <- N_i(mu_s, p_s, a, sigma, rho)
  N_f <- N_i(mu_f, p_f, a, sigma, rho)
  D <- p_s * A_i(mu_s, a, rho) + p_f * A_i(mu_f, a, rho)
  print(c(N_s = N_s, N_f = N_f, D = D)) 
  (N_s + N_f) / D
}

P_infected_sigma_and_type(ptype="slow", a=30, sigma=2, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95)

P_infected_sigma(a=30, sigma=30, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95)