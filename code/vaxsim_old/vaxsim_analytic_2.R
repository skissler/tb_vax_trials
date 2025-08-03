
# Numerator: P(infected in last sigma AND asymptomatic at age a), for given type
N_i <- function(mu, p_i, a, sigma, rho) {
  p_i * integrate(function(s) {
    rho * exp(-rho * s) * exp(-mu * (a - s))
  }, lower = a - sigma, upper = a)$value
}

# Denominator: P(asymptomatic at age a)
D_total <- function(a, rho, mu_s, mu_f, p_s) {
  p_f <- 1 - p_s
  p_uninfected <- exp(-rho * a)
  p_s * integrate(function(s) {
    rho * exp(-rho * s) * exp(-mu_s * (a - s))
  }, lower = 0, upper = a)$value +
    p_f * integrate(function(s) {
      rho * exp(-rho * s) * exp(-mu_f * (a - s))
    }, lower = 0, upper = a)$value + p_uninfected
}

# P(infected in last sigma AND type i | asymptomatic at age a)
P_infected_sigma_and_type <- function(ptype, a, sigma, rho, mu_s, mu_f, p_s) {
  if(!ptype %in% c("fast","slow")){stop("Invalid ptype")}
  if(ptype=="fast"){mu <- mu_f} else {mu <- mu_s}
  if(ptype=="fast"){p_i <- 1-p_s} else {p_i <- p_s}
  N_i_val <- N_i(mu, p_i, a, sigma, rho)
  D <- D_total(a, rho, mu_s, mu_f, p_s)
  N_i_val / D
}

# P(infected in last sigma | asymptomatic at age a)
P_infected_sigma <- function(a, sigma, rho, mu_s, mu_f, p_s) {
  p_f <- 1 - p_s
  N_s <- N_i(mu_s, p_s, a, sigma, rho)
  N_f <- N_i(mu_f, p_f, a, sigma, rho)
  D <- D_total(a, rho, mu_s, mu_f, p_s)
  (N_s + N_f) / D
}

P_infected_sigma_and_type(ptype="fast", a=30, sigma=2, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95)

P_infected_sigma(a=30, sigma=2, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95)

P_infected_sigma_and_type(ptype="fast", a=30, sigma=30, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95) / P_infected_sigma(a=30, sigma=30, rho=0.001, mu_s=0.0001, mu_f=1.5, p_s=0.95)

# Tests needed to get 50 fast progressors: 
50 / P_infected_sigma_and_type(ptype="fast", a=30, sigma=30, rho=0.025, mu_s=0.0001, mu_f=1.5, p_s=0.5)