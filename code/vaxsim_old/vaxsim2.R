library(tidyverse)

# Parameters
n <- 10000              # total individuals
a <- 30                 # age at recruitment
sigma <- 1              # infection window (years)
rho <- 0.001            # infection hazard (incidence: infections/person/year)
mu_s <- 0.0001          # slow progression rate (1/years)
mu_f <- 1.5             # fast progression rate (1/years)
p_s <- 0.95              # proportion of slow progressors
p_f <- 1 - p_s

# Vaccine effect (reduces progression rate)
rho_v <- 0.5            # vaccine reduces expected progression time by 50% (speeds rate by 1/0.5 = 2)
mu_s_v <- mu_s / rho_v
mu_f_v <- mu_f / rho_v

# Trial setup
arm <- sample(c("vaccine", "placebo"), n, replace = TRUE)

# Precompute components for posterior infection category probabilities
A_i <- function(mu) {
  exp(-rho * a) + rho * exp(-mu * a) * (exp((mu - rho) * a) - 1) / (mu - rho)
}

N_i <- function(mu, p_i) {
  infect_window <- exp(-rho * (a - sigma)) - exp(-rho * a)
  rho * exp(-mu * a) * (exp((mu - rho) * a) - exp((mu - rho) * (a - sigma))) / (mu - rho) * p_i * infect_window
}

A_s <- A_i(mu_s)
A_f <- A_i(mu_f)
D <- p_s * A_s + p_f * A_f
P_s <- N_i(mu_s, p_s) / D
P_f <- N_i(mu_f, p_f) / D

# Sample latent class
progressor_type <- sample(c("slow", "fast"), n, replace = TRUE, prob = c(P_s, P_f))

# Function to sample infection time from truncated tilted exponential
draw_infection_time <- function(mu) {
  delta <- mu - rho
  if (abs(delta) < 1e-6) {
    runif(1, a - sigma, a)
  } else {
    u <- runif(1)
    log( u * (exp(delta * a) - exp(delta * (a - sigma))) + exp(delta * (a - sigma)) ) / delta
  }
}

# Simulate infection time and progression time
infection_time <- numeric(n)
time_to_progression <- numeric(n)

for (i in 1:n) {
  type <- progressor_type[i]
  arm_i <- arm[i]
  mu <- if (type == "slow") mu_s else mu_f
  if (arm_i == "vaccine") mu <- mu / rho_v
  
  s <- draw_infection_time(if (type == "slow") mu_s else mu_f)
  infection_time[i] <- s
  time_to_progression[i] <- rexp(1, rate = mu)  # time from recruitment to symptoms
}

# Censoring at 5 years
censor_time <- 5
observed <- time_to_progression <= censor_time
follow_up_time <- pmin(time_to_progression, censor_time)

# Create data frame
trial_data <- as_tibble(data.frame(
  arm = arm,
  progressor_type = progressor_type,
  infection_time = infection_time,
  time_to_progression = time_to_progression,
  follow_up_time = follow_up_time,
  observed = observed
))

# Summarize
table(trial_data$arm, trial_data$observed)

# Kaplan-Meier curve
library(survival)
library(survminer)

surv_obj <- Surv(trial_data$follow_up_time, trial_data$observed)
fit <- survfit(surv_obj ~ trial_data$arm)

fig_km <- ggsurvplot(
  fit,
  data = trial_data,
  conf.int = TRUE,
  pval = TRUE,
  risk.table = TRUE,
  legend.labs = c("Placebo", "Vaccine"),
  xlab = "Time (years)",
  ylab = "Survival probability",
  ggtheme = theme_minimal()
)


# Bayesian comparison of progression rates
# Use Gamma(a0, b0) prior for exponential rate (uninformative: a0 = b0 = 0.001)
a0 <- 0.001
b0 <- 0.001

# Compute posterior parameters for each group
library(dplyr)
posterior_params <- trial_data %>%
  group_by(arm) %>%
  summarise(
    sum_time = sum(follow_up_time),
    num_events = sum(observed),
    a_post = a0 + num_events,
    b_post = b0 + sum_time
  )

# Extract posterior samples (gamma draws)
n_draws <- 10000
post_vaccine <- rgamma(n_draws, shape = posterior_params$a_post[posterior_params$arm == "vaccine"],
                       rate = posterior_params$b_post[posterior_params$arm == "vaccine"])
post_placebo <- rgamma(n_draws, shape = posterior_params$a_post[posterior_params$arm == "placebo"],
                       rate = posterior_params$b_post[posterior_params$arm == "placebo"])

# Posterior distribution of rate ratio
rate_ratio <- post_vaccine / post_placebo

# Probability that vaccine reduces progression rate
mean(rate_ratio < 1)  # posterior probability vaccine is effective

# Optional: plot posterior distribution of rate ratio
library(ggplot2)
fig_postratio <- ggplot(data.frame(rate_ratio = rate_ratio), aes(x = rate_ratio)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(title = "Posterior Distribution of Rate Ratio",
       x = "Vaccine / Placebo Progression Rate",
       y = "Density") +
  theme_minimal()


trial_data %>% 
  filter(observed==TRUE) %>% 
  group_by(arm, progressor_type) %>% 
  summarise(meantime=mean(follow_up_time), n=n()) 

trial_data %>% 
  mutate(yval=case_when(arm=="placebo"~0, TRUE~1)) %>% 
  filter(observed==TRUE) %>% 
  ggplot(aes(x=follow_up_time, y=yval)) + 
    geom_jitter(height=0.1) 


