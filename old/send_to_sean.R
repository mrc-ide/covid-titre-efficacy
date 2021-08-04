library(tidyverse)
library(patchwork)

# code shown to implement dose 1 titre and efficacy - as soon as the second dose is given, would replace titre with the dose 2 titre and repeat the same waning/efficacy algorithm

# parameters for user to input (example values shown here)
mu_ab_d1 <- 0.14 # mean titre dose 1
mu_ab_d2 <- 2.37 # mean titre dose 2
ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
ab_50_severe <- 0.03
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
k <- 2.94 # shape parameter of efficacy curve
max_t <- 730 # number of days to model
t_d2 <- 84 # timing of second dose relative to first
hl_s <- 108 # Half life of antibody decay - short
hl_l <- 3650 # Half life of antibody decay - long
period_s <- 250
t_period_l <- 365 # Time point at which to have switched to longest half-life
dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

nt <- NULL
t <- 0:max_t
time_to_decay <- t_period_l - period_s # time in days to reach longest half-life

# vector of decay rates over time: first we have a period of fast decay, then gradually shift to period of long decay
dr_vec <- c(rep(dr_s, period_s),
                seq(dr_s, dr_l, length.out = time_to_decay),
                rep(dr_l, (length(t) - t_period_l)))

# nab titre distribution for a given vaccine product: draw from a log-normal distribution
z1 <- rnorm(1, log10(mu_ab_d1), std10)

# initiate titre vector
nt <- rep(0, length(t))
nt[1] <- log(10^z1)

# decay antibodies over time on natural log scale
for (i in (2:length(t))){
  nt[i] <- nt[i-1] + dr_vec[i]
}
nt <- exp(nt) # return to linear scale
  
# relate titre to efficacy over time - using log-10 parameters
ef_infection <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50))))
ef_severe <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_severe))))

