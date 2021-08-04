library(tidyverse)
library(patchwork)

# parameters
mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac", "Moderna"),
                         mu_ab_d1 = c(13/94, 1/59, 28/164, ((185+273)/2)/321),
                         mu_ab_d2 = c(223/94, 32/59, 28/164, 654/158))

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

nt <- NULL
t <- 0:max_t
time_to_decay <- t_period_l - period_s # time in days to reach longest half-life
dr_s <- -log(2)/hl_s # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

# function to create titre and efficacy profiles
draw <- function(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, std10 = std10, ab_50 = ab_50, ab_50_severe, dr_s = dr_s, dr_l = dr_l, t = t, t_d2 = t_d2, k = k, time_to_decay = time_to_decay){
  
  # decay rates over time
  dr_vec_sub <- c(rep(dr_s, period_s),
                  seq(dr_s, dr_l, length.out = time_to_decay),
                  rep(dr_l, (length(t) - t_period_l)))
  
  dr_vec <- rep(0, length(t))
  dr_vec[1:(t_d2-1)] <- dr_vec_sub[1:(t_d2-1)]
  dr_vec[t_d2:length(t)] <- dr_vec_sub[1:(length(t) - t_d2 + 1) ]

  # nab titre distribution for a given vaccine product: draw from a log-normal distribution
  # here the first titre is dependent on the second but I suppose could do it the other way around
  z2 <- rnorm(1, log10(mu_ab_d2), std10)
  z1 <- log10(10^z2 * (mu_ab_d1/mu_ab_d2))

  # initiate titre vector
  nt <- rep(0, length(t))
  nt[1] <- log(10^z1)
  nt[t_d2] <- log(10^z2)
  
  # decay antibodies over time on natural log scale
  for (i in (2:(t_d2-1))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  for (i in ((t_d2+1):length(t))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  nt <- exp(nt) # return to linear scale
  
  # relate titre to efficacy over time - using log-10 parameters
  ef_infection <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50))))
  ef_severe <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_severe))))

  # output list containing titre and vaccine efficacy
  ret <- list(titre = nt, VE = ef_infection, VE_severe = ef_severe)
  return(ret)
}

plots <- NULL
plotlist <- list()
for (j in 1:4){
  mu_ab_d1 <- mu_ab_list$mu_ab_d1[j]
  mu_ab_d2 <- mu_ab_list$mu_ab_d2[j]
  
r1 <- NULL
for (i in 1:100){
  out <- draw(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, std10 = std10, ab_50 = ab_50, ab_50_severe, dr_s = dr_s, dr_l = dr_l, t = t, t_d2 = t_d2, k = k, time_to_decay = time_to_decay)
  sub <- data.frame(t = t, run = rep(i, length(t)), Titre = out$titre, Efficacy = out$VE, Efficacy_Severe = out$VE_severe)
  r1 <- rbind(r1, sub)
}

r1 <- r1 %>%
  mutate(Efficacy = Efficacy * 100,
         Efficacy_Severe = Efficacy_Severe * 100) %>%
  pivot_longer(cols = c("Titre", "Efficacy", "Efficacy_Severe"), names_to = "type") %>%
  mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe")))

r1_summary <- r1 %>%
  group_by(type, t) %>%
  summarise(median = median(value),
            upper = quantile(value, 0.975),
            lower = quantile(value, 0.025))

g1 <- ggplot(data = filter(r1, type == "Titre")) +
  geom_line(aes(x = t, y = value, group = run), col = "grey") +
  geom_line(data = filter(r1_summary, type == "Titre"), aes(x = t, y = median), size = 1) +
  geom_ribbon(data = filter(r1_summary, type == "Titre"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
  labs(x = "time (days)", y = "neutralising antibody titre (log 10 scale)", title = mu_ab_list$name[j]) +
  scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
  scale_y_log10(limits = c(1e-3,1e2)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)
  )

g2 <- ggplot(data = filter(r1, type == "Efficacy")) +
  geom_line(aes(x = t, y = value, group = run), col = "grey") +
  geom_line(data = filter(r1_summary, type == "Efficacy"), aes(x = t, y = median), size = 1) +
  geom_ribbon(data = filter(r1_summary, type == "Efficacy"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkblue") +
  lims(y = c(0, 100)) +
  labs(x = "time (days)", y = "vaccine efficacy infection (%)") +
  scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)
  )

g3 <- ggplot(data = filter(r1, type == "Efficacy_Severe")) +
  geom_line(aes(x = t, y = value, group = run), col = "grey") +
  geom_line(data = filter(r1_summary, type == "Efficacy_Severe"), aes(x = t, y = median), size = 1) +
  geom_ribbon(data = filter(r1_summary, type == "Efficacy_Severe"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.4, fill = "darkgrey") +
  lims(y = c(0, 100)) +
  labs(x = "time (days)", y = "vaccine efficacy severe (%)") +
  scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)
  )

plots[[j*3-2]] <- g1
plots[[j*3-1]] <- g2
plots[[j*3]] <- g3
}

plotout <- (plots[[1]] | plots[[2]] | plots[[3]]) / (plots[[4]] | plots[[5]] | plots[[6]]) / (plots[[7]] | plots[[8]] | plots[[9]]) / (plots[[10]] | plots[[11]] | plots[[12]])

ggsave("Fig_2dose_correlated.png", plotout, height = 12, width = 10)

