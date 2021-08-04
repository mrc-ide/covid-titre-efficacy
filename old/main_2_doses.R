library(tidyverse)
library(patchwork)

# parameters
mu_ab <- 223/94 # pfizer normalised neutralisation titre, linear scale
mu_ab <- 32/59 # az normalised neutralisation titre, linear scale
mu_ab <- 28/164 # sinovac

mu_ab_list <- data.frame(name = c("Pfizer", "AstraZeneca", "Sinovac"),
                         titre = c(223/94, 32/59, 28/164))

mu_ab_d1 <- 13/94
mu_ab_d2 <- 223/94
ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
ab_50_severe <- 0.03
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
std <- log(10^std10) # Standard deviation in natural log units
k <- 2.94
t <- 0:(365*2) # Number of days to model
t_d2 <- 84

hl <- 108 # Half life of antibody decay
hl_long <- 3650
time_to_decay <- 365 - 250 # time in days to reach longest half-life
dr <- -log(2)/hl # Corresponding decay rate in days for half life above
dr_long <- -log(2)/hl_long

nt <- NULL

draw <- function(mu_ab = mu_ab, std10 = std10, ab_50 = ab_50, ab_50_severe, rho = rho, dr = dr, t = t, dr_long = dr_long, k = k){
  
  # decay rates over time
  dr_vec_sub <- c(rep(dr, 250), seq(dr, dr_long, length.out = time_to_decay), rep(dr_long, (length(t) - time_to_decay - 250)))
  
  dr_vec <- rep(0, length(t))
  dr_vec[1:(t_d2-1)] <- dr_vec_sub[1:(t_d2-1)]
  dr_vec[t_d2:length(t)] <- dr_vec_sub[1:(length(t) - t_d2 + 1) ]

  # nab titre distribution for a given vaccine product: draw from a log-normal distribution
  z1 <- rnorm(1, log10(mu_ab_d1), std10)
  z2 <- rnorm(1, log10(mu_ab_d2), std10)

  nt <- rep(0, length(t))
  nt[1] <- log(10^z1)
  nt[t_d2] <- log(10^z2)
  

  for (i in (2:(t_d2-1))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  for (i in ((t_d2+1):length(t))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  nt <- exp(nt)
  
  # relate titre to efficacy over time
  ef <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50))))
  
  ef_severe <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_severe))))

  # output list containing titre and vaccine efficacy
  ret <- list(titre = nt, VE = ef, VE_severe = ef_severe)
  return(ret)
}

  
r1 <- NULL
for (i in 1:100){
  out <- draw(mu_ab = mu_ab, std10 = std10, ab_50 = ab_50, ab_50_severe, rho = rho, dr = dr, t = t, dr_long = dr_long, k = k)
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
  labs(x = "time (days)", y = "neutralising antibody titre (log 10 scale)") +
  scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
  scale_y_log10(limits = c(1e-3,1e2)) +
  #scale_y_continuous(trans = "log10", limits = c(0, 100)) +
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

plotout <- (plots[[1]] | plots[[2]] | plots[[3]]) / (plots[[4]] | plots[[5]] | plots[[6]]) / (plots[[7]] | plots[[8]] | plots[[9]])

plotout

ggsave("Fig.png", plotout, height = 12, width = 10)


# check distribution we are using
x <- rlnorm(10000, log(mu_ab), std)
x <- data.frame(x = x)
ggplot(x) + geom_histogram(aes(x = x)) + scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 24))

# check distribution we are using

x <- rnorm(1000, mu_ab, std_lin)
x <- data.frame(x = x)
ggplot(x) + geom_histogram(aes(x = log2(x))) + scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 24))
