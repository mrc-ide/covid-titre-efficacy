library(tidyverse)
library(patchwork)

# parameters

# ab titre params
mu_ab <- 223/94 # pfizer normalised neutralisation titre, linear scale
#mu_ab <- 32/59 # az normalised neutralisation titre, linear scale
ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
ab_50_severe <- 0.03

std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
std <- log(10^std10) # Standard deviation in natural log units
std_lin <- exp(std) # standard deviation on linear scale

k <- 3/log(10) # logistic k value in natural log units (divided by log10 to convert it to natural logarithm units)
logk <- log(k)

t <- 0:(365*4) # Number of days to model
hl <- 108 # Half life of antibody decay
hl_long <- 3650
time_to_decay <- 365 - 250 # time in days to reach longest half-life
dr <- -log(2)/hl # Corresponding decay rate in days for half life above
dr_long <- -log(2)/hl_long

nt <- NULL

draw <- function(mu_ab = mu_ab, std = std, std10 = std10, ab_50 = ab_50, ab_50_severe, rho = rho, dr = dr, t = t, dr_long = dr_long, logk = logk){
  
  # decay rates over time
  dr_vec <- c(rep(dr, 250), seq(dr, dr_long, length.out = time_to_decay), rep(dr_long, (length(t) - time_to_decay - 250)))

  # nab titre distribution for a given vaccine product: draw from a log-normal distribution
  init_titre <- rlnorm(1, log(mu_ab), std)

  # antibody titre in natural log units with decay
  nt[1] <- init_titre
  for (i in 2:length(t)){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  
  # relate titre to efficacy over time
  ef <- 1 / (1 + exp(-logk * (nt - log(ab_50))))
  ef_severe <- 1 / (1 + exp(-logk * (nt - log(ab_50_severe))))

  # output list containing titre and vaccine efficacy
  ret <- list(titre = exp(nt), VE = ef, VE_severe = ef_severe)
  return(ret)
}

r1 <- NULL
for (i in 1:500){
  out <- draw(mu_ab = mu_ab, std = std, std10 = std10, ab_50 = ab_50, ab_50_severe, rho = rho, dr = dr, t = t, dr_long = dr_long, logk = logk)
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
  labs(x = "time (days)", y = "neutralising antibody titre (log 2 scale)") +
  scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
  scale_y_continuous(trans = "log2") +
  theme_bw() +
  theme(strip.background = element_rect(fill = NA, color = "white"),
        panel.border = element_blank(),
        axis.line = element_line(),
        axis.text.x=element_text(angle=60, hjust = 1)
  )

g1

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

g2

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

g3

plots <- (g1 | g2 | g3) + plot_annotation(tag_levels = 'A')

plots

ggsave("Fig1.png", plots, height = 5, width = 13)


# check distribution we are using
x <- rlnorm(10000, log(mu_ab), std)
x <- data.frame(x = x)
ggplot(x) + geom_histogram(aes(x = x)) + scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 24))

# check distribution we are using

x <- rnorm(1000, mu_ab, std_lin)
x <- data.frame(x = x)
ggplot(x) + geom_histogram(aes(x = log2(x))) + scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 24))
