library(tidyverse)
library(patchwork)

# parameters
variant_fold_reduction <- 1
dose_3_fold_increase <- 6

mu_ab_list <- data.frame(name = c("AstraZeneca", "Pfizer", "Moderna"),
                         mu_ab_d1 = c(32/59/18.82,223/94/18.82,654/158/18.82),            
                         mu_ab_d2 = c(32/59/12,223/94/12,654/158/12)) %>%        #32/59,223/94,654/158    

  mutate(mu_ab_d1 = mu_ab_d1/variant_fold_reduction,
         mu_ab_d2 = mu_ab_d2/variant_fold_reduction) %>%
  mutate(mu_ab_d3 = mu_ab_d2 * dose_3_fold_increase)

ab_50 <- 0.2 # titre relative to convalescent required to provide 50% protection from infection, on linear scale
ab_50_severe <- 0.03
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 data
k <- 2.94 # shape parameter of efficacy curve
max_t <- 730 # number of days to model
t_d2 <- 84 # timing of second dose relative to first
t_d3 <- 240 # timing of third dose relative to second
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
draw <- function(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, mu_ab_d3, std10 = std10, ab_50 = ab_50, ab_50_severe, 
                 dr_s = dr_s, dr_l = dr_l, t = t, t_d2 = t_d2, t_d3, k = k, time_to_decay = time_to_decay){
  
  # decay rates over time
  dr_vec_sub <- c(rep(dr_s, period_s),
                  seq(dr_s, dr_l, length.out = time_to_decay),
                  rep(dr_l, (length(t) - t_period_l)))
  
  dr_vec <- rep(0, length(t))
  dr_vec[1:(t_d2-1)] <- dr_vec_sub[1:(t_d2-1)]
#  dr_vec[t_d2:(t_d3+t_d2-1)] <- dr_vec_sub[t_d2:(t_d3+t_d2-1)]
  dr_vec[t_d2:(t_d3+t_d2-1)] <- dr_vec_sub[1:((t_d3+t_d2-1)-t_d2+1)]
  dr_vec[(t_d2+t_d3):length(t)] <- dr_vec_sub[1:(length(t) - (t_d2+t_d3) + 1) ]
  
  # nab titre distribution for a given vaccine product: draw from a log-normal distribution
  # here the first titre is dependent on the second but I suppose could do it the other way around

  z2 <- rnorm(1, log10(mu_ab_d2), std10)
  z3 <- log10(10^z2 * (mu_ab_d3/mu_ab_d2))
  z1 <- log10(10^z2 * (mu_ab_d1/mu_ab_d2))
  
  # initiate titre vector
  nt <- rep(0, length(t))
  nt[1] <- log(10^z1)
  nt[t_d2] <- log(10^z2)
  nt[(t_d2+t_d3)] <- log(10^z3)
  
  # decay antibodies over time on natural log scale
  for (i in (2:(t_d2-1))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  for (i in ((t_d2+1):(t_d2+t_d3-1))){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  for (i in ((t_d2+t_d3+1):length(t))){
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

r1_summary <- NULL
plots <- NULL
plotlist <- list()
for (j in 1:3){
  mu_ab_d1 <- mu_ab_list$mu_ab_d1[j]
  mu_ab_d2 <- mu_ab_list$mu_ab_d2[j]
  mu_ab_d3 <- mu_ab_list$mu_ab_d3[j]
  
  r1 <- NULL
  for (i in 1:100){
    out <- draw(mu_ab_d1 = mu_ab_d1, mu_ab_d2 = mu_ab_d2, mu_ab_d3 = mu_ab_d3, std10 = std10, ab_50 = ab_50, ab_50_severe, 
                dr_s = dr_s, dr_l = dr_l, t = t, t_d2 = t_d2, t_d3, k = k, time_to_decay = time_to_decay)
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
  
  
  
  if(j==1){
  g1 <- ggplot(data = filter(r1, type == "Titre")) +
    geom_line(aes(x = t, y = value, group = run), col = "grey") +
    geom_line(data = filter(r1_summary, type == "Titre"), aes(x = t, y = median), size = 1,col = "black") +
    geom_ribbon(data = filter(r1_summary, type == "Titre"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[29], y = median[113]*0.13), color = "green", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[85], y = median[113]*0.1), color = "green", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[113], y = median[113]*1), color = "green", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[325], y = median[113]*0.28), color = "green", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[339], y = median[113]*1.04), color = "green", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Titre"),aes(x = t[353], y = median[113]*1.81), color = "green", size = 2) +
    labs(x = "time (days)", y = "neutralising antibody titre (log 10 scale)", title = mu_ab_list$name[j]) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_log10(limits = c(1e-3,1e2)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, color = "white"),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_text(angle=60, hjust = 1))
    } else {
    g1 <- ggplot(data = filter(r1, type == "Titre")) +
        geom_line(aes(x = t, y = value, group = run), col = "grey") +
        geom_line(data = filter(r1_summary, type == "Titre"), aes(x = t, y = median), size = 1,col = "black") +
        geom_ribbon(data = filter(r1_summary, type == "Titre"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "darkgreen") +
        labs(x = "time (days)", y = "neutralising antibody titre (log 10 scale)", title = mu_ab_list$name[j]) +
        scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
        scale_y_log10(limits = c(1e-3,1e2)) +
        theme_bw() +
        theme(strip.background = element_rect(fill = NA, color = "white"),
              panel.border = element_blank(),
              axis.line = element_line(),
              axis.text.x=element_text(angle=60, hjust = 1))}
  
  
  
  
  if(j==1){
  g2 <- ggplot(data = filter(r1, type == "Efficacy")) +
    geom_line(aes(x = t, y = value, group = run), col = "grey") +
    geom_line(data = filter(r1_summary, type == "Efficacy"), aes(x = t, y = median), size = 1) +
    geom_ribbon(data = filter(r1_summary, type == "Efficacy"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "yellow") +
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[21], y = 43.3), color = "blue", size = 2) +
    #
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+10], y = 62.7), color = "blue", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+42], y = 66.7), color = "blue", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+87], y = 59.3), color = "blue", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+122], y =52.6), color = "blue", size = 2) +
    # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+160], y =47.3), color = "blue", size = 2) +
    #
    geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[84], y =40), color = "red", size = 2) +
    lims(y = c(0, 100)) +
    labs(x = "time (days)", y = "vaccine efficacy transmission (%)") +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA, color = "white"),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.text.x=element_text(angle=60, hjust = 1))}
  else if(j==2){
    g2 <- ggplot(data = filter(r1, type == "Efficacy")) +
      geom_line(aes(x = t, y = value, group = run), col = "grey") +
      geom_line(data = filter(r1_summary, type == "Efficacy"), aes(x = t, y = median), size = 1) +
      geom_ribbon(data = filter(r1_summary, type == "Efficacy"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "yellow") +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[21], y = 51.9), color = "blue", size = 2) +
      #
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+10], y = 92.4), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+42], y = 89.8), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+87], y = 80.3), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+122], y = 73.4), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+160], y = 69.7), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[1], y = 50), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[85], y = 84), color = "blue", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+4], y = 60), color = "purple", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+19], y = 65), color = "purple", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+45], y = 48), color = "purple", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+75], y = 30), color = "purple", size = 2) +
      # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+100], y = 0), color = "purple", size = 2) +
      #
      geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[84], y =40), color = "red", size = 2) +
      lims(y = c(0, 100)) +
      labs(x = "time (days)", y = "vaccine efficacy transmission (%)") +
      scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
      theme_bw() +
      theme(strip.background = element_rect(fill = NA, color = "white"),
            panel.border = element_blank(),
            axis.line = element_line(),
            axis.text.x=element_text(angle=60, hjust = 1))}
  else if(j==3) {
      g2 <- ggplot(data = filter(r1, type == "Efficacy")) +
        geom_line(aes(x = t, y = value, group = run), col = "grey") +
        geom_line(data = filter(r1_summary, type == "Efficacy"), aes(x = t, y = median), size = 1) +
      geom_ribbon(data = filter(r1_summary, type == "Efficacy"), aes(x = t, ymin = lower, ymax = upper), alpha = 0.2, fill = "yellow") +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[21], y = 65.9), color = "blue", size = 2) +
        #
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+10], y = 95.2), color = "blue", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+42], y = 94.5), color = "blue", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+87], y = 90.3), color = "blue", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[1], y = 79), color = "blue", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[85], y = 86), color = "blue", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+4], y = 65), color = "purple", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+19], y = 80), color = "purple", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+45], y = 70), color = "purple", size = 2) +
        # geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[83+75], y = 65), color = "purple", size = 2) +
        #
        geom_point(data = filter(r1_summary, type == "Efficacy"),aes(x = t[84], y =40), color = "red", size = 2) +
        lims(y = c(0, 100)) +
        labs(x = "time (days)", y = "vaccine efficacy transmission (%)") +
        scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
        theme_bw() +
        theme(strip.background = element_rect(fill = NA, color = "white"),
              panel.border = element_blank(),
              axis.line = element_line(),
              axis.text.x=element_text(angle=60, hjust = 1))}

  plots[[j*2-1]] <- g1
  plots[[j*2]] <- g2
}

plotout <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) / (plots[[5]] | plots[[6]])
plotout
ggsave("Transmission-Blocking.png", plotout, height = 10, width = 13)