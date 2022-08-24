library('lmem.sim')
library("lme4")        # model specification / estimation
library("afex")        # anova and deriving p-values from lmer
library("broom.mixed") # extracting data from model fits 
library("faux")        # generate correlated values
library("tidyverse")   # data wrangling and visualisation
library('viridis')

set.seed(1234)

# set up the custom data simulation function
my_sim_data <- function(
    n_subj      = 1280,   # number of subjects
    n_simple  =  12,   # number of complex stimuli
    n_complex =  12,   # number of complex stimuli
    n_small_fam = 8,   #small familiarisation time
    n_medium_fam = 8,  #medium familiarisation time
    n_high_fam = 8,    #high familiarisation time
    n_lab = 40,
    
    beta_0  =  0, # intercept; i.e., the grand mean
    beta_c  =  0.3, # main effect for complexity
    beta_f = 0.3, # main effect for familiarization time
    beta_a = 0.3, # main effect for age
    
    beta_ca = 0.3,
    beta_af = 0.3,
    beta_cf = 0.3,
    
    beta_cfa = 0.3, #main effect for interaction between complexity and familiarization.
    
    item_0 =  0.2, # by-item random intercept sd
    familiarisation_0 =  0.2, # by-familiarisation random intercept sd
    age_0 =  0.2, # by-familiarisation random intercept sd
    
    subject_0   = 0.5, # by-subject random intercept sd
    
    subject_c   =  0.2, # by-subject slope complexity sd
    subject_f = 0.2, # by-subject slope familiarization sd
    subject_a = 0.2, # by-subject slope age sd
    
    subject_ca = 0.2,# by-subject slope for interaction betewen age and complexity sd
    subject_af = 0.2, # by-subject slope for interaction betewen age and familiarisation sd
    subject_cf = 0.2, # by-subject slope complexity*familiarization sd
    
    subject_cfa = 0.2, # by-subject slope for interaction betewen age, complexity and familiarisation sd

    subj_rho     =  .2, # correlations between by-subject random effects
    
    lab_0 = 0.2, # by-lab random intercept sd
    
    lab_c = 0.2, # by-lab slope complexity sd
    lab_f = 0.2, # by-lab slope familiarization sd
    lab_a = 0.2, # by-lab slope age sd
    
    lab_ca = 0.2, # by-lab slope for interaction between age and complexity sd
    lab_af = 0.2, # by-lab slope for interaction between age and familiarisation sd
    lab_cf = 0.2, # by-lab random slope complexity*familiarization sd
    
    lab_cfa = 0.2, # by-lab slope for interaction betewen age, complexity and familiarisation sd
    
    lab_rho = 0.2, # correlations between by-lab random effects
    
    sigma = 1 # residual (error) sd
    ) { # residual (standard deviation)
  
  # simulate a sample of items
  items <- data.frame(
    item_id = seq_len(n_simple + n_complex),
    category = rep(c("simple", "complex"), c(n_simple, n_complex)),
    X_c = rep(c(-0.5, 0.5), c(n_simple, n_complex)),
    O_0i = rnorm(n = n_simple + n_complex, mean = 0, sd = item_0),
    O_0c = rnorm(n = n_simple + n_complex, mean = 0, sd = item_0),
    
    familiarisation = rep(c("small", "medium", "high"), (n_simple + n_complex)/3),
    X_f = rep(c(-0.5, 0, 0.5), (n_simple + n_complex)/3),
    O_0f = rnorm(n = n_simple + n_complex, mean = 0, sd = familiarisation_0)
  )

  # simulate a sample of subjects
  subjects <- faux::rnorm_multi(
    n = n_subj, mu = 0, sd = c(subject_0, 
                               subject_c, 
                               subject_f,
                               subject_a,
                               subject_ca,
                               subject_af,
                               subject_cf,
                               subject_cfa), r = subj_rho,
    varnames = c("S_0", "S_c","S_f","S_a",
                 "S_ca","S_af", "S_cf",
                 "S_cfa")
  ) %>%
    mutate(subj_id = faux::make_id(nrow(.), "S")) %>%
    mutate(X_a = runif(n_subj, min = -0.5, max = 0.5)) #add subject age measure, sample from distribution from -0.5 to 0.5.
  #subjects$subj_id <- 1:n_subj
  
  labs <- faux::rnorm_multi(
    n = n_lab, mu = 0, sd = c(lab_0, lab_c, lab_f, lab_a,
                              lab_ca, lab_af, lab_cf,
                              lab_cfa), r = lab_rho, 
    varnames = c("L_0", "L_c","L_f","L_a",
                 "L_ca","L_af", "L_cf",
                 "L_cfa")
  ) %>%
    mutate(lab_id = faux::make_id(nrow(.), "L"))
  
  #create lab and subj nesting structure
  #Number of subjects must be a multiple of number of labs
  lab_multiplier = n_subj/n_lab
  lab_subj_dict <- data.frame(
    subj_id = subjects$subj_id,
    lab_id = rep(labs$lab_id,lab_multiplier)
  )
  
  
  # cross subject and item IDs 
  temp <- crossing(subjects, items)  %>%
    left_join(lab_subj_dict) %>%
    left_join(labs)
  
  temp %>%
    mutate(
      B_0  = beta_0 + S_0 + L_0 + O_0i,
      
      B_c  = beta_c + S_c + L_c,
      B_f  = beta_f + S_f + L_f,
      B_a = beta_a + S_a + L_a,
      
      B_ca = beta_ca + S_ca + L_ca,
      B_af = beta_af + S_af + L_af,
      B_cf = beta_cf + S_cf + L_cf,
      
      B_cfa = beta_cfa + S_cfa + L_cfa,
      
      e_si = rnorm(nrow(temp), mean = 0, sd = sigma),
      
      DV = B_0 + (B_a * X_a) + (B_c * X_c) + (B_f * X_f) + (B_cf * X_c * X_f) + (B_af * X_a * X_f) + (B_ca * X_c * X_a) + (B_cfa * X_c * X_f * X_a) + e_si
    )
}

dat_sim <- my_sim_data()

#Make plots to see what's going on:

#Main effects:

#Familiarisation:
dat_sim_plot_familiarisation <- dat_sim %>%
  group_by(X_f) %>%
  dplyr::summarise(med_DV = median(DV))

dat_sim %>%
  mutate(X_f = as.factor(X_f)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_f), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = X_f, fill = familiarisation), alpha = 0.2) +
  geom_line(aes(y = med_DV, x = as.factor(X_f), group = 1), data = dat_sim_plot_familiarisation) +
  geom_point(aes(y = med_DV, x = as.factor(X_f)), alpha = 0.8, size = 2, data = dat_sim_plot_familiarisation) +
  scale_fill_manual(values=viridis(n = 3)) +
  ggtitle('Familiarization') +
  theme_bw()

#Complexity:
dat_sim_plot_complexity <- dat_sim %>%
  group_by(X_c) %>%
  dplyr::summarise(med_DV = median(DV))

dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_c), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = X_c, fill = category), alpha = 0.2) +
  geom_line(aes(y = med_DV, x = as.factor(X_c), group = 1), data = dat_sim_plot_complexity) +
  geom_point(aes(y = med_DV, x = as.factor(X_c)), alpha = 0.8, size = 2, data = dat_sim_plot_complexity) +
  scale_fill_manual(values=viridis(n = 2)) +
  ggtitle('Complexity') +
  theme_bw()

#Age:
dat_sim %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_a), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = TRUE, aes(y = DV, x = X_a)) +
  ggtitle('Age') +
  theme_bw()

#Age*familiarisation: 
dat_sim %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_a), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = TRUE, aes(y = DV, x = X_a)) +
  facet_wrap(~X_f) +
  ggtitle('Age x Familiarization Interaction') +
  theme_bw()

#Age*Complexity
dat_sim %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_a), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = TRUE, aes(y = DV, x = X_a)) +
  facet_wrap(~X_c) +
  ggtitle('Age x Complexity Interaction') +
  theme_bw()

#Familiarisation*Complexity
dat_f_c_interaction <- dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  mutate(X_f = as.factor(X_f)) %>%
  group_by(X_f, X_c) %>%
  dplyr::summarise(med_DV = median(DV))

dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  mutate(X_f = as.factor(X_f)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_f), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_point(aes(y = med_DV, x = as.factor(X_f)), alpha = 0.8, size = 2, data = dat_f_c_interaction) +
  geom_line(aes(y = med_DV, x = as.factor(X_f), group = 1), data = dat_f_c_interaction) +
  facet_wrap(~X_c) +
  ggtitle('Familiarization x Complexity Interaction') +
  theme_bw()

run_sims <- function(filename = 'run_sims.csv') {
  
  dat_sim <- my_sim_data()
  
  mod_sim <- lmer(DV ~ 1 + X_a * X_c * X_f + (1 | lab_id / subj_id) + (1 | item_id), data=dat_sim)
  
  sim_results <- broom.mixed::tidy(mod_sim)
  
  # append the results to a file
  append <- file.exists(filename)
  write_csv(sim_results, filename, append = append)
  
  # return the tidy table
  sim_results
}

filename = 'run_sims.csv'
reps <- 100
start_time <- Sys.time()
sims <- purrr::map_df(1:reps, ~run_sims(filename))
end_time <- Sys.time()
end_time - start_time

# read saved simulation data
sims <- read_csv(filename, col_types = cols(
  # makes sure plots display in this order
  group = col_factor(ordered = TRUE),
  term = col_factor(ordered = TRUE)
))
sims

# calculate mean estimates and power for specified alpha
alpha <- 0.05

sim_stats <- sims %>% 
  filter(effect == "fixed") %>%
  group_by(term) %>%
  summarise(
    median_estimate = median(estimate),
    median_se = median(std.error),
    power = mean(p.value < alpha)
  )

# visualise estimates for fixed effects:
sims_fixed <- sims %>%
  filter(effect == "fixed") %>%
  mutate(sim = c(rep(c(1:reps), each = 8)))
                     
sims_fixed %>%
  ggplot(aes(x = sim, y = estimate, ymin = estimate-std.error, ymax = estimate+std.error)) +
  geom_pointrange(fatten = 1/2) +
  facet_wrap(~term, scales = "free_y") +
  theme_bw()

# visualise estimates for random effects:
sim_ran_stats <- sims %>%
  filter(effect == "ran_pars") %>%
  mutate(sim = c(rep(c(1:reps), each = 4)))

sim_ran_stats %>%
  ggplot(aes(x = sim, y = estimate)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~group, scales = "free_y") +
  theme_bw()


#Explore effects of realistic levels of missing data:

run_sims_missing <- function(filename = 'run_sims_missing_data.csv') {
  
  dat_sim <- my_sim_data()
  
  uneven_samples <- dat_sim %>%
    mutate(nas = rbinom(nrow(dat_sim), 1, 1 - .20)) %>%
    mutate(DV = ifelse(nas == 1, DV, NA))
  
  mod_sim <- lmer(DV ~ 1 + X_a * X_c * X_f + (1 | lab_id / subj_id) + (1 | item_id), data=uneven_samples)
  
  sim_results <- broom.mixed::tidy(mod_sim)
  
  # append the results to a file
  append <- file.exists(filename)
  write_csv(sim_results, filename, append = append)
  
  # return the tidy table
  sim_results
}

filename = 'run_sims_missing_data.csv'
reps <- 100
start_time <- Sys.time()
sims_missing <- purrr::map_df(1:reps, ~run_sims_missing(filename))
end_time <- Sys.time()
end_time - start_time


# read saved simulation data
sims <- read_csv(filename, col_types = cols(
  # makes sure plots display in this order
  group = col_factor(ordered = TRUE),
  term = col_factor(ordered = TRUE)
))
sims

# calculate mean estimates and power for specified alpha
alpha <- 0.05

sim_stats <- sims_missing %>% 
  filter(effect == "fixed") %>%
  group_by(term) %>%
  summarise(
    median_estimate = median(estimate),
    median_se = median(std.error),
    power = mean(p.value < alpha)
  )

# visualise estimates for fixed effects:
sims_fixed <- sims_missing %>%
  filter(effect == "fixed") %>%
  mutate(sim = c(rep(c(1:reps), each = 8)))

sims_fixed %>%
  ggplot(aes(x = sim, y = estimate, ymin = estimate-std.error, ymax = estimate+std.error)) +
  geom_pointrange(fatten = 1/2) +
  facet_wrap(~term, scales = "free_y") +
  theme_bw()

# visualise estimates for random effects:
sim_ran_stats <- sims_missing %>%
  filter(effect == "ran_pars") %>%
  mutate(sim = c(rep(c(1:reps), each = 4)))

sim_ran_stats %>%
  ggplot(aes(x = sim, y = estimate)) +
  geom_point(alpha = 0.7) +
  facet_wrap(~group, scales = "free_y") +
  theme_bw()
