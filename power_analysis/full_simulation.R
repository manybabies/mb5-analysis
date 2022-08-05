install.packages('lmem.sim')
devtools::install_github("debruine/lmem_sim")

library('lmem.sim')
library("lme4")        # model specification / estimation
library("afex")        # anova and deriving p-values from lmer
library("broom.mixed") # extracting data from model fits 
library("faux")        # generate correlated values
library("tidyverse")   # data wrangling and visualisation
library('viridis')

set.seed(1234)

# set all data-generating parameters
beta_0  <- 0 # intercept; i.e., the grand mean
beta_c  <-  0.3 # slope for complexity
beta_f <- 0.3 # slope for familiarization time
item_0 <-  1 # by-item random intercept sd
subject_0   <- 1 # by-subject random intercept sd
subject_c   <-  1 # by-subject random slope complexity sd
subject_f <- 1 # by-subject random slope familiarization sd
subject_cf <- 1 # by-subject random slope complexity*familiarization sd
subj_rho     <-  c(.1,.1,.1,.1,.1,.1) # correlations between by-subject random effects
lab_0   <- 1 # by-lab random intercept sd
lab_c   <-  1 # by-lab random slope complexity sd
lab_f <- 1 # by-lab random slope familiarization sd
lab_cf <- 1 # by-lab random slope complexity*familiarization sd
lab_rho     <-  c(.1,.1,.1,.1,.1,.1) # correlations between by-lab random effects
sigma   <- 1 # residual (error) sd

# set number of subjects and items
n_subj     <- 1280 # number of subjects
n_lab <- 40 #number of labs
n_simple  <-  12 # number of simple items
n_complex <-  12 # number of complex items
n_small_fam <- 8
n_medium_fam <- 8
n_high_fam <- 8


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
    beta_c  =  0.5, # main effect for complexity
    beta_f = 0.5, # main effect for familiarization time
    beta_a = 0.5, # main effect for age
    
    beta_ca = 0.5,
    beta_af = 0.5,
    beta_cf = 0.5,
    
    beta_cfa = 0.5, #main effect for interaction between complexity and familiarization.
    
    item_0 =  0.2, # by-item random intercept sd
    familiarisation_0 =  0.2, # by-familiarisation random intercept sd
    age_0 =  0.2, # by-familiarisation random intercept sd
    
    subject_0   = 0.5, # by-subject random intercept sd
    
    subject_c   =  0.5, # by-subject slope complexity sd
    subject_f = 0.5, # by-subject slope familiarization sd
    subject_a = 0.5, # by-subject slope age sd
    
    subject_ca = 0.2,# by-subject slope for interaction betewen age and complexity sd
    subject_af = 0.2, # by-subject slope for interaction betewen age and familiarisation sd
    subject_cf = 0.2, # by-subject slope complexity*familiarization sd
    
    subject_cfa = 0.5, # by-subject slope for interaction betewen age, complexity and familiarisation sd

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
    left_join(labs) %>%
    
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

dat_sim %>%
  mutate(X_f = as.factor(X_f)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_f), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = X_f, fill = familiarisation), alpha = 0.2) +
  scale_fill_manual(values=viridis(n = 3)) +
  theme_bw()

dat_sim %>%
  mutate(X_f = as.factor(category)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = category), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = category, fill = category), alpha = 0.2) +
  scale_fill_manual(values=viridis(n = 2)) +
  theme_bw()

dat_sim %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_a), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = TRUE, aes(y = DV, x = X_a)) +
  theme_bw()

dat_sim %>%
  mutate(lab_id = as.factor(lab_id)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = lab_id), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = lab_id, fill = lab_id), alpha = 0.2) +
  scale_fill_manual(values=viridis(n = n_lab)) +
  theme_bw()

dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_a), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_smooth(method = "lm", se = TRUE, aes(y = DV, x = X_a)) +
  facet_wrap(~lab_id) +
  theme_bw()

dat_interaction <- dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  group_by(X_f, X_c) %>%
  dplyr::summarise(med_DV = median(DV))

dat_sim %>%
  mutate(X_c = as.factor(X_c)) %>%
  ggplot() +
  geom_point(aes(y = DV, x = X_c), position = "jitter", alpha = 0.2, size = 0.2) +
  geom_violin(aes(y = DV, x = X_c, fill = X_c), alpha = 0.2) +
  
  geom_point(aes(y = med_DV, x = X_c), alpha = 0.8, size = 2, data = dat_interaction) +
  geom_line(aes(y = med_DV, x = X_c, group = X_f), data = dat_interaction) +
  #scale_fill_manual(values=viridis(n = X_c), data = dat_sim) +
  facet_wrap(~X_f) +
  theme_bw()
  
dat_sim <- my_sim_data()
m1 <- lmer(DV ~ 1 + X_a * X_c * X_f + (1 | lab_id / subj_id) + (1 | item_id), data=dat_sim)
summary(m1, corr = TRUE)
coef(m1)

single_run <- function(filename = NULL, ...) {
  # ... is a shortcut that forwards any additional arguments to my_sim_data()
  dat_sim <- my_sim_data(...)
  mod_sim <- lmer(DV ~ X_i + (1 | item_id) + (1 + X_i | subj_id),
                  dat_sim)
  
  sim_results <- broom.mixed::tidy(mod_sim)
  
  # append the results to a file if filename is set
  if (!is.null(filename)) {
    append <- file.exists(filename) # append if the file exists
    write_csv(sim_results, filename, append = append)
  }
  
  # return the tidy table
  sim_results
}

filename <- "sims.csv" # change for new analyses
if (!file.exists(filename)) {
  # run simulations and save to a file
  reps <- 5
  sims <- purrr::map_df(1:reps, ~single_run(filename))
}

# read saved simulation data
sims <- read_csv(filename)
