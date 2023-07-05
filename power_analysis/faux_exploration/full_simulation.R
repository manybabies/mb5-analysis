library(lmem.sim)
library("lme4")        # model specification / estimation
library("afex")        # anova and deriving p-values from lmer
library("broom.mixed") # extracting data from model fits 
library("faux")        # generate correlated values
library("tidyverse")   # data wrangling and visualisation
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


# set up the custom data simulation function
my_sim_data <- function(
    n_subj      = 1200,   # number of subjects
    n_simple  =  12,   # number of complex stimuli
    n_complex =  12,   # number of complex stimuli
    n_lab = 40,
    beta_0  =  0, # intercept; i.e., the grand mean
    beta_c  =  0.3, # slope for complexity
    beta_f = 0.3, # slope for familiarization time
    item_0 =  1, # by-item random intercept sd
    subject_0   = 1, # by-subject random intercept sd
    subject_c   =  1, # by-subject random slope complexity sd
    subject_f = 1, # by-subject random slope familiarization sd
    subject_cf = 1, # by-subject random slope complexity*familiarization sd
    subj_rho     =  c(.1,.1,.1,.1,.1,.1), # correlations between by-subject random effects
    lab_0   = 1, # by-lab random intercept sd
    lab_c   =  1, # by-lab random slope complexity sd
    lab_f = 1, # by-lab random slope familiarization sd
    lab_cf = 1, # by-lab random slope complexity*familiarization sd
    lab_rho     =  c(.1,.1,.1,.1,.1,.1), # correlations between by-lab random effects
    sigma   = 1 # residual (error) sd
    ) { # residual (standard deviation)
  
  # simulate a sample of items
  items <- data.frame(
    item_id = seq_len(n_simple + n_complex),
    category = rep(c("simple", "complex"), c(n_simple, n_complex)),
    X_i = rep(c(-0.5, 0.5), c(n_simple, n_complex)),
    O_0i = rnorm(n = n_simple + n_complex, mean = 0, sd = item_0)
  )
  
  # simulate a sample of subjects
  subjects <- faux::rnorm_multi(
    n = n_subj, mu = 0, sd = c(subject_0, subject_c,subject_f,subject_cf), r = subj_rho, 
    varnames = c("S_0", "S_c","S_f","S_cf")
  ) %>%
    mutate(subj_id = faux::make_id(nrow(.), "S"))
  #subjects$subj_id <- 1:n_subj
  
  labs <- faux::rnorm_multi(
    n = n_lab, mu = 0, sd = c(lab_0, lab_c,lab_f,lab_cf), r = lab_rho, 
    varnames = c("L_0", "L_c","L_f","L_cf")
  ) %>%
    mutate(lab_id = faux::make_id(nrow(.), "L"))
  
  #create lab and subj nesting structure
  #Number of subjects must be a multiple of number of labs
  lab_multiplier = n_subj/n_lab
  lab_subj_dict <- data_frame(
    subj_id = subjects$subj_id,
    lab_id = rep(labs$lab_id,lab_multiplier)
  )
  
  # cross subject and item IDs 
  temp <- crossing(subjects, items)  %>%
    left_join(lab_subj_dict) %>%
    left_join(labs) %>%
    
  ## LEFT OFF HERE
    mutate(
      e_si = rnorm(nrow(.), mean = 0, sd = sigma),
      DV = beta_0 + T_0s + O_0i + (beta_1 + T_1s) * X_i + e_si
    ) %>%
    select(subj_id, item_id, category, X_i, DV)
}

dat_sim <- my_sim_data()

# fit a linear mixed-effects model to data
mod_sim <- lmer(DV ~ 1 + X_i + (1 | item_id) + (1 + X_i | subj_id),
                data = dat_sim)

summary(mod_sim, corr = FALSE)

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
