library(lmem.sim)
library("lme4")        # model specification / estimation
library("afex")        # anova and deriving p-values from lmer
library("broom.mixed") # extracting data from model fits 
library("faux")        # generate correlated values
library("tidyverse")   # data wrangling and visualisation
set.seed(1234)
# set all data-generating parameters
beta_0  <- 0 # intercept; i.e., the grand mean
beta_1  <-  0.3 # slope; i.e, effect of category
omega_0 <-  1 # by-item random intercept sd
tau_0   <- 1 # by-subject random intercept sd
tau_1   <-  1 # by-subject random slope sd
rho     <-  .2 # correlation between intercept and slope ##explore
sigma   <- 1 # residual (error) sd

# set number of subjects and items
n_subj     <- 1280 # number of subjects
n_simple  <-  12 # number of simple items
n_complex <-  12 # number of complex items


# set up the custom data simulation function
my_sim_data <- function(
    n_subj      = 1280,   # number of subjects
    n_simple  =  12,   # number of complex stimuli
    n_complex =  12,   # number of complex stimuli
    beta_0     = 0,   # grand mean
    beta_1     =  0.3,   # effect of complexity
    omega_0    =  1,   # by-item random intercept sd
    tau_0      = 1,   # by-subject random intercept sd
    tau_1      =  1,   # by-subject random slope sd
    rho        = 0.2,   # correlation between intercept and slope
    sigma      = 1) { # residual (standard deviation)
  
  # simulate a sample of items
  items <- data.frame(
    item_id = seq_len(n_simple + n_complex),
    category = rep(c("simple", "complex"), c(n_simple, n_complex)),
    X_i = rep(c(-0.5, 0.5), c(n_simple, n_complex)),
    O_0i = rnorm(n = n_simple + n_complex, mean = 0, sd = omega_0)
  )
  
  # simulate a sample of subjects
  subjects <- faux::rnorm_multi(
    n = n_subj, mu = 0, sd = c(tau_0, tau_1), r = rho, 
    varnames = c("T_0s", "T_1s")
  )
  subjects$subj_id <- 1:n_subj
  
  # cross subject and item IDs 
  crossing(subjects, items)  %>%
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
