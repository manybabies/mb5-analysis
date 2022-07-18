Sys.setenv(DOWNLOAD_STATIC_LIBV8 = 1) # only necessary for Linux without the nodejs library / headers

install.packages('Rcpp', type = "source")
library(Rcpp)

install.packages('ggplot2', version = "3.3.5")
library(ggplot2)

options(mc.cores = parallel::detectCores())

install.packages('tidybayes')
library(tidybayes)

install.packages('pacman')

pacman::p_load(
  #general
  tidyverse,
  here,
  dplyr,
  moments,
  tidybayes,
  #visualization
  ggplot2,
  ggridges,
  phonTools,
  phonR,
  factoextra,
  plyr,
  ellipse,
  viridis,
  ggtext,
  #modelfitting
  brms,
  mclust,
  #correlationplot
  corrplot,
  qgraph,
  tidyverse,
  corrr,
  igraph,
  ggraph)

## Simulate Datasets
# To use the rbeta function, we need to supply two shape functions that correspond to the mean and the 
# variance of the data. To get those shape parameters, we will use the get.ab function below from 
# Mijke Rhemtulla/ Minhajuddin et al. (2004)

get.ab <- function(mu, var){
  v <- var
  w <- mu/(1-mu)
  b <- ((w/ (v*(w^2 + 2*w + 1))) - 1)  / (w + 1)
  a <- b*w
  return(list(a = a, b = b))
}

set.seed(123) # reproducible sampling

generate_dataset <- function(n_labs=40, n_per_lab=32, 
                             effect_sizes=list(type = .3, 
                                               familiarization = .1,
                                               age = .1, "age*type"=.1,
                                               "age*familiarization"=0, "type*familiarization"=0,
                                               "type*age*familiarization"=.1)) { 
  # critical test is the 3-way interaction?
  
  # rewrite to use expand.grid ?
  labID = rep(as.character(1:n_labs), each=n_per_lab)
  subjID = 1:(n_labs*n_per_lab)
  
  familiarization_times = c(-0.5,0,0.5) # or maybe we expect linear effect on log(fam_time)?
  fam_times_sc = c(-0.5,0,0.5) # scaled
  #fam_times_sc = log(familiarization_times) # tried this: yields much lower power for main effects, only slight benefit for interactions
  stimulus_types = c(rep("high",4), rep("low",4)) # stimulus complexity
  # trials each subject gets (but randomly ordered)
  fam_by_stim = expand.grid(fam_time = fam_times_sc, stimulus_type = stimulus_types)
  
  # assume each lab uses one procedure
  lab_procedure = sample(c("IC","FD"), n_labs, replace=T, prob=c(.5,.5)) # 50/50 IC / FD procedures?
  procedure = rep(lab_procedure, each=n_per_lab)
  
  test_order = rep(1:4, n_per_lab/4*n_labs) 
  
  # per-subject data
  simd <- tibble(subjID, labID, procedure, test_order) %>%
    mutate(subjInt = rnorm(length(subjID), mean=0, sd=1))
  
  # add lab random intercept
  simd$labInt = 0.0
  for(lab in unique(labID)) {
    labInd = which(simd$labID==lab)
    simd[labInd,]$labInt = rnorm(1, mean=0, sd=1) # could increase per-lab variability ..
  }
  
  # uniform random vars
  simd$age_mos = runif(nrow(simd), min=3.0, max=15.0)
  simd$age = scale(simd$age_mos, center=T, scale=T)[,1]
  
  # generate per-subject data, put in long (row per-trial) df
  
  siml <- tibble()
  for(i in 1:nrow(simd)) {
    # randomized trial order (but maybe should be done according to preset pseudorandom orders?)
    tmp_sdat <- fam_by_stim[sample(1:nrow(fam_by_stim), size=nrow(fam_by_stim), replace=F),]
    # let's assume prop_novel is normally-distributed
    stimulus_type = with(tmp_sdat, ifelse(stimulus_type=="high", .5, -.5)) 
    error_term = rnorm(nrow(tmp_sdat), 0, sd=1) + simd[i,]$labInt + simd[i,]$subjInt # add random slopes? (e.g. by age..)
    # rescale error to be >0
    # ToDo: scale familiarization time ?
    age_effect_subj = effect_sizes$age * rep(simd[i,]$age, nrow(tmp_sdat))
    
    # can we assume these are z-scored proportions of novel looking? maybe truncate them?
    # ToDo: check if problems when effect sizes are 0?
    tmp_sdat$dv_zscore = effect_sizes$type * stimulus_type + # main
      age_effect_subj +  # main
      effect_sizes$familiarization * tmp_sdat$fam_time +  # main
      effect_sizes$`age*type` * stimulus_type * effect_sizes$type * age_effect_subj + 
      effect_sizes$`age*familiarization` * age_effect_subj * tmp_sdat$fam_time * effect_sizes$familiarization + 
      effect_sizes$`type*familiarization` * tmp_sdat$fam_time * stimulus_type * effect_sizes$type + 
      effect_sizes$`type*age*familiarization` * stimulus_type * effect_sizes$type * age_effect_subj * tmp_sdat$fam_time * effect_sizes$familiarization + 
      error_term
    
    # since DV has SD~.2, and must be in the range [0,1], let's make floor=-2.5 and ceiling=2.5. # the next four lines generate a     # strange distribution: 
    #min_ind = which(tmp_sdat$dv_zscore< -2.5)
    #max_ind = which(tmp_sdat$dv_zscore > 2.5)
    #if(length(min_ind)>0) tmp_sdat[min_ind,]$dv_zscore = -2.5
    #if(length(max_ind)>0) tmp_sdat[max_ind,]$dv_zscore = 2.5
    
    siml <- siml %>% 
      bind_rows(tmp_sdat %>% mutate(subjID = simd[i,]$subjID,
                                    labID = simd[i,]$labID,
                                    age = simd[i,]$age,
                                    age_mos = simd[i,]$age_mos,
                                    subjInt = simd[i,]$subjInt,
                                    labInt = simd[i,]$labInt,
                                    trial_num = 1:nrow(tmp_sdat)))
    #novel_looking_time = rnorm(n = nrow(tmp_sdat), mean=0, sd=1), # = .05
    #familiar_looking_time = rnorm(n = nrow(tmp_sdat), mean=0, sd=1), # = .05
    #prop_novel = novel_looking_time / (novel_looking_time + familiar_looking_time), # use beta distribution?
    #prop_novel = rbeta(n=nrow(tmp_sdat), shape1=??, shape2=??)
    # mean_beta = .5 + familiarization_time*age*type
    # how to choose beta parameters: more non-central = more of a novelty/familiarity effect
  }
  
  siml$trial_num_sc = scale(siml$trial_num, center=T, scale=T) 
  
  siml$subjID = as.factor(siml$subjID)
  # switch from dummy-code to effects code 
  siml$stimulus_type = as.factor(siml$stimulus_type)
  contrasts(siml$stimulus_type) = c(0.5, -0.5)
  
  return(siml)
}


effect_size_pt3 = list(type = .3, 
                       familiarization = .3, 
                       age = .3, 
                       "age*type"=.3, 
                       "age*familiarization"=.3, 
                       "type*familiarization"=.3, 
                       "type*age*familiarization"=.3)

siml = generate_dataset(effect_sizes = effect_size_pt3)

siml %>%
  ggplot() +
  geom_density(aes(dv_zscore)) +
  theme_bw()

#save generated dataset
write.csv(siml, "example_power_analysis_data_0.3.csv")

model_formula <- bf(dv_zscore ~ 1 + stimulus_type * fam_time * age + (1 | subjID) + (1 | labID))

get_prior(model_formula,
          data = siml, 
          family = gaussian)

priors1 <- c(prior(normal(0, 0.5), class = Intercept),
             prior(normal(0, 0.5), class = b),
             prior(normal(1, 1), class = sd),
             prior(normal(1, 1), class = sigma))

fit <- brm(data = siml,
           family = gaussian,
           model_formula,
           prior = priors1,
           sample_prior = "yes",
           iter = 3000,
           warmup = 500,
           cores = 64,
           chains = 2,
           save_pars = save_pars(all = TRUE))

pp_check(fit, ndraws = 100)
plot(conditional_effects(fit), points = T)


#Sample the parameters of interest:
Posterior_m1 <- as_draws_df(fit)
#Plot the prior-posterior update plot for the intercept:
ggplot(Posterior_m1) +
  geom_density(aes(prior_Intercept), fill="steelblue", color="black",alpha=0.6) +
  geom_density(aes(b_Intercept), fill="#FC4E07", color="black",alpha=0.6) + 
  theme_classic()

ggplot(Posterior_m1) +
  geom_density(aes(prior_sigma), fill="steelblue", color="black",alpha=0.6) +
  geom_density(aes(sigma), fill="#FC4E07", color="black",alpha=0.6) + 
  theme_classic()

ggplot(Posterior_m1) +
  geom_density(aes(prior_sd_Subject), fill="steelblue", color="black",alpha=0.6) +
  geom_density(aes(sd_Subject), fill="#FC4E07", color="black",alpha=0.6) + 
  geom_density(aes(sd_labID), fill="#FC4E07", color="black",alpha=0.6) + 
  theme_classic()

ggplot(Posterior_m1) +
  geom_density(aes(prior_b), fill="steelblue", color="black",alpha=0.6) +
  geom_density(aes(b_age), fill="#FC4E07", color="black",alpha=0.6) + 
  geom_density(aes(b_fam_time), fill="#FC4E07", color="black",alpha=0.6) + 
  geom_density(aes(b_stimulus_type), fill="#FC4E07", color="black",alpha=0.6) + 
  theme_classic()


siml = generate_dataset(effect_sizes = effect_size_pt3)

sim_d_and_fit <- function(seed) {
  
  set.seed(seed)
  
  siml = generate_dataset(effect_sizes = effect_size_pt3)
  
  update(fit,
         newdata = siml,
         seed = seed) %>% 
    fixef() %>% 
    data.frame()
}

glimpse(siml)

n_sim = 30

s3 <-
  tibble(seed = 1:n_sim) %>%
  mutate(b1 = purrr::map(seed, sim_d_and_fit)) %>%
  unnest(b1)

fixedeff <- fixef(fit)

model_parameters <- s3 %>%
  mutate(parameters = rep(rownames(fixedeff), n_sim)) %>%
  print(n=50)

write.csv(model_parameters, 'model_parameters_thirty_sims.csv')

model_parameters <- read.csv('/Users/au620441/Downloads/model_parameters_thirty_sims.csv')

unique(model_parameters$parameters)

intercept_analysis <- model_parameters %>% 
  filter(parameters == "Intercept") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) + 
  ggtitle('intercept') + theme_bw()
intercept_analysis<- intercept_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
intercept_analysis

stimulus_type1_analysis <- model_parameters %>% 
  filter(parameters == "stimulus_type1") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('stim_type') + theme_bw()
stimulus_type1_analysis <- stimulus_type1_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
stimulus_type1_analysis

fam_time_analysis <- model_parameters %>% 
  filter(parameters == "fam_time") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('fam_time') + theme_bw()
fam_time_analysis <- fam_time_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
fam_time_analysis

age_analysis <- model_parameters %>% 
  filter(parameters == "age") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('age') + theme_bw()
age_analysis <- age_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
age_analysis

stim_fam_analysis <- model_parameters %>% 
  filter(parameters == "stimulus_type1:fam_time") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('stim_type:fam_time') + theme_bw()
stim_fam_analysis <- stim_fam_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
stim_fam_analysis


stim_age_analysis <- model_parameters %>% 
  filter(parameters == "stimulus_type1:age") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('stim_type:age') + theme_bw()
stim_age_analysis <- stim_age_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
stim_age_analysis

fam_age_analysis <- model_parameters %>% 
  filter(parameters == "fam_time:age") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('fam_time:age') + theme_bw()
fam_age_analysis <- fam_age_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
fam_age_analysis

stim_fam_age_analysis <- model_parameters %>% 
  filter(parameters == "stimulus_type1:fam_time:age") %>%
  ggplot(aes(x = reorder(seed, Q2.5), y = Estimate, ymin = Q2.5, ymax = Q97.5)) +
  geom_pointrange(fatten = 1/2) +
  labs(x = "simulated dataset",
       y = expression(beta[1])) +
  ggtitle('stim_type:fam_time:age') + theme_bw()
stim_fam_age_analysis <- stim_fam_age_analysis +
  theme(plot.title = element_text(hjust = 0.5, size=20))
stim_fam_age_analysis

siml %>%
  ggplot() +
  geom_density(aes(age)) +
  theme_bw()

precision_analysis <- cowplot::plot_grid(intercept_analysis, stimulus_type1_analysis,
                                         fam_time_analysis, age_analysis, 
                                         stim_fam_analysis, stim_age_analysis, fam_age_analysis,
                                         stim_fam_age_analysis,
                                         ncol=2)

unique(model_parameters$parameters)

model_parameters_width <- model_parameters %>% 
  mutate(width = Q97.5 - Q2.5) %>%
  filter(parameters == "stimulus_type1:fam_time")

model_parameters_width %>% 
  mutate(check = ifelse(width < .1, 1, 0)) %>% 
  summarise(`proportion below 0.7` = mean(check),
            `average width`        = mean(width))


ggsave(plot = precision_analysis, file = "precision_analysis.png", width = 15, height = 10)