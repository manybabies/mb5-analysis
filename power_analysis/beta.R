library(matrixStats)

# The mean of the Oakes dataset is .58 and the SD is .18
.18^2 # .03 variance

# Now we know the mean (.59) and variance (.03), we can simulate the data! 

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

# define parameters
n = 72 # 72 infants in observed data (Oakes dataset)
num_trials = 8
ab <- get.ab(mu = .59, var = .03) # variance refers to the population within-subject variance; mu refers to the population mean. 
alpha <- ab[[1]] #alpha shape parameter for beta distribution
beta <- ab[[2]] #beta shape parameter for beta distribution

# simulate data:
sim_dat <- matrix(rbeta(n = n*num_trials, shape1 = alpha, shape2 = beta), n, num_trials)
hist(sim_dat)
# *Note that this simulates data with zero within-subject consistency (no reliability), but the Oakes data 
# shows zero intercept variance so this might be okay. If we think there is reliability in the pref scores, we would have
# to use a diff. beta data generating model

# impose missing data completely randomly to get closer to number of observations in Oakes data (442)
# note that we don't know if the data is MCAR, but we will assume this for now.

# determine the percent of missing data in the Oakes data
1 - 442/(8*72) # about 23% of the data is missing in Oakes
# Create a random matrix of values ranging from 0-1. Because these values are random, on average 
# 23% of the time the numbers should be less than .23
rand_trials <- matrix(runif(n*num_trials), n, num_trials)
sim_dat[][rand_trials[] < .23] <- NA  
# It's also possible that doing this could wipe an entire baby out, so we'd have to be 
# careful about that 

# Check out the mean and SD of the simulated data:
summary_dat <- data.frame(variance = colVars(sim_dat, na.rm = T), 
                          means = colMeans(sim_dat, na.rm = T), 
                          trial = 1:num_trials)
summary_dat$sd <- sqrt(summary_dat$variance)

hist(sim_dat, main = "simulated data") 
plot(sim_dat) # couple "extreme scores", good

