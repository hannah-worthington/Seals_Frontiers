### Code to generate multi-cohort capture-recapture data

# Name: sim.data
# Objective: To generate multiple cohorts of capture-recapture data using the given parameter values
# Inputs: n - number of individuals in a cohort (constant)
#         C - number of cohorts
#         K - number of occasions for cohort 0 (assume consecutive cohorts)
#         r - list of parameters for the recruitment distributions, distribution dependent
#           - for normal (means, sds, w)
#           - for log-normal (means, sds, w)
#           - for Poisson (means, w)
#         p - capture probability (constant) 
#         phi - survival probability (constant) 
#         arr.dist - arrival distribution (normal, lognormal, poisson)
# Outputs: X - list of capture history matrices across the cohorts

sim.data <- function(n, C, K, r, p, phi, arr.dist)  {
  
  # create storage variables
  X <- list()
  
  # loop over cohorts
  for (c in 1:C)  {
    
    # number of capture occasions for cohort
    K.c <- K - c + 1
    
    # storage variable
    attendance <- matrix(0, nrow = n, ncol = K.c) 
    
    # generate recruitment probabilities
    beta_fun <- match.fun(paste('two', arr.dist, 'betas', sep= ''))
    if (arr.dist == 'normal' | arr.dist == 'lognormal')  {
      beta <- beta_fun(r$means[,c], r$sds[,c], r$w[c], K.c)
    } else if (arr.dist == 'poisson')
      beta <- beta_fun(r$means[,c], r$w[c], K.c)
    
    # generate recruitment occasions
    recruit.freq <- rmultinom(1, n, beta)
    recruit <- rep(1:K.c, times = recruit.freq)
    
    # loop over individuals
    for (i in 1:n)  {
      
      # record years of attendance
      attendance[i, recruit[i]] <- 1
      if (recruit[i] < K.c)  {
        for (k in (recruit[i] + 1):K.c)  {
          if (attendance[i, (k-1)] == 1)  {
            attendance[i,k] <- rbinom(1, 1, phi)
          }
        }
      }
    }
    
    # determine capture histories
    capture <- ifelse(attendance == 1, rbinom(1, 1, p), 0)
    
    # add capture matrix to list
    X <- append(X, capture)
  }
  
  # return
  return(list('X' = X))
}



