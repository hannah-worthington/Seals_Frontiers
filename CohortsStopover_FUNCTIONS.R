### Likelihood function

### Function to evaluate the negative log-likelihood of the multiple cohort stopover model in the HMM framework

# Name: likelihood_cohort
# Objective: To evaluate the negative log-likelihood of the multiple cohort HMM given a set of parameter values and capture histories
# Inputs: param - model specific, will be passed to other functions for structuring
#         X - capture histories, a matrix for each cohort stored in a list
#         model - model to be fitted
# Outputs: lik - negative log-likelihood value

likelihood_cohort<- function(param, X, model)  {
  
  # define constants
  n.cohorts <- length(X)  # number of cohorts
  n <- rep(0, n.cohorts)  # number of individuals in each cohort
  K <- rep(0, n.cohorts)  # number of capture occasions for each cohort
  for (c in 1:n.cohorts)  {
    n[c] <- length(data[[c]][,1])
    K[c] <- length(data[[c]][1,])
  }
  
  # unpack the parameter vector
  param_function <- match.fun(model)
  params <- param_function(param, n.cohorts, n, K)
  HMM <- HMM.str(params, n.cohorts, K)
  
  # storage variables
  lik.cohort <- rep(0, n.cohorts)
  
  # loop over cohorts evaluating likelihood for single cohort
  for (c in 1:n.cohorts)  {
    lik.partial <- 0
    # loop over individuals
    for (i in 1:n[c])  {
      contrib <- HMM$pione[[c]]
      # loop over occasions
      for (k in 1:K[c])  {
        contrib <- contrib %*% HMM$gamma[[c]][,,k]
        contrib <- contrib %*% HMM$pmat[[c]][,,X[[c]][i,k]+1, k]
      }
      # individual contribution
      contrib <- contrib %*% matrix(1, nrow = 3, ncol = 1)
      # add to likelihood for this cohort
      lik.partial <- lik.partial + log(contrib)
    }
    # store likelihood value for cohort
    lik.cohort[c] <- lik.partial
  }
  
  # evaluate full likelihood
  lik <- -sum(lik.cohort)
  
  # return
  return(lik)
}
  
  
  
  
### Function to unpack a parameter vector - normal distribution on arrivals, cohort mean and sd, constant capture, constant survival

# Name: onenormalalldiff_constp_constphi
# Objective: To unpack and transform a vector of parameter values for the multiple cohort stopover model HMM
# Inputs: param - vector of parameter values
#               - [n.cohorts] - log mu for arrival distribution, one per cohort
#               - [n.cohorts] - log sd for arrival distribution, one per cohort
#               - [1] - logit p, constant capture probability, shared across all cohorts 
#               - [1] - logit phi, survival probability, shared across all cohorts
#         n.cohorts - number of cohorts
#         n - vector of observed number of individuals in each cohort
#         K - vector of the number of capture occasions for each cohort
#         min.age - minimum age of return
# Outputs: mu - means for each cohort
#          sd - sd for each cohort
#          beta - arrival probabilities for each cohort
#          betastar - conditional arrival probabilities for each cohort
#          p - capture probabilities
#          phi - survival probabilities

onenormalalldiff_constp_constphi <- function(param, n.cohorts, n, K, min.age = 3)  {
  
  # means
  normal.means <- exp(param[1:n.cohorts])
  
  # sds
  normal.sds <- exp(param[(n.cohorts+1):(2*n.cohorts)])
  
  # arrival probabilities
  beta <- list()
  for (c in 1:n.cohorts)  {
    beta[[c]] <- onenormalbetas(normal.means[c], normal.sds[c], K[c], min.age)
  }
  
  # conditional arrival probabilities
  betastar
  }
}