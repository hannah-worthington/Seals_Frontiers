### Code to run simulation study

# Name: sim.study
# Objective: To generate data and compare the multiple-cohort to single-cohort approach for estimating recruitment
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
#         n.sim - number of simulations to run
# Outputs: multiple - parameter estimates for the multiple-cohort model
#          single - parameter estimates for the single-cohort models

sim_study <- function(n, C, K, r, p, phi, arr.dist, n.sim) {
  
  # storage
  beta.all <- list()
  p.all <- rep(0, n.sim)
  phi.all <- rep(0, n.sim)
  beta.1 <- list()
  p.1 <- list()
  phi.1 <- list()
  
  # constants
  str.param <- list('cohort', 'cohort')
  param.all <- c(rep(0, 2*C), rep(0, 2*C), rep(0, C), 0, 0)
  param.1 <- c(rep(0, 2), rep(0, 2), 0, 0, 0)
  
  # loop over iterations
  for (sim in 1: n.sim)  {
  
    # simulate data
    data <- sim_data(n, C, K, r, p, phi, arr.dist)
    
    # multiple-cohort approach
    opt.all <- nlm(likelihood_allcohorts, param.all, X = data, arr.dist = arr.dist, arr.str = str.param)
    while (opt.all$iterations == 100)  {
      opt.all <- nlm(likelihood_allcohorts, opt.all$estimate, X = data, arr.dist = arr.dist, arr.str = str.param)
    }
    res.all <- unpack_param(opt.all$estimate, arr.dist, str.param, min.age = 3, n.cohorts, K)
    for (c in 1:n.cohorts)  {
      beta.all[[c]] <- rbind(beta.all[[c]], res.all$beta[[c]])
    } 
    p.all[sim] <- res.all$values$p
    phi.all[sim] <- res.all$values$phi
  
    # single-cohort approach
    for (c in 1:n.cohorts)  {
      opt.1 <- nlm(likelihood_singlecohort, param.1, X = data[[c]], arr.dist = arr.dist, arr.str = str.param)
      while (opt.1$iterations == 100)  {
        opt.1 <- nlm(likelihood_singlecohort, opt.1$estimate, X = data[[c]], arr.dist = arr.dist, arr.str = str.param)
      }
      res.1 <- unpack_param(opt.1$estimate, arr.dist, str.param, min.age = 3, n.cohorts, K)
      beta.1[[c]] <- rbind(beta.1[[c]], res.1$beta)
      p.1[[c]] <- rbind(p.1[[c]], res.1$values$p)
      phi.1[[c]] <- rbind(phi.1[[c]], res.1$values$phi)
    }
  }
  
  # return the parameter estimates for both approaches
  return(list('multiple' = list('beta' = beta.all, 'p' = p.all, 'phi' = phi.all), 'single' = list('beta' = beta.1, 'p' = p.1, 'phi' = phi.1)))
}
  


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

sim_data <- function(n, C, K, r, p, phi, arr.dist)  {
  
  # create storage variables
  X <- list()
  
  # loop over cohorts
  for (c in 1:C)  {
    
    # number of capture occasions for cohort
    K.c <- K - c + 1
    
    # storage variable
    attendance <- matrix(0, nrow = n, ncol = K.c) 
    
    # generate recruitment probabilities
    beta_fun <- match.fun(paste('two', arr.dist, 'betas', sep = ''))
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



### Likelihood function

# Name: likelihood_allcohorts
# Objective: To evaluate the negative log-likelihood of the multiple cohort HMM given a set of parameter values and capture histories
# Inputs: param - model specific, will be passed to unpack_param function for structuring
#         X - capture histories, a matrix for each cohort stored in a list
#         arr.dist - distribution on arrivals, 'normal' or 'lognormal'
#         arr.str - structure for arrivals (shared or cohort), list over parameters of chosen distribution
#         min.age - minimum age of return (default = 3)
# Outputs: lik - negative log-likelihood value

likelihood_allcohorts <- function(param, X, arr.dist, arr.str, min.age = 3)  {
  
  # print(param)
  
  # define constants
  n.cohorts <- length(X)  # number of cohorts
  n <- rep(0, n.cohorts)  # number of individuals in each cohort
  K <- rep(0, n.cohorts)  # number of capture occasions for each cohort
  for (c in 1:n.cohorts)  {
    n[c] <- length(X[[c]][,1])
    K[c] <- length(X[[c]][1,])
  }
  
  # unpack the parameter vector
  params <- unpack_param(param, arr.dist, arr.str, min.age, n.cohorts, K)
  HMM <- HMM_str(params, n.cohorts, K)
  
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
      # add to likelihood for this cohort if nonzero
      if (contrib > 0)  {
        lik.partial <- lik.partial + log(contrib)
      } else if (contrib == 0)  {
        lik.partial <- lik.partial - 10000
      }
    }
    # store likelihood value for cohort, checking for nonzero
    lik.cohort[c] <- lik.partial
  }
  
  # evaluate full likelihood
  lik <- -sum(lik.cohort)
  
  # return
  return(lik)
}



### Likelihood function

# Name: likelihood_singlecohort
# Objective: To evaluate the negative log-likelihood of a single cohort HMM given a set of parameter values and capture histories
# Inputs: param - model specific, will be passed to unpack_param function for structuring
#         X - capture histories, a matrix
#         arr.dist - distribution on arrivals, 'normal' or 'lognormal'
#         arr.str - structure for arrivals (shared or cohort), list over parameters of chosen distribution
#         min.age - minimum age of return (default = 3)
# Outputs: lik - negative log-likelihood value

likelihood_singlecohort <- function(param, X, arr.dist, arr.str, min.age = 3)  {
  
  # print(param)
  
  # define constants
  n <- length(X[,1])  # number of individuals in each cohort
  K <- length(X[1,]) # number of capture occasions for each cohort
  
  # unpack the parameter vector
  params <- unpack_param(param, arr.dist, arr.str, min.age, n.cohorts = 1, K)
  HMM <- HMM_str(params, n.cohorts = 1, K)
  
  lik <- 0
  # loop over individuals
  for (i in 1:n)  {
    contrib <- HMM$pione[[1]]
    # loop over occasions
    for (k in 1:K)  {
      contrib <- contrib %*% HMM$gamma[[1]][,,k]
      contrib <- contrib %*% HMM$pmat[[1]][,,X[i,k]+1, k]
    }
    # individual contribution
    contrib <- contrib %*% matrix(1, nrow = 3, ncol = 1)
    # add to likelihood if nonzero
    if (contrib > 0)  {
      lik <- lik + log(contrib)
    } else if (contrib == 0)  {
      lik <- lik - 10000
    }
  }
  
  # return
  return(lik)
}



### Function to return the HMM structures

# Name: HMM_str
# Objective: To take a list of parameters returned from the arr.dist functions and return the HMM structures for use in the likelihood
# Inputs: param - list of parameters generated by the unpack_param function
#         n.cohorts - the number of cohorts
#         K - number of capture occasions for each cohort
# Outputs: pione - initial state distributions for each cohort
#          gamma - transition probability matrices for each cohort
#          pmat - observation matrices for each cohort

HMM_str <- function(param, n.cohorts, K)  {
  
  # initial state probabilities
  # state probabilities at time 0, year of tagging
  pione <- list()
  for (c in 1:n.cohorts)  {
    pione[[c]] <- c(1, 0, 0)
  }
  
  # transition probability matrices (state, state, occasion)
  gamma <- list()
  for (c in 1:n.cohorts)  {
    gamma[[c]] <- array(0, dim = c(3, 3, K[c]))
    gamma[[c]][1,1,] <- 1 - param$betastar[[c]]
    gamma[[c]][1,2,] <- param$betastar[[c]]
    gamma[[c]][2,2,] <- param$phi
    gamma[[c]][2,3,] <- 1 - param$phi
    gamma[[c]][3,3,] <- 1
  }
  
  # observation matrices (state, state, capture/not, occasion)
  pmat <- list()
  for (c in 1:n.cohorts)  {
    pmat[[c]] <- array(0, dim = c(3, 3, 2, K[c]))
    pmat[[c]][1,1,1,] <- 1
    pmat[[c]][2,2,1,] <- 1-param$p[[c]]
    pmat[[c]][3,3,1,] <- 1
    pmat[[c]][2,2,2,] <- param$p[[c]]
  }
  
  # return
  return(list('pione' = pione, 'gamma' = gamma, 'pmat' = pmat))
}



### Function to unpack the parameter vector
### Options: normal distribution(s) on arrivals, constant capture, constant survival
### Options: log-normal distribution(s) on arrivals, constant capture, constant survival

# Name: unpack_param
# Objective: To unpack and transform a vector of parameter values for either the cohort or single stopover HMM model
# Inputs: param - vector of parameter values, model dependent:
#               - For normal distribution(s): arrival distribution parameters (mean(s) (log), sd(s) (log), w (mixture parameter(s), logit))
#               - For log-normal distribution(s): arrival distribution parameters (mean(s), sd(s) (log), w (mixture parameter(s), logit))
#               - capture probability (logit)
#               - survival probability (logit)
#         arr.dist - distribution, 'normal' or 'lognormal'
#         arr.str - arrivals structure, shared or cohort specific (list with vector for each of mean, sd, w (when required) for both normal and log-normal)
#         min.age - minimum age of return
#         n.cohorts - number of cohorts
#         K - number of capture occasions for each cohort
# Outputs: mu - mean(s) of arrival distribution for each cohort
#          sd - sd(s) of arrival distribution for each cohort
#          w - mixture proportion for first mixture (1 for single distribution)
#          beta - arrival probabilities for each cohort
#          betastar - conditional arrival probabilities for each cohort
#          p - capture probability, time structured to permit structural zeros
#          phi- survival probability

unpack_param <- function(param, arr.dist, arr.str, min.age, n.cohorts, K)  {
  
  # number of mixtures taken from length of first object in arr.str list
  n.mixtures <- length(arr.str[[1]])
  
  # storage
  means <- matrix(0, nrow = n.mixtures, ncol = n.cohorts)
  sds <- matrix(0, nrow = n.mixtures, ncol = n.cohorts)
  w <- matrix(0, nrow = n.mixtures - 1, ncol = n.cohorts)
  res <- list()
  names <- c()
  
  # mean(s) of arrival distributions for normal distribution(s)
  if (arr.dist == 'normal')  {
    for (m in 1:n.mixtures)  {
      names <- c(names, paste('mu', m, sep = ''))
      if (arr.str[[1]][m] == 'shared')  {
        mean <- exp(param[1])
        means[m,] <- rep(mean, n.cohorts)
        param <- param[-1]
      } else if (arr.str[[1]][m] == 'cohort')  {
        mean <- exp(param[1:n.cohorts])
        means[m,] <- mean
        param <- param[-(1:n.cohorts)]
      }
      res <- c(res, list(mean))
    }
  }
  
  # mean(s) of arrival distributions for log-normal distribution(s)
  if (arr.dist == 'lognormal')  {
    for (m in 1:n.mixtures)  {
      names <- c(names, paste('mu', m, sep = ''))
      if (arr.str[[1]][m] == 'shared')  {
        mean <- param[1]
        means[m,] <- rep(mean, n.cohorts)
        param <- param[-1]
      } else if (arr.str[[1]][m] == 'cohort')  {
        mean <- param[1:n.cohorts]
        means[m,] <- mean
        param <- param[-(1:n.cohorts)]
      }
      res <- c(res, list(mean))
    }
  }
  
  # sd of arrival distributions
  for (m in 1:n.mixtures)  {
    names <- c(names, paste('sd', m, sep = ''))
    if (arr.str[[2]][m] == 'shared')  {
      sd <- exp(param[1])
      sds[m,] <- rep(sd, n.cohorts)
      param <- param[-1]
    } else if (arr.str[[2]][m] == 'cohort')  {
      sd <- exp(param[1:n.cohorts])
      sds[m,] <- sd
      param <- param[-(1:n.cohorts)]
    }
    res <- c(res, list(sd))
  }
  
  # mixture proportions if needed
  if (n.mixtures >= 2)  {
    for (m in 1:(n.mixtures - 1))  {
      names <- c(names, paste('w', m, sep = ''))
      if (arr.str[[3]][m] == 'shared')  {
        mix <- 1/(1 + exp(-param[1]))
        w[m,] <- rep(mix, n.cohorts)
        param <- param[-1]
      } else if (arr.str[[3]][m] == 'cohort')  {
        mix <- 1/(1 + exp(-param[1:n.cohorts]))
        w[m,] <- mix
        param <- param[-(1:n.cohorts)]
      }
      res <- c(res, list(mix))
    }
  } else if (n.mixtures == 1)  {
    w <- 1
  }
  
  # arrival probabilities
  beta <- list()
  text.cohorts <- ifelse(n.mixtures == 1, 'one', 'two')
  beta_fun <- match.fun(paste(text.cohorts, arr.dist, 'betas', sep=''))
  for (c in 1:n.cohorts)  {
    if (n.mixtures == 1)  {
      beta[[c]] <- beta_fun(mu = means[,c], sd = sds[,c], K = K[c], min.age = min.age)
    } else if (n.mixtures == 2)  {
      beta[[c]] <- beta_fun(mu = means[,c], sd = sds[,c], w = w[,c], K = K[c], min.age = min.age)
    }
  }
  
  # conditional arrival probabilities
  betastar <- list()
  for (c in 1:n.cohorts)  {
    betastar[[c]] <- rep(0, K[c])
    for (k in 1:K[c])  {
      if (sum(beta[[c]][k:K[c]]) > 0)  {
        betastar[[c]][k] <- beta[[c]][k]/sum(beta[[c]][k:K[c]])
      }
    }
  }
  
  # capture probability (constant)
  pnonzero <- 1/(1 + exp(-param[1]))
  param <- param[-1]
  names <- c(names, 'p')
  res <- c(res, list(pnonzero))
  p <- list()
  for (c in 1:n.cohorts)  {
    p[[c]] <- rep(0, K[c])
    p[[c]][min.age:K[c]] <- pnonzero
  }
  
  # retention probability (constant)
  phi <- 1/(1 + exp(-param))
  names <- c(names, 'phi')
  res <- c(res, list(phi))
  
  # add names to parameter list
  names(res) <- names
  
  # return all parameters
  return(c(list('mu' = means, 'sd' = sds, 'w' = w, 'beta' = beta, 'betastar' = betastar, 'p' = p, 'phi' = phi), 'values' = list(res)))
}




### Functions to calculate arrival probabilities

### Function to calculate beta probabilities from a mixture of two normal distributions over the plausible recruitment ages

# Name: twonormalbetas
# Objective: To calculate the beta probabilities from a mixture of two truncated normal distributions
# Inputs: mu - means vector for the two normals that form the arrival distribution
#         sd - sds vector for the two normals that form the arrival distribution
#         w - mixture proportion for distribution 1
#         K - number of occasions
#         min.age - minimum age of return
# Outputs: beta - set of beta parameters

twonormalbetas <- function(mu, sd, w, K, min.age)  {
  
  # find the separate truncated mixture distributions
  mix.1 <- onenormalbetas(mu[1], sd[1], K, min.age)
  mix.2 <- onenormalbetas(mu[2], sd[2], K, min.age)
  
  # create mixture
  beta <- w*mix.1 + (1-w)*mix.2
  
  # check the probabilities sum to one
  total <- sum(beta)
  if (total > 0)  {
    beta <- beta/total
  } 
  
  # return the beta probabilities
  return(beta)
}



### Function to calculate beta probabilities from a mixture of two log-normal distributions over the plausible recruitment ages

# Name: twolognormalbetas
# Objective: To calculate the beta probabilities from a mixture of two truncated log-normal distributions
# Inputs: mu - means vector for the two log-normals that form the arrival distribution
#         sd - sds vector for the two log-normals that form the arrival distribution
#         w - mixture proportion for distribution 1
#         K - number of occasions
#         min.age - minimum age of return
# Outputs: beta - set of beta parameters

twolognormalbetas <- function(mu, sd, w, K, min.age)  {
  
  # find the separate truncated mixture distributions
  mix.1 <- onelognormalbetas(mu[1], sd[1], K, min.age)
  mix.2 <- onelognormalbetas(mu[2], sd[2], K, min.age)
  
  # create mixture
  beta <- w*mix.1 + (1-w)*mix.2
  
  # check the probabilities sum to one
  total <- sum(beta)
  if (total > 0)  {
    beta <- beta/total
  }
  
  # return the beta probabilities
  return(beta)
}


