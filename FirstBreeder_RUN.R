# source functions 
source('FirstBreeder_FUNCTIONS.R')

# constants
n.cohorts <- 4
min.age <- 3

# read in data
data <- list()
for (c in 1:n.cohorts)  {
  data[[c]] <- read.csv(paste(getwd(), '/Data/199', c, 'mark.csv', sep=''), 
                        header = T)
}

# number of individuals and capture occasions in each cohort
n <- rep(0, n.cohorts)
K <- rep(0, n.cohorts)
for (c in 1:n.cohorts)  {
  n[c] <- length(data[[c]][,1])
  K[c] <- length(data[[c]][1,])
}



### One normal distribution for recruitment, constant p, constant first-time breeder phi, constant return breeder phi

### model 1
### one normal distribution, mean and sd shared across cohorts
param.1b <- c(0, 0, 0, rep(0, 2))
str.1b <- list('shared', 'shared')
opt.1b <- nlm(likelihood_breeder, param.1b, X = data, arr.dist = 'normal', arr.str = str.1b)
res.1b <- unpack_param_breeder(opt.1b$estimate, 'normal', str.1b, min.age, n.cohorts, K)
AIC.1b <- 2*opt.1b$minimum + 2*length(param.1b) 


### model 2
### one normal distribution, mean cohort dependent, sd shared across cohorts
param.2b <- c(rep(0, n.cohorts), 0, 0, rep(0, 2))
str.2b <- list('cohort', 'shared')
opt.2b <- nlm(likelihood_breeder, param.2b, X = data, arr.dist = 'normal', arr.str = str.2b)
res.2b <- unpack_param_breeder(opt.2b$estimate, 'normal', str.2b, min.age, n.cohorts, K)
AIC.2b <- 2*opt.2b$minimum + 2*length(param.2b)


### model 3
### one normal distribution, mean shared across cohorts, sd cohort dependent
param.3b <- c(0, rep(0, n.cohorts), 0, rep(0, 2))
str.3b <- list('shared', 'cohort')
opt.3b <- nlm(likelihood_breeder, param.3b, X = data, arr.dist = 'normal', arr.str = str.3b)
res.3b <- unpack_param_breeder(opt.3b$estimate, 'normal', str.3b, min.age, n.cohorts, K)
AIC.3b <- 2*opt.3b$minimum + 2*length(param.3b)


### model 4
### one normal distribution, mean and sd cohort dependent
param.4b <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, rep(0, 2))
param.4b <- opt.4b$estimate  # max iterations exceeded, restart at current values
str.4b <- list('cohort', 'cohort')
opt.4b <- nlm(likelihood_breeder, param.4b, X = data, arr.dist = 'normal', arr.str = str.4b)
res.4b <- unpack_param_breeder(opt.4b$estimate, 'normal', str.4b, min.age, n.cohorts, K)
AIC.4b <- 2*opt.4b$minimum + 2*length(param.4b)



### Mixture of two normals for recruitment, constant p, constant first-time breeder phi, constant return breeder phi

### model 5
### two normal distributions, mean, sd and mixture proportion all shared
param.5b <- c(rep(0, 2), rep(0, 2), 0, 0, rep(0, 2))
str.5b <- list(c('shared', 'shared'), c('shared', 'shared'), 'shared')
opt.5b <- nlm(likelihood_breeder, param.5b, X = data, arr.dist = 'normal', arr.str = str.5b)
res.5b <- unpack_param_breeder(opt.5b$estimate, 'normal', str.5b, min.age, n.cohorts, K)
AIC.5b <- 2*opt.5b$minimum + 2*length(param.5b)


### model 6
### two normal distributions, one mean cohort dependent
param.6b <- c(rep(0, n.cohorts), 0, rep(0, 2), 0, 0, rep(0,2))
str.6b <- list(c('cohort', 'shared'), c('shared', 'shared'), 'shared')
opt.6b <- nlm(likelihood_breeder, param.6b, X = data, arr.dist = 'normal', arr.str = str.6b)
res.6b <- unpack_param_breeder(opt.6b$estimate, 'normal', str.6b, min.age, n.cohorts, K)
AIC.6b <- 2*opt.6b$minimum + 2*length(param.6b)


### model 7
### two normal distributions, one sd cohort dependent
param.7b <- c(rep(0, 2), rep(0, n.cohorts), 0, 0, 0, rep(0, 2))
str.7b <- list(c('shared', 'shared'), c('cohort', 'shared'), 'shared')
opt.7b <- nlm(likelihood_breeder, param.7b, X = data, arr.dist = 'normal', arr.str = str.7b)
res.7b <- unpack_param_breeder(opt.7b$estimate, 'normal', str.7b, min.age, n.cohorts, K)
AIC.7b <- 2*opt.7b$minimum + 2*length(param.7b)


### model 8
### two normal distributions, mixing proportion cohort dependent
param.8b <- c(rep(0, 2), rep(0, 2), rep(0, n.cohorts), 0, rep(0, 2))
param.8b <- opt.8b$estimate  # max iterations exceeded, restart at current values
str.8b <- list(c('shared', 'shared'), c('shared', 'shared'), 'cohort')
opt.8b <- nlm(likelihood_breeder, param.8b, X = data, arr.dist = 'normal', arr.str = str.8b)
res.8b <- unpack_param_breeder(opt.8b$estimate, 'normal', str.8b, min.age, n.cohorts, K)
AIC.8b <- 2*opt.8b$minimum + 2*length(param.8b)


### model 9
### two normal distributions, both means cohort dependent
param.9b <- c(rep(0, 2*n.cohorts), rep(0, 2), 0, 0, rep(0, 2))
param.9b <- opt.9b$estimate  # max iterations exceeded, restart at current values
str.9b <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'shared')
opt.9b <- nlm(likelihood_breeder, param.9b, X = data, arr.dist = 'normal', arr.str = str.9b)
res.9b <- unpack_param_breeder(opt.9b$estimate, 'normal', str.9b, min.age, n.cohorts, K)
AIC.9b <- 2*opt.9b$minimum + 2*length(param.9b)


### model 10
### two normal distributions, one normal (mean and sd) cohort dependent
param.10b <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0, 0, rep(0, 2))
str.10b <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'shared')
opt.10b <- nlm(likelihood_breeder, param.10b, X = data, arr.dist = 'normal', arr.str = str.10b)
res.10b <- unpack_param_breeder(opt.10b$estimate, 'normal', str.10b, min.age, n.cohorts, K)
AIC.10b <- 2*opt.10b$minimum + 2*length(param.10b)


### model 11
### two normal distributions, one mean and mixture proportions cohort dependent
param.11b <- c(rep(0, n.cohorts), 0, rep(0, 2), rep(0, n.cohorts), 0, rep(0, 2))
param.11b <- opt.11b$estimate  # max iterations exceeded, restart at current values
str.11b <- list(c('cohort', 'shared'), c('shared', 'shared'), 'cohort')
opt.11b <- nlm(likelihood_breeder, param.11b, X = data, arr.dist = 'normal', arr.str = str.11b)
res.11b <- unpack_param_breeder(opt.11b$estimate, 'normal', str.11b, min.age, n.cohorts, K)
AIC.11b <- 2*opt.11b$minimum + 2*length(param.11b)


### model 12
### two normal distributions, both sds cohort dependent
param.12b <- c(rep(0, 2), rep(0, 2*n.cohorts), 0, 0, rep(0, 2))
param.12b <- opt.12b$estimate  # max iterations exceeded, restart at current values
str.12b <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'shared')
opt.12b <- nlm(likelihood_breeder, param.12b, X = data, arr.dist = 'normal', arr.str = str.12b)
res.12b <- unpack_param_breeder(opt.12b$estimate, 'normal', str.12b, min.age, n.cohorts, K)
AIC.12b <- 2*opt.12b$minimum + 2*length(param.12b)


### model 13
### two normal distributions, one sd and mixture proportions cohort dependent
param.13b <- c(rep(0, 2), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, 2))
param.13b <- opt.13b$estimate  # max iterations exceeded, restart at current values
str.13b <- list(c('shared', 'shared'), c('cohort', 'shared'), 'cohort')
opt.13b <- nlm(likelihood_breeder, param.13b, X = data, arr.dist = 'normal', arr.str = str.13b)
res.13b <- unpack_param_breeder(opt.13b$estimate, 'normal', str.13b, min.age, n.cohorts, K)
AIC.13b <- 2*opt.13b$minimum + 2*length(param.13b)


### model 14
### two normal distributions, both means and one sd cohort dependent
param.14b <- c(rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0, 0, rep(0, 2))
str.14b <- list(c('cohort', 'cohort'), c('cohort', 'shared'), 'shared')
opt.14b <- nlm(likelihood_breeder, param.14b, X = data, arr.dist = 'normal', arr.str = str.14b)
res.14b <- unpack_param_breeder(opt.14b$estimate, 'normal', str.14b, min.age, n.cohorts, K)
AIC.14b <- 2*opt.14b$minimum + 2*length(param.14b)


### model 15
### two normal distributions, both means and mixture proportions cohort dependent
param.15b <- c(rep(0, 2*n.cohorts), rep(0, 2), rep(0, n.cohorts), 0, rep(0, 2))
param.15b <- opt.15b$estimate  # max iterations exceeded, restart at current values
str.15b <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'cohort')
opt.15b <- nlm(likelihood_breeder, param.15b, X = data, arr.dist = 'normal', arr.str = str.15b)
res.15b <- unpack_param_breeder(opt.15b$estimate, 'normal', str.15b, min.age, n.cohorts, K)
AIC.15b <- 2*opt.15b$minimum + 2*length(param.15b)


### model 16
### two normal distributions, one mean and both sds cohort dependent
param.16b <- c(rep(0, n.cohorts), 0, rep(0, 2*n.cohorts), 0, 0, rep(0, 2))
param.16b <- opt.16b$estimate  # max iterations exceeded, restart at current values
str.16b <- list(c('cohort', 'shared'), c('cohort', 'cohort'), 'shared')
opt.16b <- nlm(likelihood_breeder, param.16b, X = data, arr.dist = 'normal', arr.str = str.16b)
res.16b <- unpack_param_breeder(opt.16b$estimate, 'normal', str.16b, min.age, n.cohorts, K)
AIC.16b <- 2*opt.16b$minimum + 2*length(param.16b)


### model 17
### two normal distributions, one normal (mean and sd) and mixture proportions cohort dependent
param.17b <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, 2))
param.17b <- opt.17b$estimate  # max iterations exceeded, restart at current values
str.17b <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'cohort')
opt.17b <- nlm(likelihood_breeder, param.17b, X = data, arr.dist = 'normal', arr.str = str.17b)
res.17b <- unpack_param_breeder(opt.17b$estimate, 'normal', str.17b, min.age, n.cohorts, K)
AIC.17b <- 2*opt.17b$minimum + 2*length(param.17b)


### model 18
### two normal distributions, both sds and mixture proportions cohort dependent
param.18b <- c(rep(0, 2), rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, rep(0, 2))
param.18b <- opt.18b$estimate  # max iterations exceeded, restart at current values
str.18b <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'cohort')
opt.18b <- nlm(likelihood_breeder, param.18b, X = data, arr.dist = 'normal', arr.str = str.18b)
res.18b <- unpack_param_breeder(opt.18b$estimate, 'normal', str.18b, min.age, n.cohorts, K)
AIC.18b <- 2*opt.18b$minimum + 2*length(param.18b)



### One log-normal distribution for recruitment, constant p, constant first-time breeder phi, constant return breeder phi

### model 19
### one log-normal distribution, mean and sd shared across cohorts
param.19b <- c(0, 0, 0, rep(0, 2))
str.19b <- list('shared', 'shared')
opt.19b <- nlm(likelihood_breeder, param.19b, X = data, arr.dist = 'lognormal', arr.str = str.19b)
res.19b <- unpack_param_breeder(opt.19b$estimate, 'lognormal', str.19b, min.age, n.cohorts, K)
AIC.19b <- 2*opt.19b$minimum + 2*length(param.19b)


### model 20
### one log-normal distribution, mean cohort dependent, sd shared across cohorts
param.20b <- c(rep(0, n.cohorts), 0, 0, rep(0, 2))
str.20b <- list('cohort', 'shared')
opt.20b <- nlm(likelihood_breeder, param.20b, X = data, arr.dist = 'lognormal', arr.str = str.20b)
res.20b <- unpack_param_breeder(opt.20b$estimate, 'lognormal', str.20b, min.age, n.cohorts, K)
AIC.20b <- 2*opt.20b$minimum + 2*length(param.20b)


### model 21
### one log-normal distribution, mean shared across cohorts, sd cohort dependent
param.21b <- c(0, rep(0, n.cohorts), 0, rep(0, 2))
str.21b <- list('shared', 'cohort')
opt.21b <- nlm(likelihood_breeder, param.21b, X = data, arr.dist = 'lognormal', arr.str = str.21b)
res.21b <- unpack_param_breeder(opt.21b$estimate, 'lognormal', str.21b, min.age, n.cohorts, K)
AIC.21b <- 2*opt.21b$minimum + 2*length(param.21b)


### model 22
### one log-normal distribution, mean and sd cohort dependent
param.22b <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, rep(0, 2))
str.22b <- list('cohort', 'cohort')
opt.22b <- nlm(likelihood_breeder, param.22b, X = data, arr.dist = 'lognormal', arr.str = str.22b)
res.22b <- unpack_param_breeder(opt.22b$estimate, 'lognormal', str.22b, min.age, n.cohorts, K)
AIC.22b <- 2*opt.22b$minimum + 2*length(param.22b)



### Two log-normal distributions for recruitment, constant p, constant first-time breeder phi, constant return breeder phi

### model 23
### two log-normal distributions, means, sds and mixture proportions all shared
param.23b <- c(rep(0, 2), rep(0, 2), 0, 0, rep(0, 2))
param.23b <- opt.23b$estimate  # max iterations exceeded, restart at current values
str.23b <- list(c('shared', 'shared'), c('shared', 'shared'), 'shared')
opt.23b <- nlm(likelihood_breeder, param.23b, X = data, arr.dist = 'lognormal', arr.str = str.23b)
res.23b <- unpack_param_breeder(opt.23b$estimate, 'lognormal', str.23b, min.age, n.cohorts, K)
AIC.23b <- 2*opt.23b$minimum + 2*length(param.23b)


### model 24
### two log-normal distributions, one mean cohort dependent
param.24b <- c(rep(0, n.cohorts), 0, rep(0, 2), 0, 0, rep(0, 2))
str.24b <- list(c('cohort', 'shared'), c('shared', 'shared'), 'shared')
opt.24b <- nlm(likelihood_breeder, param.24b, X = data, arr.dist = 'lognormal', arr.str = str.24b)
res.24b <- unpack_param_breeder(opt.24b$estimate, 'lognormal', str.24b, min.age, n.cohorts, K)
AIC.24b <- 2*opt.24b$minimum + 2*length(param.24b)


### model 25
### two log-normal distributions, one sd cohort dependent
param.25b <- c(rep(0, 2), rep(0, n.cohorts), 0, 0, 0, rep(0, 2))
str.25b <- list(c('shared', 'shared'), c('cohort', 'shared'), 'shared')
opt.25b <- nlm(likelihood_breeder, param.25b, X = data, arr.dist = 'lognormal', arr.str = str.25b)
res.25b <- unpack_param_breeder(opt.25b$estimate, 'lognormal', str.25b, min.age, n.cohorts, K)
AIC.25b <- 2*opt.25b$minimum + 2*length(param.25b)


### model 26
### two log-normal distributions, mixture proportion cohort dependent
param.26b <- c(rep(0, 2), rep(0, 2), rep(0, n.cohorts), 0, rep(0, 2))
param.26b <- opt.26b$estimate  # max iterations exceeded, restart at current values
str.26b <- list(c('shared', 'shared'), c('shared', 'shared'), 'cohort')
opt.26b <- nlm(likelihood_breeder, param.26b, X = data, arr.dist = 'lognormal', arr.str = str.26b)
res.26b <- unpack_param_breeder(opt.26b$estimate, 'lognormal', str.26b, min.age, n.cohorts, K)
AIC.26b <- 2*opt.26b$minimum + 2*length(param.26b)
CIs.26b <- bootstrap_fn_breeder(999, param.26b, data, 'lognormal', str.26b)


### model 27
### two log-normal distributions, both means cohort dependent
param.27b <- c(rep(0, 2*n.cohorts), rep(0, 2), 0, 0, rep(0, 2))
str.27b <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'shared')
opt.27b <- nlm(likelihood_breeder, param.27b, X = data, arr.dist = 'lognormal', arr.str = str.27b)
res.27b <- unpack_param_breeder(opt.27b$estimate, 'lognormal', str.27b, min.age, n.cohorts, K)
AIC.27b <- 2*opt.27b$minimum + 2*length(param.27b)


### model 28
### two log-normal distributions, one distribution (mean and sd) cohort dependent
param.28b <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0, 0, rep(0, 2))
str.28b <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'shared')
opt.28b <- nlm(likelihood_breeder, param.28b, X = data, arr.dist = 'lognormal', arr.str = str.28b)
res.28b <- unpack_param_breeder(opt.28b$estimate, 'lognormal', str.28b, min.age, n.cohorts, K)
AIC.28b <- 2*opt.28b$minimum + 2*length(param.28b)


### model 29
### two log-normal distributions, one mean and mixture proportions cohort dependent
param.29b <- c(rep(0, n.cohorts), 0, rep(0, 2), rep(0, n.cohorts), 0, rep(0, 2))
str.29b <- list(c('cohort', 'shared'), c('shared', 'shared'), 'cohort')
opt.29b <- nlm(likelihood_breeder, param.29b, X = data, arr.dist = 'lognormal', arr.str = str.29b)
res.29b <- unpack_param_breeder(opt.29b$estimate, 'lognormal', str.29b, min.age, n.cohorts, K)
AIC.29b <- 2*opt.29b$minimum + 2*length(param.29b)


### model 30
### two log-normal distributions, both sds cohort dependent
param.30b <- c(rep(0, 2), rep(0, 2*n.cohorts), 0, 0, rep(0, 2))
param.30b <- opt.30b$estimate  # max iterations exceeded, restart at current values
str.30b <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'shared')
opt.30b <- nlm(likelihood_breeder, param.30b, X = data, arr.dist = 'lognormal', arr.str = str.30b)
res.30b <- unpack_param_breeder(opt.30b$estimate, 'lognormal', str.30b, min.age, n.cohorts, K)
AIC.30b <- 2*opt.30b$minimum + 2*length(param.30b)


### model 31
### two log-normal distributions, one sd and mixture proportions cohort dependent
param.31b <- c(rep(0, 2), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, 2))
str.31b <- list(c('shared', 'shared'), c('cohort', 'shared'), 'cohort')
opt.31b <- nlm(likelihood_breeder, param.31b, X = data, arr.dist = 'lognormal', arr.str = str.31b)
res.31b <- unpack_param_breeder(opt.31b$estimate, 'lognormal', str.31b, min.age, n.cohorts, K)
AIC.31b <- 2*opt.31b$minimum + 2*length(param.31b)
