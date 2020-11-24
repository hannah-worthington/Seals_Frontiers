# source functions
source('CohortsStopover_FUNCTIONS.R')

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



### One normal distribution for recruitment, constant p, constant phi

### model 1
### one normal distribution, mean and sd shared across cohorts
param.1 <- c(0, 0, 0, 0)
str.1 <- list('shared', 'shared')
opt.1 <- nlm(likelihood_cohort, param.1, X = data, arr.dist = 'normal', arr.str = str.1)
res.1 <- unpack_param(opt.1$estimate, 'normal', str.1, min.age, n.cohorts, K)
AIC.1 <- 2*opt.1$minimum + 2*length(param.1) 


### model 2
### one normal distribution, mean cohort dependent, sd shared across cohorts
param.2 <- c(rep(0, n.cohorts), 0, 0, 0)
str.2 <- list('cohort', 'shared')
opt.2 <- nlm(likelihood_cohort, param.2, X = data, arr.dist = 'normal', arr.str = str.2)
res.2 <- unpack_param(opt.2$estimate, 'normal', str.2, min.age, n.cohorts, K)
AIC.2 <- 2*opt.2$minimum + 2*length(param.2)


### model 3
### one normal distribution, mean shared across cohorts, sd cohort dependent
param.3 <- c(0, rep(0, n.cohorts), 0, 0)
str.3 <- list('shared', 'cohort')
opt.3 <- nlm(likelihood_cohort, param.3, X = data, arr.dist = 'normal', arr.str = str.3)
res.3 <- unpack_param(opt.3$estimate, 'normal', str.3, min.age, n.cohorts, K)
AIC.3 <- 2*opt.3$minimum + 2*length(param.3)


### model 4
### one normal distribution, mean and sd cohort dependent
param.4 <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, 0)
param.4 <- opt.4$estimate # max iteration exceeded, restart at current values
str.4 <- list('cohort', 'cohort')
opt.4 <- nlm(likelihood_cohort, param.4, X = data, arr.dist = 'normal', arr.str = str.4)
res.4 <- unpack_param(opt.4$estimate, 'normal', str.4, min.age, n.cohorts, K)
AIC.4 <- 2*opt.4$minimum + 2*length(param.4)



### Mixture of two normals for recruitment, constant p, constant phi

### model 5
### two normal distributions, mean, sd and mixture proportion all shared
param.5 <- c(rep(0, 2), rep(0, 2), 0, 0, 0)
str.5 <- list(c('shared', 'shared'), c('shared', 'shared'), 'shared')
opt.5 <- nlm(likelihood_cohort, param.5, X = data, arr.dist = 'normal', arr.str = str.5)
res.5 <- unpack_param(opt.5$estimate, 'normal', str.5, min.age, n.cohorts, K)
AIC.5 <- 2*opt.5$minimum + 2*length(param.5)


### model 6
### two normal distributions, one mean cohort dependent
param.6 <- c(rep(0, n.cohorts), 0, rep(0, 2), 0, 0, 0)
str.6 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'shared')
opt.6 <- nlm(likelihood_cohort, param.6, X = data, arr.dist = 'normal', arr.str = str.6)
res.6 <- unpack_param(opt.6$estimate, 'normal', str.6, min.age, n.cohorts, K)
AIC.6 <- 2*opt.6$minimum + 2*length(param.6)


### model 7
### two normal distributions, one sd cohort dependent
param.7 <- c(rep(0, 2), rep(0, n.cohorts), 0, 0, 0, 0)
param.7 <- opt.7$estimate # max iterations exceeded, restart at current values
str.7 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'shared')
opt.7 <- nlm(likelihood_cohort, param.7, X = data, arr.dist = 'normal', arr.str = str.7)
res.7 <- unpack_param(opt.7$estimate, 'normal', str.7, min.age, n.cohorts, K)
AIC.7 <- 2*opt.7$minimum + 2*length(param.7)


### model 8
### two normal distributions, mixing proportion cohort dependent
param.8 <- c(rep(0, 2), rep(0, 2), rep(0, n.cohorts), 0, 0)
param.8 <- opt.8$estimate # max iterations exceeded, restart at current values
str.8 <- list(c('shared', 'shared'), c('shared', 'shared'), 'cohort')
opt.8 <- nlm(likelihood_cohort, param.8, X = data, arr.dist = 'normal', arr.str = str.8)
res.8 <- unpack_param(opt.8$estimate, 'normal', str.8, min.age, n.cohorts, K)
AIC.8 <- 2*opt.8$minimum + 2*length(param.8)


### model 9
### two normal distributions, both means cohort dependent
param.9 <- c(rep(0, 2*n.cohorts), rep(0, 2), 0, 0, 0)
param.9 <- opt.9$estimate  # max iterations exceeded, restart at current values
str.9 <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'shared')
opt.9 <- nlm(likelihood_cohort, param.9, X = data, arr.dist = 'normal', arr.str = str.9)
res.9 <- unpack_param(opt.9$estimate, 'normal', str.9, min.age, n.cohorts, K)
AIC.9 <- 2*opt.9$minimum + 2*length(param.9)


### model 10
### two normal distributions, one normal (mean and sd) cohort dependent
param.10 <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0, 0, 0)
str.10 <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'shared')
opt.10 <- nlm(likelihood_cohort, param.10, X = data, arr.dist = 'normal', arr.str = str.10)
res.10 <- unpack_param(opt.10$estimate, 'normal', str.10, min.age, n.cohorts, K)
AIC.10 <- 2*opt.10$minimum + 2*length(param.10)


### model 11
### two normal distributions, one mean and mixture proportions cohort dependent
param.11 <- c(rep(0, n.cohorts), 0, rep(0, 2), rep(0, n.cohorts), 0, 0)
param.11 <- opt.11$estimate  # max iterations exceeded, restart at current values
str.11 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'cohort')
opt.11 <- nlm(likelihood_cohort, param.11, X = data, arr.dist = 'normal', arr.str = str.11)
res.11 <- unpack_param(opt.11$estimate, 'normal', str.11, min.age, n.cohorts, K)
AIC.11 <- 2*opt.11$minimum + 2*length(param.11)


### model 12
### two normal distributions, both sds cohort dependent
param.12 <- c(rep(0, 2), rep(0, 2*n.cohorts), 0, 0, 0)
param.12 <- opt.12$estimate  # max iterations exceeded, restart at current values
str.12 <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'shared')
opt.12 <- nlm(likelihood_cohort, param.12, X = data, arr.dist = 'normal', arr.str = str.12)
res.12 <- unpack_param(opt.12$estimate, 'normal', str.12, min.age, n.cohorts, K)
AIC.12 <- 2*opt.12$minimum + 2*length(param.12)


### model 13
### two normal distributions, one sd and mixture proportions cohort dependent
param.13 <- c(rep(0, 2), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
param.13 <- opt.13$estimate  # max iterations exceeded, restart at current values
str.13 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'cohort')
opt.13 <- nlm(likelihood_cohort, param.13, X = data, arr.dist = 'normal', arr.str = str.13)
res.13 <- unpack_param(opt.13$estimate, 'normal', str.13, min.age, n.cohorts, K)
AIC.13 <- 2*opt.13$minimum + 2*length(param.13)


### model 14
### two normal distributions, both means and one sd cohort dependent
param.14 <- c(rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0, 0, 0)
str.14 <- list(c('cohort', 'cohort'), c('cohort', 'shared'), 'shared')
opt.14 <- nlm(likelihood_cohort, param.14, X = data, arr.dist = 'normal', arr.str = str.14)
res.14 <- unpack_param(opt.14$estimate, 'normal', str.14, min.age, n.cohorts, K)
AIC.14 <- 2*opt.14$minimum + 2*length(param.14)


### model 15
### two normal distributions, both means and mixture proportions cohort dependent
param.15 <- c(rep(0, 2*n.cohorts), rep(0, 2), rep(0, n.cohorts), 0, 0)
param.15 <- opt.15$estimate  # max iterations exceeded, restart at current values
str.15 <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'cohort')
opt.15 <- nlm(likelihood_cohort, param.15, X = data, arr.dist = 'normal', arr.str = str.15)
res.15 <- unpack_param(opt.15$estimate, 'normal', str.15, min.age, n.cohorts, K)
AIC.15 <- 2*opt.15$minimum + 2*length(param.15)


### model 16
### two normal distributions, one mean and both sds cohort dependent
param.16 <- c(rep(0, n.cohorts), 0, rep(0, 2*n.cohorts), 0, 0, 0)
param.16 <- opt.16$estimate  # max iterations exceeded, restart at current values
str.16 <- list(c('cohort', 'shared'), c('cohort', 'cohort'), 'shared')
opt.16 <- nlm(likelihood_cohort, param.16, X = data, arr.dist = 'normal', arr.str = str.16)
res.16 <- unpack_param(opt.16$estimate, 'normal', str.16, min.age, n.cohorts, K)
AIC.16 <- 2*opt.16$minimum + 2*length(param.16)


### model 17
### two normal distributions, one normal (mean and sd) and mixture proportions cohort dependent
param.17 <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
param.17 <- opt.17$estimate  # max iterations exceeded, restart at current values
str.17 <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'cohort')
opt.17 <- nlm(likelihood_cohort, param.17, X = data, arr.dist = 'normal', arr.str = str.17)
res.17 <- unpack_param(opt.17$estimate, 'normal', str.17, min.age, n.cohorts, K)
AIC.17 <- 2*opt.17$minimum + 2*length(param.17)


### model 18
### two normal distributions, both sds and mixture proportions cohort dependent
param.18 <- c(rep(0, 2), rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0)
param.18 <- opt.18$estimate  # max iterations exceeded, restart at current values
str.18 <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'cohort')
opt.18 <- nlm(likelihood_cohort, param.18, X = data, arr.dist = 'normal', arr.str = str.18)
res.18 <- unpack_param(opt.18$estimate, 'normal', str.18, min.age, n.cohorts, K)
AIC.18 <- 2*opt.18$minimum + 2*length(param.18)



### One log-normal distribution for recruitment, constant p, constant phi

### model 19
### one log-normal distribution, mean and sd shared across cohorts
param.19 <- c(0, 0, 0, 0)
str.19 <- list('shared', 'shared')
opt.19 <- nlm(likelihood_cohort, param.19, X = data, arr.dist = 'lognormal', arr.str = str.19)
res.19 <- unpack_param(opt.19$estimate, 'lognormal', str.19, min.age, n.cohorts, K)
AIC.19 <- 2*opt.19$minimum + 2*length(param.19)


### model 20
### one log-normal distribution, mean cohort dependent, sd shared across cohorts
param.20 <- c(rep(0, n.cohorts), 0, 0, 0)
str.20 <- list('cohort', 'shared')
opt.20 <- nlm(likelihood_cohort, param.20, X = data, arr.dist = 'lognormal', arr.str = str.20)
res.20 <- unpack_param(opt.20$estimate, 'lognormal', str.20, min.age, n.cohorts, K)
AIC.20 <- 2*opt.20$minimum + 2*length(param.20)


### model 21
### one log-normal distribution, mean shared across cohorts, sd cohort dependent
param.21 <- c(0, rep(0, n.cohorts), 0, 0)
str.21 <- list('shared', 'cohort')
opt.21 <- nlm(likelihood_cohort, param.21, X = data, arr.dist = 'lognormal', arr.str = str.21)
res.21 <- unpack_param(opt.21$estimate, 'lognormal', str.21, min.age, n.cohorts, K)
AIC.21 <- 2*opt.21$minimum + 2*length(param.21)


### model 22
### one log-normal distribution, mean and sd cohort dependent
param.22 <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, 0)
str.22 <- list('cohort', 'cohort')
opt.22 <- nlm(likelihood_cohort, param.22, X = data, arr.dist = 'lognormal', arr.str = str.22)
res.22 <- unpack_param(opt.22$estimate, 'lognormal', str.22, min.age, n.cohorts, K)
AIC.22 <- 2*opt.22$minimum + 2*length(param.22)



### Two log-normal distributions for recruitment, constant p, constant phi

### model 23
### two log-normal distributions, means, sds and mixture proportions all shared
param.23 <- c(rep(0, 2), rep(0, 2), 0, 0, 0)
str.23 <- list(c('shared', 'shared'), c('shared', 'shared'), 'shared')
opt.23 <- nlm(likelihood_cohort, param.23, X = data, arr.dist = 'lognormal', arr.str = str.23)
res.23 <- unpack_param(opt.23$estimate, 'lognormal', str.23, min.age, n.cohorts, K)
AIC.23 <- 2*opt.23$minimum + 2*length(param.23)


### model 24
### two log-normal distributions, one mean cohort dependent
param.24 <- c(rep(0, n.cohorts), 0, rep(0, 2), 0, 0, 0)
str.24 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'shared')
opt.24 <- nlm(likelihood_cohort, param.24, X = data, arr.dist = 'lognormal', arr.str = str.24)
res.24 <- unpack_param(opt.24$estimate, 'lognormal', str.24, min.age, n.cohorts, K)
AIC.24 <- 2*opt.24$minimum + 2*length(param.24)


### model 25
### two log-normal distributions, one sd cohort dependent
param.25 <- c(rep(0, 2), rep(0, n.cohorts), 0, 0, 0, 0)
param.25 <- opt.25$estimate  # max iterations exceeded, restart at current values
str.25 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'shared')
opt.25 <- nlm(likelihood_cohort, param.25, X = data, arr.dist = 'lognormal', arr.str = str.25)
res.25 <- unpack_param(opt.25$estimate, 'lognormal', str.25, min.age, n.cohorts, K)
AIC.25 <- 2*opt.25$minimum + 2*length(param.25)


### model 26
### two log-normal distributions, mixture proportion cohort dependent
param.26 <- c(rep(0, 2), rep(0, 2), rep(0, n.cohorts), 0, 0)
param.26 <- opt.26$estimate  # max iterations exceeded, restart at current values
str.26 <- list(c('shared', 'shared'), c('shared', 'shared'), 'cohort')
opt.26 <- nlm(likelihood_cohort, param.26, X = data, arr.dist = 'lognormal', arr.str = str.26)
res.26 <- unpack_param(opt.26$estimate, 'lognormal', str.26, min.age, n.cohorts, K)
AIC.26 <- 2*opt.26$minimum + 2*length(param.26)
CIs.26 <- bootstrap_fn(999, param.26, data, 'lognormal', str.26)


### model 27
### two log-normal distributions, both means cohort dependent
param.27 <- c(rep(0, 2*n.cohorts), rep(0, 2), 0, 0, 0)
str.27 <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'shared')
opt.27 <- nlm(likelihood_cohort, param.27, X = data, arr.dist = 'lognormal', arr.str = str.27)
res.27 <- unpack_param(opt.27$estimate, 'lognormal', str.27, min.age, n.cohorts, K)
AIC.27 <- 2*opt.27$minimum + 2*length(param.27)


### model 28
### two log-normal distributions, one distribution (mean and sd) cohort dependent
param.28 <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0, 0, 0)
param.28 <- opt.28$estimate  # max iterations exceeded, restart at current values
str.28 <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'shared')
opt.28 <- nlm(likelihood_cohort, param.28, X = data, arr.dist = 'lognormal', arr.str = str.28)
res.28 <- unpack_param(opt.28$estimate, 'lognormal', str.28, min.age, n.cohorts, K)
AIC.28 <- 2*opt.28$minimum + 2*length(param.28)


### model 29
### two log-normal distributions, one mean and mixture proportions cohort dependent
param.29 <- c(rep(0, n.cohorts), 0, rep(0, 2), rep(0, n.cohorts), 0, 0)
param.29 <- opt.29$estimate  # max iterations exceeded, restart at current values
str.29 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'cohort')
opt.29 <- nlm(likelihood_cohort, param.29, X = data, arr.dist = 'lognormal', arr.str = str.29)
res.29 <- unpack_param(opt.29$estimate, 'lognormal', str.29, min.age, n.cohorts, K)
AIC.29 <- 2*opt.29$minimum + 2*length(param.29)


### model 30
### two log-normal distributions, both sds cohort dependent
param.30 <- c(rep(0, 2), rep(0, 2*n.cohorts), 0, 0, 0)
param.30 <- opt.30$estimate  # max iterations exceeded, restart at current values
str.30 <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'shared')
opt.30 <- nlm(likelihood_cohort, param.30, X = data, arr.dist = 'lognormal', arr.str = str.30)
res.30 <- unpack_param(opt.30$estimate, 'lognormal', str.30, min.age, n.cohorts, K)
AIC.30 <- 2*opt.30$minimum + 2*length(param.30)


### model 31
### two log-normal distributions, one sd and mixture proportions cohort dependent
param.31 <- c(rep(0, 2), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
str.31 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'cohort')
opt.31 <- nlm(likelihood_cohort, param.31, X = data, arr.dist = 'lognormal', arr.str = str.31)
res.31 <- unpack_param(opt.31$estimate, 'lognormal', str.31, min.age, n.cohorts, K)
AIC.31 <- 2*opt.31$minimum + 2*length(param.31)
