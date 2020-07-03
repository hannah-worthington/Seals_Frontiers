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

# storage for AIC values
AIC <- data.frame('arrivals model' = numeric(0), 'AIC' = numeric (0))



### One normal distribution for recruitment, constant p, constant phi

### model 1
### one normal distribution, mean and sd shared across cohorts
param.1 <- c(0, 0, 0, 0)
str.1 <- list('shared', 'shared')
opt.1 <- nlm(likelihood_cohort, param.1, X = data, arr.dist = normal_arr, arr.str = str.1)
res.1 <- normal_arr(opt.1$estimate, str.1, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared mean and sd', 'AIC' = 2*opt.1$minimum + 2*length(param.1))) 


### model 2
### one normal distribution, mean cohort dependent, sd shared across cohorts
param.2 <- c(rep(0, n.cohorts), 0, 0, 0)
str.2 <- list('cohort', 'shared')
opt.2 <- nlm(likelihood_cohort, param.2, X = data, arr.dist = normal_arr, arr.str = str.2)
res.2 <- normal_arr(opt.2$estimate, str.2, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared sd', 'AIC' = 2*opt.2$minimum + 2*length(param.2)))


### model 3
### one normal distribution, mean shared across cohorts, sd cohort dependent
param.3 <- c(0, rep(0, n.cohorts), 0, 0)
str.3 <- list('shared', 'cohort')
opt.3 <- nlm(likelihood_cohort, param.3, X = data, arr.dist = normal_arr, arr.str = str.3)
res.3 <- normal_arr(opt.3$estimate, str.3, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared mean', 'AIC' = 2*opt.3$minimum + 2*length(param.3)))


### model 4
### one normal distribution, mean and sd cohort dependent
param.4 <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, 0)
param.4 <- opt.4$estimate # max iteration exceeded, restart at current values
str.4 <- list('cohort', 'cohort')
opt.4 <- nlm(likelihood_cohort, param.4, X = data, arr.dist = normal_arr, arr.str = str.4)
res.4 <- normal_arr(opt.4$estimate, c('cohort', 'cohort'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal cohort mean and sd', 'AIC' = 2*opt.4$minimum + 2*length(param.4)))



### Mixture of two normals for recruitment, constant p, constant phi

### model 5
### two normal distributions, mean, sd and mixture proportion all shared
param.5 <- c(rep(0, 2), rep(0, 2), 0, 0, 0)
str.5 <- list(c('shared', 'shared'), c('shared', 'shared'), 'shared')
opt.5 <- nlm(likelihood_cohort, param.5, X = data, arr.dist = normal_arr, arr.str = str.5)
res.5 <- normal_arr(opt.5$estimate, str.5, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared mean, sd and mixture proportion', 'AIC' = 2*opt.5$minimum + 2*length(param.5)))


### model 6
### two normal distributions, one mean cohort dependent
param.6 <- c(rep(0, n.cohorts), 0, rep(0, 2), 0, 0, 0)
str.6 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'shared')
opt.6 <- nlm(likelihood_cohort, param.6, X = data, arr.dist = normal_arr, arr.str = str.6)
res.6 <- normal_arr(opt.6$estimate, str.6, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared mean, sd and mixture proportion', 'AIC' = 2*opt.6$minimum + 2*length(param.6)))


### model 7
### two normal distributions, one sd cohort dependent
param.7 <- c(rep(0, 2), rep(0, n.cohorts), 0, 0, 0, 0)
param.7 <- opt.7$estimate # max iterations exceeded, restart at current values
str.7 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'shared')
opt.7 <- nlm(likelihood_cohort, param.7, X = data, arr.dist = normal_arr, arr.str = str.7)
res.7 <- normal_arr(opt.7$estimate, str.7, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared means, one sd and mixture proportion', 'AIC' = 2*opt.7$minimum + 2*length(param.7)))


### model 8
### two normal distributions, mixing proportion cohort dependent
param.8 <- c(rep(0, 2), rep(0, 2), rep(0, n.cohorts), 0, 0)
param.8 <- opt.8$estimate # max iterations exceeded, restart at current values
str.8 <- list(c('shared', 'shared'), c('shared', 'shared'), 'cohort')
opt.8 <- nlm(likelihood_cohort, param.8, X = data, arr.dist = normal_arr, arr.str = str.8)
res.8 <- normal_arr(opt.8$estimate, str.8, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared means and sds', 'AIC' = 2*opt.8$minimum + 2*length(param.8)))
