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
opt.1 <- nlm(likelihood_cohort, param.1, X = data, arr.dist = normal_one, arr.str = c('shared', 'shared'))
res.1 <- normal_one(opt.1$estimate, c('shared', 'shared'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared mean and sd', 'AIC' = 2*opt.1$minimum + 2*length(param.1))) 


### model 2
### one normal distribution, mean cohort dependent, sd shared across cohorts
param.2 <- c(rep(0, n.cohorts), 0, 0, 0)
opt.2 <- nlm(likelihood_cohort, param.2, X = data, arr.dist = normal_one, arr.str = c('cohort', 'shared'))
res.2 <- normal_one(opt.2$estimate, c('cohort', 'shared'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared sd', 'AIC' = 2*opt.2$minimum + 2*length(param.2)))


### model 3
### one normal distribution, mean shared across cohorts, sd cohort dependent
param.3 <- c(0, rep(0, n.cohorts), 0, 0)
opt.3 <- nlm(likelihood_cohort, param.3, X = data, arr.dist = normal_one, arr.str = c('shared', 'cohort'))
res.3 <- normal_one(opt.3$estimate, c('shared', 'cohort'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal shared mean', 'AIC' = 2*opt.3$minimum + 2*length(param.3)))


### model 4
### one normal distribution, mean and sd cohort dependent
param.4 <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, 0)
opt.4 <- nlm(likelihood_cohort, param.4, X = data, arr.dist = normal_one, arr.str = c('cohort', 'cohort'))
res.4 <- normal_one(opt.4$estimate, c('cohort', 'cohort'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal cohort mean and sd', 'AIC' = 2*opt.4$minimum + 2*length(param.4)))



### Mixture of two normals for recruitment, constant p, constant phi

### model 5
### two normal distributions, mean, sd and mixture proportion all shared
param.4 <- c(rep(0, n.cohorts), rep(0, n.cohorts), 0, 0)
opt.4 <- nlm(likelihood_cohort, param.4, X = data, arr.dist = normal_one, arr.str = c('cohort', 'cohort'))
res.4 <- normal_one(opt.4$estimate, c('cohort', 'cohort'), min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'one normal cohort mean and sd', 'AIC' = 2*opt.4$minimum + 2*length(param.4)))