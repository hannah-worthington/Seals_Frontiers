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


### model 9
### two normal distributions, both means cohort dependent
param.9 <- c(rep(0, 2*n.cohorts), rep(0, 2), 0, 0, 0)
param.9 <- opt.9$estimate  # max iterations exceeded, restart at current values
str.9 <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'shared')
opt.9 <- nlm(likelihood_cohort, param.9, X = data, arr.dist = normal_arr, arr.str = str.9)
res.9 <- normal_arr(opt.9$estimate, str.9, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared sds and mixture proportion', 'AIC' = 2*opt.9$minimum + 2*length(param.9)))


### model 10
### two normal distributions, one normal (mean and sd) cohort dependent
param.10 <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0, 0, 0)
str.10 <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'shared')
opt.10 <- nlm(likelihood_cohort, param.10, X = data, arr.dist = normal_arr, arr.str = str.10)
res.10 <- normal_arr(opt.10$estimate, str.10, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared normal and mixture proportion', 'AIC' = 2*opt.10$minimum + 2*length(param.10)))


### model 11
### two normal distributions, one mean and mixture proportions cohort dependent
param.11 <- c(rep(0, n.cohorts), 0, rep(0, 2), rep(0, n.cohorts), 0, 0)
param.11 <- opt.11$estimate  # max iterations exceeded, restart at current values
str.11 <- list(c('cohort', 'shared'), c('shared', 'shared'), 'cohort')
opt.11 <- nlm(likelihood_cohort, param.11, X = data, arr.dist = normal_arr, arr.str = str.11)
res.11 <- normal_arr(opt.11$estimate, str.11, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared mean and both sds', 'AIC' = 2*opt.11$minimum + 2*length(param.11)))


### model 12
### two normal distributions, both sds cohort dependent
param.12 <- c(rep(0, 2), rep(0, 2*n.cohorts), 0, 0, 0)
param.12 <- opt.12$estimate  # max iterations exceeded, restart at current value
str.12 <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'shared')
opt.12 <- nlm(likelihood_cohort, param.12, X = data, arr.dist = normal_arr, arr.str = str.12)
res.12 <- normal_arr(opt.12$estimate, str.12, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared means and mixture proportion', 'AIC' = 2*opt.12$minimum + 2*length(param.12)))


### model 13
### two normal distributions, one sd and mixture proportions cohort dependent
param.13 <- c(rep(0, 2), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
param.13 <- opt.13$estimate  # max iterations exceeded, restart at current values
str.13 <- list(c('shared', 'shared'), c('cohort', 'shared'), 'cohort')
opt.13 <- nlm(likelihood_cohort, param.13, X = data, arr.dist = normal_arr, arr.str = str.13)
res.13 <- normal_arr(opt.13$estimate, str.13, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared means and one sd', 'AIC' = 2*opt.13$minimum + 2*length(param.13)))


### model 14
### two normal distributions, both means and one sd cohort dependent
param.14 <- c(rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0, 0, 0)
str.14 <- list(c('cohort', 'cohort'), c('cohort', 'shared'), 'shared')
opt.14 <- nlm(likelihood_cohort, param.14, X = data, arr.dist = normal_arr, arr.str = str.14)
res.14 <- normal_arr(opt.14$estimate, str.14, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared sd and mixture proportion', 'AIC' = 2*opt.14$minimum + 2*length(param.14)))


### model 15
### two normal distributions, both means and mixture proportions cohort dependent
param.15 <- c(rep(0, 2*n.cohorts), rep(0, 2), rep(0, n.cohorts), 0, 0)
param.15 <- opt.15$estimate  # max iterations exceeded, restart at current values
str.15 <- list(c('cohort', 'cohort'), c('shared', 'shared'), 'cohort')
opt.15 <- nlm(likelihood_cohort, param.15, X = data, arr.dist = normal_arr, arr.str = str.15)
res.15 <- normal_arr(opt.15$estimate, str.15, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared sds', 'AIC' = 2*opt.15$minimum + 2*length(param.15)))


### model 16
### two normal distributions, one mean and both sds cohort dependent
param.16 <- c(rep(0, n.cohorts), 0, rep(0, 2*n.cohorts), 0, 0, 0)
param.16 <- opt.16$estimate  # max iterations exceeded, restart at current values
str.16 <- list(c('cohort', 'shared'), c('cohort', 'cohort'), 'shared')
opt.16 <- nlm(likelihood_cohort, param.16, X = data, arr.dist = normal_arr, arr.str = str.16)
res.16 <- normal_arr(opt.16$estimate, str.16, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared mean and mixture proportions', 'AIC' = 2*opt.16$minimum + 2*length(param.16)))


### model 17
### two normal distributions, one normal (mean and sd) and mixture proportions cohort dependent
param.17 <- c(rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
param.17 <- opt.17$estimate  # max iterations exceeded, restart at current values
str.17 <- list(c('cohort', 'shared'), c('cohort', 'shared'), 'cohort')
opt.17 <- nlm(likelihood_cohort, param.17, X = data, arr.dist = normal_arr, arr.str = str.17)
res.17 <- normal_arr(opt.17$estimate, str.17, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared normal', 'AIC' = 2*opt.17$minimum + 2*length(param.17)))


### model 18
### two normal distributions, both sds and mixture proportions cohort dependent
param.18 <- c(rep(0, 2), rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0)
param.18 <- opt.18$estimate  # max iterations exceeded, restart at current values
str.18 <- list(c('shared', 'shared'), c('cohort', 'cohort'), 'cohort')
opt.18 <- nlm(likelihood_cohort, param.18, X = data, arr.dist = normal_arr, arr.str = str.18)
res.18 <- normal_arr(opt.18$estimate, str.18, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared means', 'AIC' = 2*opt.18$minimum + 2*length(param.18)))


### model 19
### two normal distributions, both normals (means and sds) cohort dependent
param.19 <- c(rep(0, 2*n.cohorts), rep(0, 2*n.cohorts), 0, 0, 0)
param.19 <- opt.19$estimate  # max iterations exceeded, restart at current values
str.19 <- list(c('cohort', 'cohort'), c('cohort', 'cohort'), 'shared')
opt.19 <- nlm(likelihood_cohort, param.19, X = data, arr.dist = normal_arr, arr.str = str.19)
res.19 <- normal_arr(opt.19$estimate, str.19, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals shared mixture proportions', 'AIC' = 2*opt.19$minimum + 2*length(param.19)))


### model 20
### two normal distributions, both means, one sd and mixture proportions cohort dependent
param.20 <- c(rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, rep(0, n.cohorts), 0, 0)
param.20 <- opt.20$estimate  # max iterations exceeded, restart at current values
str.20 <- list(c('cohort', 'cohort'), c('cohort', 'shared'), 'cohort')
opt.20 <- nlm(likelihood_cohort, param.20, X = data, arr.dist = normal_arr, arr.str = str.20)
res.20 <- normal_arr(opt.20$estimate, str.20, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared sd', 'AIC' = 2*opt.20$minimum + 2*length(param.20)))


### model 21
### two normal distributions, one mean, both sds and mixture proportions cohort dependent
param.21 <- c(rep(0, n.cohorts), 0, rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0)
param.21 <- opt.21$estimate  # max iterations exceeded, restart at current values
str.21 <- list(c('cohort', 'shared'), c('cohort', 'cohort'), 'cohort')
opt.21 <- nlm(likelihood_cohort, param.21, X = data, arr.dist = normal_arr, arr.str = str.21)
res.21 <- normal_arr(opt.21$estimate, str.21, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals one shared mean', 'AIC' = 2*opt.21$minimum + 2*length(param.21)))


### model 22
### two normal distributions, both normals (means and sds) and mixture proportions cohort dependent
param.22 <- c(rep(0, 2*n.cohorts), rep(0, 2*n.cohorts), rep(0, n.cohorts), 0, 0)
param.22 <- opt.22$estimate  # max iterations exceeded, restart at current values
str.22 <- list(c('cohort', 'cohort'), c('cohort', 'cohort'), 'cohort')
opt.22 <- nlm(likelihood_cohort, param.22, X = data, arr.dist = normal_arr, arr.str = str.22)
res.22 <- normal_arr(opt.22$estimate, str.22, min.age, n.cohorts, K)
AIC <- rbind(AIC, data.frame('arrivals model' = 'two normals arrivals cohort dependent', 'AIC' = 2*opt.22$minimum + 2*length(param.22)))
