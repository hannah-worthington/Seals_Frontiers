# source functions
source('CohortsStopover_FUNCTIONS.R')

# constants
n.cohorts <- 4

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
opt.1 <- nlm(likelihood_cohort, param.1, X = data, arr.dist = normal_one, arr.str = c('shared', 'shared'))
res.1 <- normal_one(opt.1$estimate, c('shared', 'shared'), 3, 4, K)





param.test <- c(21.13785, 104.39981, -196.93070, 402.59331)
test <- normal_one(param.1, c('shared', 'shared'), 3, 4, K)
test2 <- onenormalbetas(exp(param.1[1]), exp(param.1[2]), K[1], 3)
test3 <- HMM.str(test, n.cohorts, K)
test4 <- likelihood_cohort(param.1, data, arr.dist = normal_one, arr.str = c('shared', 'shared'))