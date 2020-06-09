# libraries
library(ggplot2)



### produce histograms of observed age at first recapture

# constants
n_cohorts <- 4

# read in data
data <- list()
for (c in 1:n_cohorts)  {
  data[[c]] <- read.csv(paste(getwd(), '/Data/199', c, 'mark.csv', sep=''), 
                        header = T)
}

# number of individuals and capture occasions in each cohort
n <- rep(0, n_cohorts)
K <- rep(0, n_cohorts)
for (c in 1:n_cohorts)  {
  n[c] <- length(data[[c]][,1])
  K[c] <- length(data[[c]][1,])
}

# find age at first breeding observation
age_first <- list()
for (c in 1:n_cohorts)  {
  age <- rep(0, n[c])
  for (i in 1:n[c])  {
    age[i] <- which(data[[c]][i,] == 1)[1]
  }
  age_first[[c]] <- age
}
age_first <- unlist(age_first)

# create data.frame for plotting
age_plot <- data.frame('age' = age_first, 'cohort' = as.factor(rep(1991:1994, times = n)))

# plot histograms of age of first recapture
ggplot(data = age_plot) +
  geom_histogram(mapping = aes(x = age, 
                               fill = cohort),
                 colour = 'black',
                 binwidth = 1) +
  facet_wrap( ~ cohort) +
  geom_vline(mapping = aes(xintercept = 3)) +
  annotate(geom = 'text',
           x = 3,
           y = 5.5,
           label = 'age 3') +
  geom_vline(mapping = aes(xintercept = 12)) +
  annotate(geom = 'text',
           x = 12,
           y = 5.5,
           label = 'age 12') +
  scale_fill_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold',
                                  size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
  

