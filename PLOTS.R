# libraries
library(ggplot2)

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



### produce histograms of observed age at first recapture

# find age at first breeding observation
age.first <- list()
for (c in 1:n.cohorts)  {
  age <- rep(0, n[c])
  for (i in 1:n[c])  {
    age[i] <- which(data[[c]][i,] == 1)[1]
  }
  age.first[[c]] <- age
}
age.first <- unlist(age.first)

# create data.frame for plotting
age.plot <- data.frame('age' = age.first, 'cohort' = as.factor(rep(1991:1994, times = n)))

# plot histograms of age of first recapture
ggplot(data = age.plot) +
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
  



### Produce plots from CohortsStopover_RUN model results

# model number from CohortsStopover_RUN
model <- 5

# create data for plotting
res <- eval(as.name(paste('res', model, sep='.')))
res.plot <- data.frame('beta' = unlist(res$beta), 'cohort' = as.factor(rep(1991:1994, times = K)), 'age' = sequence(K))

# plot arrival distributions with observed ages underneath
ggplot() +
  # add the observed ages
  geom_histogram(data = age.plot,
                 mapping = aes(x = age, 
                               y = ..density..,
                               fill = cohort),
                 colour = 'black',
                 binwidth = 1,
                 alpha = 0.5) +
  geom_line(data = res.plot, 
            mapping = aes(x = age, 
                          y = beta,
                          col = cohort),
            size = 2)  +
  facet_wrap( ~ cohort)  +
  scale_color_brewer(palette = 'Set1') +
  scale_fill_brewer(palette = 'Set1') +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold',
                                  size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
  
