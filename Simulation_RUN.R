# libraries
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# source functions
source('Simulation_FUNCTIONS.R')

# simulation 1
n <- 20
C <- 4
K <- 20
logmeans <- c(2, 2.5)
logsds <- c(0.25, 0.05)
means <- exp(logmeans + (logsds^2)/2)
sds <- exp(2*logmeans + logsds^2) * (exp(logsds^2) - 1)
w <- c(1, 0.7, 0.3, 0)
r <- list('means' = matrix(logmeans, nrow = 2, ncol = C), 'sds' = matrix(logsds, nrow = 2, ncol = C), 'w' = w)
p <- 0.7
phi <- 0.8
arr.dist <- 'lognormal'
n.sim <- 1000

sim.res.1 <- sim_study(n, C, K, r, p, phi, arr.dist, n.sim)


# simulation 2
n <- 50
sim.res.2 <- sim_study(n, C, K, r, p, phi, arr.dist, n.sim)


# simulation 3
n <- 100
sim.res.3 <- sim_study(n, C, K, r, p, phi, arr.dist, n.sim)



### Produce plots from sim.res results

# simulation number
sim.no <- 3
sim <- eval(as.name(paste('sim.res', sim.no, sep = '.')))

# create data for plotting
sim.res.plot.multiple <- data.frame('beta' = unlist(sim$multiple$beta),
                                    'cohort' = as.factor(rep(0:(C-1), times = n.sim*(K:(K-C+1)))),
                                    'age' = rep(sequence(K:(K-C+1)), each = n.sim),
                                    'simulation' = rep(1:n.sim, times = sum(K:(K-C+1))),
                                    'model' = rep('multiple', times = n.sim*sum(K:(K-C+1))))
sim.res.plot.single <- data.frame('beta' = unlist(sim$single$beta),
                                  'cohort' = as.factor(rep(0:(C-1), times = n.sim*(K:(K-C+1)))),
                                  'age' = rep(sequence(K:(K-C+1)), each = n.sim),
                                  'simulation' = rep(1:n.sim, times = sum(K:(K-C+1))),
                                  'model' = rep('single', times = n.sim*sum(K:(K-C+1))))
sim.res.plot <- rbind(sim.res.plot.multiple, sim.res.plot.single)

sim.p.phi.plot <- data.frame('p' = c(sim$multiple$p, as.vector(sim$single$p)),
                         'phi' = c(sim$multiple$phi, as.vector(sim$single$phi)),
                         'cohort' = as.factor(c(rep('all', n.sim), rep(0:(C-1), times = n.sim))),
                         'model' = c(rep('multiple', n.sim), rep('single', n.sim*C)))

# create truth for plotting
beta.true <- list()
for (c in 1:C)  {
  beta.true[[c]] <- twolognormalbetas(logmeans, logsds, w[c], K = K-c+1, min.age = 3)
}
res.plot.true <- data.frame('beta' = unlist(beta.true),
                            'cohort' = as.factor(rep(0:(C-1), times = K:(K-C+1))),
                            'age' = sequence(K:(K-C+1)))

# define colour palette for plotting to match across cohorts
my.palette <- brewer.pal(5, 'Set2')


# plot simulation results with truth overlaid
beta.plot <- ggplot() +
  geom_line(data = sim.res.plot,
            mapping = aes(x = age,
                          y = beta,
                          col = cohort,
                          group = simulation),
            alpha = 0.7,
            size = 0.5) +
  geom_line(data = res.plot.true,
            mapping = aes(x = age,
                          y = beta),
            colour = 'black',
            size = 1.5) +
  facet_wrap( ~ cohort + model, nrow = 2) +
  scale_color_manual(values = my.palette[1:4]) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold',
                                  size = 12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

p.plot <- ggplot() +
  geom_boxplot(data = sim.p.phi.plot,
               mapping = aes(x = p,
                             fill = cohort)) +
  geom_vline(xintercept = p,
             size = 1.5) +
  facet_wrap( ~ model) +
  scale_fill_manual(values = my.palette[c(1, 2, 3, 4, 5)]) +
  theme_minimal() +
  theme(strip.text = element_text(face = 'bold',
                                  size = 12),
        axis.title = element_text(size = 12),
        axis.text.y = element_blank())

phi.plot <- ggplot() +
  geom_boxplot(data = sim.p.phi.plot,
               mapping = aes(x = phi,
                             fill = cohort)) +
  geom_vline(xintercept = phi,
             size = 1.5) +
  facet_wrap( ~ model) +
  scale_fill_manual(values = my.palette[c(1, 2, 3, 4, 5)]) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(face = 'bold',
                                  size = 12),
        axis.title = element_text(size = 12),
        axis.text.y = element_blank())

beta.plot + {
  p.plot + phi.plot
} + plot_layout(ncol = 1, heights = c(2, 1))
