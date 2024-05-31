library(tidyverse)
setwd("C:/Users/User/Documents/Overleaf/IP and DP/code")

n <- 100
epsilon <- 0.01
t <- seq(from = -n, to = 2*n, by = 0.1)


upperBound <- epsilon/2*exp(epsilon*(n/2-abs(t-n/2))) 
lowerBound <- (t < 0)*epsilon/2*exp(-epsilon*(n-t)) +
  (t >= 0 & t <= n)*(epsilon/2*exp(-epsilon*n)) +
  (t > n)*epsilon/2*exp(-epsilon*t)

plotData <- tibble(t = t, upperBound = upperBound, lowerBound = lowerBound) |>
  pivot_longer(cols = c(upperBound, lowerBound), names_to = "Bound", values_to = "y")

ggplot(plotData, aes(x = t, y = y, colour = Bound)) +
  geom_line() +
  ylab(expression("Marginal Likelihood p(t| " * theta * ")")) +
  scale_x_continuous(labels = c("-n", "0", "n", "2n"))

ggsave("../figs/LaplaceMarginalLikelihood.png")
       