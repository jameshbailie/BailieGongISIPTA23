#Some code looking at the two definitions of convexity

library("tidyverse")

dPMM <- function(P,Q) {
  (Q-P)/(1-Q)
}
dLV <- function(P,Q) {
  (Q-P)/Q
}

P <- 0.5
QOne <- 0.5
QTwo <- 0.7


convexCombination <- function(alpha, distance) {
  do.call(distance, list(P, QOne))*alpha + do.call(distance, list(P, QTwo))*(1-alpha)
}

mixtureDistance <- function(alpha, distance) {
  Q <- alpha*QOne + (1-alpha)*QTwo
  do.call(distance, list(P,Q))
}

makePlot <- function(distance) {
  alphas <- seq(from = 0, to = 1, by = 0.01)
  convexValues <- convexCombination(alphas, distance)
  mixtureValues <- mixtureDistance(alphas, distance)
  
  tibble(alphas = alphas,
         convexValues = convexValues,
         mixtureValues = mixtureValues) %>%
    pivot_longer(cols = !alphas) %>% 
    ggplot() +
    geom_line(aes(x = alphas, y = value, colour = name))
}

makePlot(dPMM)
makePlot(dLV)
