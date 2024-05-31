
library(tidyverse)
library(ggplot2)
library(jmuOutlier) # for laplace

library(ggplot2)

library(RColorBrewer)
library(xtable)
options(scipen = 10)

mypalette <- brewer.pal(3, 'Set1')
col1 <- "#c6494a"
plot.figures <- T


#Laplace density bounds

upper_ml <- function(x, eps, n){
  #eps/2*exp(eps*(n/2 - abs(x - n/2)))
  dlaplace(x = x, n/2, sqrt(2)/eps)*exp(eps*n/2)
}

lower_ml <- function(x, eps, n){
  #eps/2*exp(eps*(n/2 - abs(x - n/2)))
  pmin(dlaplace(x = x, n/2, sqrt(2)/eps),
       eps/2*exp(-eps*n/2))*exp(-eps*n/2)
}


plot_ml <- function(nobs, epsilon, color = 'blue'){
  p_ <- ggplot(data.frame(x = seq(-nobs*2, nobs*3, by = 0.1)), aes(x)) +
    geom_function(fun = lower_ml, args = list(eps = epsilon, n = nobs),
                  colour = color) +
    theme_bw() +
    labs(x = 't', y = 'density')#,title = 'density bounds for p( t | \u03B8 )',  
         #caption = paste0('n = ', nobs, ', \u03B5 = ', epsilon))
  
  
  p_both <- p_ + geom_function(fun = upper_ml, args = list(eps = epsilon, n = nobs), colour = color) #+ labs(subtitle = 'Upper and lower densities')
  
  p_lower <- p_ + labs(subtitle = 'Lower density only')
  
  return(list(both=p_both, lower=p_lower))
}


plot_ml(nobs = 10, epsilon = 0.1, color = col1)$both
if (plot.figures){
  ggsave('posterBailie23Figs/marginal1.pdf', width = 5, height = 5/0.99*0.47*6/7*1.22, device = cairo_pdf)
}

plot_ml(nobs = 10, epsilon = 0.25, color = col1)$both
if (plot.figures){
  ggsave('posterBailie23Figs/marginal2.pdf', width = 5, height = 5/0.99*0.47*6/7*1.22, device = cairo_pdf)
}








# Local DP randomized response

lower_rr <- function(t, eps){
  (exp(eps) + 1)^(-t)
}

upper_rr <- function(t, eps){
  (exp(eps) + 1)^(-t)*exp(t*eps)
}


epsilon <- 1
nt = 10
dat2 <- data.frame(t = seq(1, nt, by = 1)) %>%
  mutate(lower = lower_rr(t=t, eps = epsilon)) %>%
  mutate(upper = upper_rr(t=t, eps = epsilon)) %>%
  pivot_longer(!t)
#head(dat2)

p_rr <- ggplot(dat2, aes(x = t, y = value, color = name)) +
  geom_line() + geom_point() +
  theme_bw() + scale_color_manual(values = c(col1, col1), guide = 'none') +
  scale_x_continuous(breaks = seq(2,10,2)) +
  labs( x = '|t|', y = 'density')#, title = 'density bounds for p( t | \u03B8 )',
       #caption = paste0('max |t| = ', nt, ', \u03B5 = ', epsilon))


p_rr
if (plot.figures){
  ggsave('posterBailie23Figs/marginal_rr.pdf', width = 5, height = 5/0.99*0.47*6/7*1.5, device = cairo_pdf)
}
