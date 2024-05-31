library(tidyverse)
library(ggplot2)
library(jmuOutlier) # for laplace
library(pracma) # for incomplete gamma function
library(Matching)
library(reshape2)
library(rstan)
library(RColorBrewer)
library(xtable)
library(tikzDevice)
library(graphics)
library(cowplot)
options(scipen = 10)

mypalette <- brewer.pal(3, 'Set1')

common_theme <- function(...) {
  theme_bw(...)
}
title_theme <- calc_element("axis.text.x", common_theme())


############################################################################
#                          laplace - Figure 1.                             #
############################################################################
# Figure 1 in Bailie and Gong 2024+

#For adding the curly braces:
library(grid)
library(pBrackets) 


bracketsGrob <- function(...){
  l <- list(...)
  e <- new.env()
  e$l <- l
  grid:::recordGrob(  {
    do.call(grid.brackets, l)
  }, e)
}

one_density <- function(x, eps, mu){
  dlaplace(x = x, mu, sqrt(2)/eps)
}

second_density <- function(x, eps, mu){
  dlaplace(x = x, mu+1, sqrt(2)/eps)
}

#Making the plot data:
eps <- 2
mu <- 10
t <- seq(8,13,0.01)
firstDensity <- one_density(x = t, eps = eps, mu = mu)
secondDensity <- second_density(x = t, eps = eps, mu = mu)
numPoints <- length(t)

plotData <- tibble(t = rep(t, 2), px = c(firstDensity, secondDensity),
                   type = c(rep(0,numPoints), rep(1, numPoints)))


#Variables for formatting the plot:
col1 <- mypalette[2] 
col2 <- mypalette[2] 
bracketCol <- "grey30"
vlineCol <- bracketCol
b1 <- bracketsGrob(5/11+0.004, 0.05, 0.591*10/11+0.004, 0.05, h=0.05, lwd=1, col=bracketCol)

#Making the plot
p <- ggplot(plotData, aes(x = t, y = px, linetype=factor(type)))

# Add the lines
p <- p + geom_line(colour = col1) +
  scale_linetype_manual(values=c("solid", "dashed"), 
                        labels = c("$p_x(t)$",
                                   "$p_{x'}(t)$"))

# p <- p + scale_x_continuous(breaks = 5:16,
#                             labels = c(rep("", 4),"$q(x)-\\Delta(q)$", "$q(x)$", "$q(x')$", rep("", 5)))

#Add the ratio of densities
p <- p + geom_segment(x = 8.5975, xend = 8.5975, y = 0, yend = one_density(x = 8.5975, eps = eps, mu = mu),#-0.001,
                      linetype = "solid", color = "grey30", linewidth = 0.4) +
  geom_segment(x = 8.6975, xend = 8.6975, y = 0, yend = second_density(x = 8.6975, eps = eps, mu = mu),#-0.001,
               linetype = "solid", color = "grey30", linewidth = 0.4)

p <- p + draw_label(x = 8.65, y = 0.1,#one_density(x = 9, eps = eps, mu = mu) + 0.02, 
                   label = "$.05\\exp(\\varepsilon)$", 
                   color = "grey30",
                   hjust = 0, 
                   fontfamily = title_theme$family,
                   fontface = title_theme$face,
                   size = title_theme$size) +
  draw_label(x = 8.75, y = 0.03,#one_density(x = 9, eps = eps, mu = mu) + 0.02,
            label = list("$.05$"), 
            color = "grey30",
            hjust = 0,
            fontfamily = title_theme$family,
            fontface = title_theme$face,
            size = title_theme$size)

# Add the curly bracket
p <- p + annotation_custom(b1)

# Add the label "\Delta(q)=1"
p <- p + draw_label(x = 10.5, y = 0.03,
                   label = "$\\Delta(q)=1$",
                   fontfamily = title_theme$family,
                   fontface = title_theme$face,
                   size = title_theme$size,
                   hjust = 0.5, vjust = -0.5,
                   color = bracketCol)

# Move the legend inside the plot
p <- p + common_theme() + 
  labs(y = 'Density', x = '$t$', linetype = NULL) + 
  theme(legend.position = c(.8,.775),
        legend.text.align = 0,
        legend.box.background = element_rect(colour = "grey30"),
        legend.margin=margin(t=-0.19, r=0.1, b=0, l=0, unit="cm"))


print(p)
tikz('bailie24Figs/laplaceTikZ.tex', width = 5, height = 3.5)
p
dev.off()



############################################################################

############################################################################
#                       marginal1 & 2 - Figure 2                           #
############################################################################

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
    labs(title = 'density bounds for p( t | \u03B8 )',  x = 't', y = 'density',
         caption = paste0('n = ', nobs, ', \u03B5 = ', epsilon))
  
  
  p_both <- p_ + geom_function(fun = upper_ml, args = list(eps = epsilon, n = nobs), colour = color) #+ labs(subtitle = 'Upper and lower densities')
  
  p_lower <- p_ + labs(subtitle = 'Lower density only')
  
  return(list(both=p_both, lower=p_lower))
}

plot_ml(nobs = 10, epsilon = 0.1, color = mypalette[2])$both
ggsave('bailie24Figs/marginal1.pdf', width = 4, height = 3.5, device = cairo_pdf)

plot_ml(nobs = 10, epsilon = 0.25, color = mypalette[2])$both
ggsave('bailie24Figs/marginal2.pdf', width = 4, height = 3.5, device = cairo_pdf)

############################################################################

############################################################################
#                        marginal_rr - Figure 3                            #
############################################################################

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
  theme_bw() + scale_color_manual(values = mypalette[c(2, 2)], guide = 'none') +
  scale_x_continuous(breaks = seq(2,10,2)) +
  labs(title = 'density bounds for p( t | \u03B8 )',  x = '|t|', y = 'density',
       caption = paste0('max |t| = ', nt, ', \u03B5 = ', epsilon))


p_rr
ggsave('bailie24Figs/marginal_rr.pdf', width = 5, height = 3.5, device = cairo_pdf)


############################################################################

############################################################################
#                         posterior - Figure 4                             #
############################################################################


clamp <- function(x, lower, upper){
  y = x
  y[y < lower] <- lower
  y[y > upper] <- upper
  return(y)
}

post_abc_clamp <- function(y, nsample = 1e3, a = 1, b = 1, 
                           eps = 1, a0 = 0, a1 = 100){
  # Gamma-Pois DP with clamping, ABC algorithm
  # lower, upper: clamping range
  sens <- a1 - a0
  s_ <- sqrt(2)*sens/eps
  temp_ <- function(nsim){
    l_ <- rgamma(nsim, shape = a, rate = b)
    x_ <- rpois(nsim, l_) %>% clamp(lower = a0, upper = a1)
    probs <- dlaplace(y-x_, mean = 0, sd = s_)/dlaplace(0, mean = 0, sd = s_)
    ind <- which(runif(nsim) < probs)
    return(list(lambda = l_[ind], yield = max(1, length(ind))/nsim))
  }
  tp_ <- temp_(nsample)
  lambda <- tp_$lambda; yield = tp_$yield
  cat('Yield is about', yield)
  if (length(lambda) < nsample){
    lambda <- c(lambda, temp_(ceil(nsample/yield))$lambda)
  }
  return(list(lambda, yield))
}


prior_predictive <- function(nsim, a, b, a0, a1, eps){
  l_ <- rgamma(nsim, shape = a, rate = b)
  x_ <- rpois(nsim, l_) %>% clamp(lower = a0, upper = a1)
  t_ <- x_ + rlaplace(nsim, mean = 0, sd = sqrt(2)*(a1-a0)/eps)
  return(list(x = x_, t = t_))
}

alpha <- 3; beta <- 1; # hyperparameters
a0 <- 0; a1 <- 6
epsilon <- 1

set.seed(10)
nruns = 10
pp0 <- prior_predictive(nruns, a = alpha, b = beta, a0 = a0, a1 = a1, eps = epsilon)

set.seed(100)
draws = list()
for (i in 1:nruns){
  run_ <- post_abc_clamp(y = pp0$t[i], nsample = 2e3, 
                         a = alpha, b = beta, eps = epsilon, a0 = a0, a1 = a1)
  draws[[i]] <- run_[[1]]
}


draws.dat <- data.frame(sapply(draws, function(x) x[1:max(lengths(draws))])) %>%
  pivot_longer(everything())

ll <- seq(range(draws.dat$value, na.rm = T)[1], range(draws.dat$value, na.rm = T)[2], by = 0.05)

dat <- data.frame(lambda = ll) %>%
  mutate(prior = dgamma(lambda, shape = alpha, rate = beta)) %>%
  mutate(upper = prior*exp(epsilon)) %>%
  mutate(lower = prior*exp(-epsilon)) %>%
  pivot_longer(!lambda)
#head(dat)

p <- ggplot(dat, aes(x = lambda, y = value, linetype = name)) +
  geom_density(data = draws.dat, aes(group = name, x = value), color = 'gray80', inherit.aes = F) +
  geom_line(color = mypalette[2]) +
  scale_linetype_manual(values = c('solid', 'dashed', 'solid'), guide = 'none') +
  theme_bw() +
  labs(x = '\u03B8', y = 'posterior density', title = 'density bounds for posterior p( \u03B8 | t )',
       caption = paste0('prior for \u03B8 is Gamma(', alpha, ',', beta, '), clamping range is [', a0, ',', a1, '], \u03B5 = ', epsilon))

p
ggsave('bailie24Figs/posterior.pdf', width = 5, height = 3.5, device = cairo_pdf)

############################################################################

