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
library(magrittr)
options(scipen = 10)

#setwd("~/Documents/git/BailieGongISIPTA23")
library("here")
setwd(here())


mypalette <- brewer.pal(3, 'Set1')

common_theme <- function(...) {
  theme_bw(...)
}

panelBorderColour <- calc_element('panel.border', common_theme())$colour
panelBorderWidth <- 0.8#calc_element('panel.border', common_theme())$linewidth

common_theme <- function(...) {
  theme_bw(...) +
  theme(panel.border = element_rect(linewidth = panelBorderWidth))
}
title_theme <- calc_element("axis.text.x", common_theme())

marginal12RelSize <- 0.49


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
t <- seq(5,16,0.01)
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


xAxisLimits <- c(8,13)
xAxisBreaks <- xAxisLimits[1]:xAxisLimits[2]
expandConstant <- 0.06*(xAxisLimits[2]-xAxisLimits[1])
xAxisLimits <- c(xAxisLimits[1]-expandConstant,xAxisLimits[2]+expandConstant)
xAxisLabels <- c("$q(x)-2$", "$q(x)-1$", "$q(x)$", "$q(x')$", "$q(x')+1$", "$q(x')+2$")

b1Start <- (10-xAxisLimits[1])/(xAxisLimits[2]-xAxisLimits[1])
b1End <- (11-xAxisLimits[1])/(xAxisLimits[2]-xAxisLimits[1])
b1 <- bracketsGrob(b1Start, 0.05, b1End, 0.05, h=0.05, lwd=1, col=bracketCol)

firstSegmentX <- 9.49
secondSegmentX <- 9.51
distXLabel <- 0.03
distYLabel <- -0.01

#Making the plot
p <- ggplot(plotData, aes(x = t, y = px, linetype=factor(type)))

# Add the lines
p <- p + geom_line(colour = col1) +
  scale_linetype_manual(values=c("solid", "dashed"), 
                        labels = c("$p_x(t)$",
                                   "$p_{x'}(t)$"))

p <- p + scale_x_continuous(limits = xAxisLimits,
                            breaks = xAxisBreaks,
                            expand = c(0,0),
                            labels = xAxisLabels)

#Add the ratio of densities
p <- p + geom_segment(x = firstSegmentX, xend = firstSegmentX, y = 0, yend = one_density(x = firstSegmentX, eps = eps, mu = mu),#-0.001,
                      linetype = "solid", color = panelBorderColour,#"grey30"
                      linewidth = panelBorderWidth*0.352777778#0.4 #Change the units by 0.35 because https://stackoverflow.com/questions/17311917/ggplot2-the-unit-of-size
                      ) +
  geom_segment(x = secondSegmentX, xend = secondSegmentX, y = 0, yend = second_density(x = secondSegmentX, eps = eps, mu = mu),#-0.001,
               linetype = "solid", color = panelBorderColour,#"grey30"
               linewidth = panelBorderWidth*0.352777778#0.4
  )

p <- p + draw_label(x = firstSegmentX+distXLabel, y = one_density(x = firstSegmentX, eps = eps, mu = mu)+distYLabel,#one_density(x = 9, eps = eps, mu = mu) + 0.02, 
                   label = "$.05\\exp(\\varepsilon)$", 
                   color = "grey30",
                   hjust = 0, 
                   fontfamily = title_theme$family,
                   fontface = title_theme$face,
                   size = title_theme$size) +
  draw_label(x = secondSegmentX+distXLabel, y = second_density(x = secondSegmentX, eps = eps, mu = mu)+distYLabel,#one_density(x = 9, eps = eps, mu = mu) + 0.02,
            label = list("$.05$"), 
            color = "grey30",
            hjust = 0,
            fontfamily = title_theme$family,
            fontface = title_theme$face,
            size = title_theme$size)

# Add the curly bracket
p <- p + annotation_custom(b1)

# Add the label "\Delta(q)=1"
p <- p + draw_label(x = 10.5, y = 0.07,
                   label = "$\\Delta(q)=1$",
                   fontfamily = title_theme$family,
                   fontface = title_theme$face,
                   size = title_theme$size,
                   hjust = 0.5, vjust = -0.5,
                   color = bracketCol)

# Move the legend inside the plot
p <- p + common_theme() + 
  labs(y = 'Density', x = '$t$', linetype = NULL) + 
  theme(legend.position = "inside",
        legend.position.inside = c(.8,.775),
        legend.text = element_text(hjust=0),
        legend.box.background = element_rect(colour = panelBorderColour,
                                             linewidth = panelBorderWidth))
        #legend.box.background = element_rect(colour = "grey30"))
        #legend.margin=margin(t=-0.19, r=0.1, b=0, l=0, unit="cm"))

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
  xAxisLimits <- c(-nobs*2,nobs*3)
  #xAxisBreaks <- xAxisLimits[1]:xAxisLimits[2]
  expandConstant <- 0.05*(xAxisLimits[2]-xAxisLimits[1])
  xAxisLimits <- c(xAxisLimits[1]-expandConstant,xAxisLimits[2]+expandConstant)
  
  p_ <- ggplot(data.frame(x = seq(xAxisLimits[1], xAxisLimits[2], by = 0.1)), aes(x)) +
    geom_function(fun = lower_ml, args = list(eps = epsilon, n = nobs),
                  colour = color) +
    scale_x_continuous(expand = c(0,0),
                       breaks = c(-nobs*1.5, nobs*0.5, nobs*2.5)) +
    scale_y_continuous(n.breaks = 4) +
    common_theme() +
    labs(#title = 'density bounds for p( t | \u03B8 )',  
         x = '$t$', y = 'Density')
         #caption = paste0('n = ', nobs, ', \u03B5 = ', epsilon))
  
  
  p_both <- p_ + geom_function(fun = upper_ml, args = list(eps = epsilon, n = nobs), colour = color) #+ labs(subtitle = 'Upper and lower densities')
  
  p_lower <- p_ + labs(subtitle = 'Lower density only')
  
  return(list(both=p_both, lower=p_lower))
}

p <- plot_ml(nobs = 10, epsilon = 0.1, color = mypalette[2])$both
print(p)
tikz('bailie24Figs/marginal1TikZ.tex', width = 5*marginal12RelSize, height = 3.5*marginal12RelSize)
p
dev.off()

p <- plot_ml(nobs = 10, epsilon = 0.25, color = mypalette[2])$both
print(p)
tikz('bailie24Figs/marginal2TikZ.tex', width = 5*marginal12RelSize, height = 3.5*marginal12RelSize)
p
dev.off()

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

xAxisLimits <- c(1,10)
xAxisBreaks <- seq(2,10,2)
expandConstant <- 0.05*(xAxisLimits[2]-xAxisLimits[1])
xAxisLimits <- c(xAxisLimits[1]-expandConstant,xAxisLimits[2]+expandConstant)

dat2 %<>% add_row(t = xAxisLimits[2], 
                  value = c(
                    lower_rr(nt,epsilon)*(11-xAxisLimits[2])+lower_rr(nt+1,epsilon)*(xAxisLimits[2]-10),
                    upper_rr(nt,epsilon)*(11-xAxisLimits[2])+upper_rr(nt+1,epsilon)*(xAxisLimits[2]-10)
                  ),
                  name = c("lower", "upper"))


p_rr <- ggplot(dat2, aes(x = t, y = value, color = name)) +
  geom_line() + geom_point(data = dat2[dat2$t <= nt,]) +
  common_theme() + scale_color_manual(values = mypalette[c(2, 2)], guide = 'none') +
  scale_x_continuous(limits = xAxisLimits,
                     breaks = xAxisBreaks,
                     expand = c(0,0)) +
  labs(#title = 'density bounds for p( t | \u03B8 )',  
       x = '$|t|$', y = 'Density')
       #caption = paste0('max |t| = ', nt, ', \u03B5 = ', epsilon))


p_rr
tikz('bailie24Figs/marginal_rrTikZ.tex', width = 5, height = 3.5)
p_rr
dev.off()


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
  common_theme() +
  labs(x = '$\\theta$', y = 'Posterior density')#, title = 'density bounds for posterior p( \u03B8 | t )',
       #caption = paste0('prior for \u03B8 is Gamma(', alpha, ',', beta, '), clamping range is [', a0, ',', a1, '], \u03B5 = ', epsilon))

p

tikz('bailie24Figs/posteriorTikZ.tex', width = 5, height = 3.5)
p
dev.off()

############################################################################

