#' ------------------
#'  Differential Privacy and Interval of Measures
#'  author: "Robin Gong"
#'  date: "11/11/2018"
#'  -----------------

setwd('~/Documents/various_projects/DP-IoM/')
library(ggplot2)

#'  ---
#'  Example 2, Lavine (1991).
#'  ---

gk <- function(theta, k, x = 1, xobs = 1, e = 0.1, lower.tail = F){
  # g(\theta) as a function of k;  x: desired interval boundary 
  # lower.tail = T: [0, x]; F: (x, +\infty),
  # xobs: observed data, 
  # e: epsilon in the DP-IoM: U = exp(e)*P, L = exp(-e)*P
  theta <- as.vector(theta)
  pe <- pexp(x, theta, lower.tail = lower.tail)
  q <- pmin(pe*exp(e), (1-exp(-e))+pe*exp(-e))
  #q <- pmin(pe/e, (1-e)+pe*e)
  gk_ <- (q - k)*dexp(x = xobs, theta)*exp(e)^ifelse((q-k) > 0, 1, -1)
  return(gk_)
}


## Some explorations.

if (F){
k0 = 1; e0 = 0.1; xobs = 1
vt = seq(0.1, 10, by = 0.1)
temp_ <- gk(vt, k=k0, x = 1, xobs=xobs, e=e0, lower.tail = F)
plot(vt, temp_, type = 'l')
set.seed(1)
nsim = 1e5; e0 = 0.1; xobs = 1
stheta <- rexp(nsim, 1)
vk = seq(0, 1, by = 0.005)
#vk = seq(1e-3, 1, by = 1e-3)
gk_est <- sapply(vk, function(kk){mean(gk(theta = stheta, k=kk, 
                                          x=1, xobs=xobs, e=e0, 
                                          lower.tail = T))})
plot(vk, gk_est, type = 'l')
(ind <- sum(gk_est>0)); gk_est[ind]
}

## Bisection algorithm to find the suitable $k$ value.
## For now suppose the prior $\pi_m$ is precise, that is, 
## $\theta \sim Exp(1)$ and the IoM parameter $\delta = 0$.

phi.bisect <- function(x, xobs, e, left = 0, right = 1, lower = F,
                       tol = 1e-7, max.iter = 5e3, nsim = 1e5){
  
  # bisection method to find phi.lower and phi.upper
  # using Bayesian precise prior.
  # left/right: boundary k values
  
  stheta <- rexp(nsim, rate = 1)
  ck.mc <- function(kk, t){
    # calculate \overbar{c} = \underbar{c}(k) w/ monte carlo
    mean(gk(theta = t, k=kk, x=x, xobs=xobs, e=e, lower.tail = lower))
  }
  c.left <- ck.mc(kk = left, t = stheta)
  c.right <- ck.mc(kk = right, t = stheta)
  
  if(c.left*c.right >= 0){
    cat('check initial conditions')
    break
  } 
  
  iter <- 0
  while ((abs(c.left - c.right) > tol)){
    k_ <- (left + right)/2
    c_ <-  ck.mc(kk = k_, t = stheta)
    if (c_*c.right > 0){
      right <- k_
      c.right <- c_
    } else {
      left <- k_
      c.left <- c_
    }
    
    iter <- iter + 1
    if (iter > max.iter){
      cat('max iteration reached')
      break
    }
  }
  
  if (lower == T){
    phi = 1 - max(left, right)
  } else {
    phi = max(left, right)
  }
  
  rt <- list(phi = phi, c.estimate = c.right, n.iter = iter, 
             params = list(x = x, xobs = xobs, epsilon = e, 
                           delta = 0, max.iter = max.iter, 
                           lower.tail = lower))
  return(rt)
}

(value <- phi.bisect(x = 1, xobs = 1, e = 0.1, lower = F))
(value <- phi.bisect(x = 1, xobs = 1, e = 0.1, lower = T))

## Calculate for a range of values.

vepsilon = c(0, 0.2, 0.5)
params <- expand.grid(x  = c(0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2:4),
                      T_obs = c(0.5, 2),
                      epsilon = vepsilon, 
                      phi.lower = NA,
                      phi.upper = NA)

for (i in 1:nrow(params)){
  pars <- params[i, ]
  temp2_ <- phi.bisect(x = pars$x, xobs = pars$T_obs, e = pars$epsilon, lower = F,
                       nsim = 1e4)
  temp_ <- phi.bisect(x = pars$x, xobs = pars$T_obs, e = pars$epsilon, lower = T,
                      nsim = 1e4)
  params$phi.lower[i] <- temp_$phi
  params$phi.upper[i] <- temp2_$phi
}

params$epsilon <- factor(params$epsilon, labels = vepsilon)

#data = subset(params, subset = (x_obs == 2))
p2 <- ggplot(params, 
            aes(x = x, fill = epsilon, color = epsilon)) +
  facet_grid(T_obs~epsilon, labeller = "label_both") +
  geom_ribbon(aes(ymin = phi.lower, ymax = phi.upper), alpha = 0.3) +
  labs(x = 't', 
       y = 'P(T > t)',
       title = expression(paste('Estimated ', underline(phi),' and ' , bar(phi),': posterior predictive upper and lower probabilities with Exp(1) prior' ))) +
  theme_bw() + 
  scale_fill_discrete(guide = 'none') +
 # scale_fill_discrete(expression(paste(epsilon,': sampling IoM'))) +
  scale_color_discrete(guide = 'none')

pdf('dp_IoM_lavine.pdf', width = 9, height = 5.5)
p2
dev.off()

#if (F){
#pdf('Lavine1991_ex2.pdf', width = 9, height = 3)
#p
#dev.off()
#}




phi.bisect2 <- function(x, xobs, e, left = 0, right = 1, lower = F,
                       tol = 1e-7, max.iter = 5e3, nsim = 1e5){
  
  # bisection method to find phi.lower and phi.upper
  # using Bayesian precise prior: gamma(2,2)
  # left/right: boundary k values
  
  #stheta <- rexp(nsim, 1)
  #stheta <- rexp(nsim, rate = 5)
  stheta <- rgamma(nsim, shape = 6, rate = 2)
  ck.mc <- function(kk, t){
    # calculate \overbar{c} = \underbar{c}(k) w/ monte carlo
    mean(gk(theta = t, k=kk, x=x, xobs=xobs, e=e, lower.tail = lower))
  }
  c.left <- ck.mc(kk = left, t = stheta)
  c.right <- ck.mc(kk = right, t = stheta)
  
  if(c.left*c.right >= 0){
    cat('check initial conditions')
    break
  } 
  
  iter <- 0
  while ((abs(c.left - c.right) > tol)){
    k_ <- (left + right)/2
    c_ <-  ck.mc(kk = k_, t = stheta)
    if (c_*c.right > 0){
      right <- k_
      c.right <- c_
    } else {
      left <- k_
      c.left <- c_
    }
    
    iter <- iter + 1
    if (iter > max.iter){
      cat('max iteration reached')
      break
    }
  }
  
  if (lower == T){
    phi = 1 - max(left, right)
  } else {
    phi = max(left, right)
  }
  
  rt <- list(phi = phi, c.estimate = c.right, n.iter = iter, 
             params = list(x = x, xobs = xobs, epsilon = e, 
                           delta = 0, max.iter = max.iter, 
                           lower.tail = lower))
  return(rt)
}

value <- phi.bisect2(x = 1, xobs = 2, e = exp(0.2), lower = F)
value2 <- phi.bisect2(x = 1, xobs = 2, e = exp(0.2), lower = T)
value$phi
value2$phi
value <- phi.bisect(x = 1, xobs = 2, e = exp(0.2), lower = F)
value2 <- phi.bisect(x = 1, xobs = 2, e = exp(0.2), lower = T)
value$phi
value2$phi

## Calculate for a range of values.

vepsilon = c(0, 0.2, 0.5)
params2 <- expand.grid(x  = c(0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2:4),
                      T_obs = c(0.5, 2),
                      epsilon = vepsilon, 
                      phi.lower = NA,
                      phi.upper = NA)


for (i in 1:nrow(params2)){
  pars <- params2[i, ]
  temp2_ <- phi.bisect2(x = pars$x, xobs = pars$T_obs, e = pars$epsilon, lower = F,
                        nsim = 2e4)
  temp_ <- phi.bisect2(x = pars$x, xobs = pars$T_obs, e = pars$epsilon, lower = T,
                       nsim = 2e4)
  params2$phi.lower[i] <- temp_$phi
  params2$phi.upper[i] <- temp2_$phi
}

params2$epsilon <- factor(params2$epsilon, labels = vepsilon)

p3 <- ggplot(data = params2, 
             aes(x = x, fill = epsilon, color = epsilon)) +
  facet_grid(T_obs~epsilon, labeller = "label_both") +
  geom_ribbon(aes(ymin = phi.lower, ymax = phi.upper), alpha = 0.3) +
  labs(x = 't', 
       y = 'P(T > t)',
       title = expression(paste('Estimated ', underline(phi),' and ' , bar(phi),': posterior predictive upper and lower probabilities with Gamma(6, 2) prior' ))) +
  theme_bw() + 
  scale_fill_discrete(guide = 'none') +
  # scale_fill_discrete(expression(paste(epsilon,': sampling IoM'))) +
  scale_color_discrete(guide = 'none')

pdf('dp_IoM_lavine_gamma62.pdf', width = 9, height = 5.5)
p3
dev.off()


params$prior = 'Exp(1)'
params2$prior = 'Gamma(6, 2)'
param <- rbind(params, params2)

p4 <- ggplot(data = subset(param, subset = (T_obs == 2)), 
             aes(x = x, fill = epsilon, color = epsilon)) +
  facet_grid(prior~epsilon, labeller = "label_both") +
  geom_ribbon(aes(ymin = phi.lower, ymax = phi.upper), alpha = 0.3) +
  labs(x = 't', 
       y = 'P(T > t)',
       title = expression(paste('Estimated ', underline(phi),' and ' , bar(phi),': posterior predictive upper and lower probabilities (', T[obs], ' = 2)' ))) +
  theme_bw() + 
  scale_fill_discrete(guide = 'none') +
  # scale_fill_discrete(expression(paste(epsilon,': sampling IoM'))) +
  scale_color_discrete(guide = 'none')

pdf('dp_IoM_lavine_both_priors_flat.pdf', width = 9, height = 4)
p4
dev.off()

###### JUNK ######

if (F){
  params <- expand.grid(x  = c(0.05, 0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2:4),
                        x_obs = c(.5, 1, 2),
                        epsilon = c(0.1, 0.5, 1), 
                        phi = NA)
  for (i in 1:nrow(params)){
    pars <- params[i, ]
    temp_ <- phi.bisect(x = pars$x, xobs = pars$x_obs, e = pars$epsilon)
    params$phi[i] <- temp_$phi
  }
  
  p <- ggplot(data = subset(params, subset = (x <= 4)), 
              aes(x = x, y = phi, color = as.factor(epsilon))) +
    facet_grid(.~x_obs, labeller = "label_both") +
    geom_line() +
    labs(x = 't', 
         y = expression(paste('estimated ',phi,' = P(x > t | ', x[obs], ')')),
         title = 'Lavine (1991) Example 2') +
    theme_bw() +
    scale_color_discrete(expression(paste(epsilon,': sampling IoM')))
  p
}

## --- ##
