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
options(scipen = 10)

mypalette <- brewer.pal(3, 'Set1')

############################################################################
#                               paperLaplaceEx                             #
############################################################################
# Figure 1 in Bailie and Gong 24+

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
  #eps/2*exp(eps*(n/2 - abs(x - n/2)))
  dlaplace(x = x, mu, sqrt(2)/eps)
}

second_density <- function(x, eps, mu){
  #eps/2*exp(eps*(n/2 - abs(x - n/2)))
  dlaplace(x = x, mu+1, sqrt(2)/eps)
}

#Making the plot data:
eps <- 1
mu <- 10
t <- seq(5,16,0.01)
firstDensity <- one_density(x = t, eps = eps, mu = mu)
secondDensity <- second_density(x = t, eps = eps, mu = mu)
numPoints <- length(t)

plotData <- tibble(t = rep(t, 2), px = c(firstDensity, secondDensity),
                   type = c(rep(0,numPoints), rep(1, numPoints)))


#Variables for formatting the plot:
#text <- "x\'"
col1 <- mypalette[2] #"blue" #"#c6494a" #"#cc0033" 
col2 <- mypalette[2] #"blue" #"#518b7b" #"#478C5C" 
bracketCol <- "grey30" #"#bac075" #"#BACC81"
b1 <- bracketsGrob(5/11+0.004, 0.05, 0.591*10/11+0.004, 0.05, h=0.05, lwd=1, col=bracketCol)


#Making the plot
p <- ggplot(plotData, aes(x = t, y = px, linetype=factor(type)))

# Add vlines for q(x) and q(x) true values
p <- p + geom_vline(xintercept = 10, linetype = "solid", color = "grey80") #This must come first so it goes underneath the other lines
p <- p + geom_vline(xintercept = 11, linetype = "solid", color = "grey80")

# Add the lines
p <- p + geom_line(colour = col1) +
  scale_linetype_manual(values=c("solid", "twodash"), 
                        labels = list("$p_x$",
                                      "$p_{x'}$")) #NOTE LABEL MUST BE LIST (APPARENTLY) (at least for latex2exp)
# labels = c(expression(p[x]~with~q(x) == 10),
#            bquote(p[.(text)] * " with q(" * .(text) * ") = 11")))


p <- p + scale_x_continuous(breaks = c(10,11),# seq(min(t)+1, max(t)-1, 2), 
                            labels = c("$q(x)$",
                                       "$q(x')$")) +
  scale_y_continuous(minor_breaks = seq(0, 0.5, 0.1))
#scale_y_continuous(labels = NULL, breaks = NULL)#limits = c(0,0.6))

#Add the ratio of densities
p <- p + geom_segment(x = 8.5975, xend = 8.5975, y = 0, yend = one_density(x = 8.5975, eps = eps, mu = mu),#-0.001,
                      linetype = "solid", color = "grey30", linewidth = 0.4) +
  geom_segment(x = 8.6975, xend = 8.6975, y = 0, yend = second_density(x = 8.6975, eps = eps, mu = mu),#-0.001,
               linetype = "solid", color = "grey30", linewidth = 0.4)

p <- p + geom_text(x = 8.65, y = 0.1,#one_density(x = 9, eps = eps, mu = mu) + 0.02, 
                   label = list("$.05\\exp(\\epsilon)$"), 
                   show.legend = F, color = "grey30", size = 4,
                   hjust = 0) +
  geom_text(x = 8.75, y = 0.03,#one_density(x = 9, eps = eps, mu = mu) + 0.02,
            label = list("$.05$"), 
            show.legend = F, color = "grey30", size = 4,
            hjust = 0)

# Add the curly bracket
p <- p + annotation_custom(b1)

# Add the label "\delta(q)"
p <- p + geom_text(x = 10.5, y = 0.023,
                   label = list("$\\Delta(q)$"),
                   size = 3, hjust = 0.5, vjust = -0.5,
                   show.legend = FALSE,
                   color = bracketCol)

# Move the legend inside the plot
p <- p + theme_bw() + labs(y = 'density', x = '$t$', linetype = NULL, #Add some x and y labels
                           caption = "$\\epsilon = 1$",
                           title = "Two densities of the Laplace mechanism") + 
  theme(legend.position = c(.8,.775),
        legend.text.align = 0,
        legend.box.background = element_rect(colour = "grey30"),
        legend.margin=margin(t=-0.19, r=0.1, b=0, l=0, unit="cm"),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())


#print(p)
tikz('paperLaplaceExTikZ.tex', width = 5,height = 3.5)
p
dev.off()


############################################################################