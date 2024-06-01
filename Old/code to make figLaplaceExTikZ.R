
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
eps <- 1
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
b1 <- bracketsGrob(5/11+0.004, 0.05, 0.591*10/11+0.004, 0.05, h=0.05, lwd=1, col=bracketCol)

#Making the plot
p <- ggplot(plotData, aes(x = t, y = px, linetype=factor(type)))

# Add vlines for q(x) and q(x) true values
# p <- p + geom_vline(xintercept = 10, linetype = "twodash", color = vlineCol, linewidth = 0.4)
# p <- p + geom_vline(xintercept = 11, linetype = "twodash", color = vlineCol, linewidth = 0.4)

# Add the lines
p <- p + geom_line(colour = col1) +
  scale_linetype_manual(values=c("solid", "dashed"), 
                        labels = c("$p_x(t)$",
                                   "$p_{x'}(t)$"))

p <- p + scale_x_continuous(breaks = 5:16,
                            labels = c(rep("", 4),"$q(x)-\\Delta(q)$", "$q(x)$", "$q(x')$", rep("", 5)))

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

# Add the label "\delta(q)"
p <- p + draw_label(x = 10.5, y = 0.03,
                    label = "$\\Delta(q)$",
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
        legend.margin=margin(t=-0.19, r=0.1, b=0, l=0, unit="cm"),
        panel.grid.minor.x = element_blank())


print(p)
tikz('bailie24Figs/laplaceTikZ.tex', width = 5, height = 3.5)
p
dev.off()