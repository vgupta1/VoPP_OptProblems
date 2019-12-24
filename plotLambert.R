#### 
# Simple Plot of Lambert-W and Chaz bounds
library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)
library(lamW)
library(purrr)

sz = 2

dat <- tibble(x=seq(-1/exp(1), 5, length.out=100))
dat <-dat %>% mutate(lam0 = lambertW0(x), 
                     lam1 = lambertWm1(x))

g <- ggplot(aes(x=x, y=lam1), data=dat, size=3) + 
  geom_line() + 
  geom_line(aes(y=lam0), linetype="dashed") + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  ylab("$W(x)$") +
  annotate("text", label = "$W_{-1}(x)$", 
           x = .4, y = -5, size=1.75) + 
  annotate("text", label = "$W_{0}(x)$", 
           x = .4, y = -.2, size=1.75)

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/lamWFuncs.tex", 
                 width=sz, height=4/6 * sz )
g
dev.off()

dat <- tibble(x=seq(0, 1, length.out=100))
dat <- dat %>% mutate(lam = lambertWm1(-x/exp(1)), 
                      ub = -1 - sqrt( 2 * log(1/x)) - 2/3 * log(1/x), 
                      lb = -1 - sqrt( 2 * log(1/x)) - log(1/x))

g <- ggplot(data=dat, aes(x, lam)) + 
  geom_line() + 
  geom_line(aes(y=ub), linetype="dotted") + 
  geom_line(aes(y=lb), linetype="dotted") + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  ylab("$W_1(-x/e)$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/lamWFunWithBounds.tex", 
                 width=sz, height=4/6 * sz )
g
dev.off()
