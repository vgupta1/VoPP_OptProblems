## Plots of Tight Distributions for MAD CASE
###

library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)
library(lamW)
library(purrr)

####Plotting the various tight distributions
dat <- read_tsv("Results/tightDistMAD__4_100.tab")

sz = 1.5

#Low Dev
#pull out the jump point
pt1 <- filter(dat, x == 1) %>% select(Low) %>% .[[1]]
pt2 <- 0.18771074557664427 / log(4)  #hardcoded based on experiment

g <- ggplot(dat, aes(x, Low)) + 
  geom_line(data=filter(dat, x <= 1)) + 
  geom_line(data=filter(dat, x > 1)) + 
  annotate("point", x = 1, y = pt1, size=1.75) +
  annotate("point", x = 1, y = pt2, size=1.75, shape=21) + 
  annotate("text", x = 2.75, y = .85, size = 1.75,
           label="Low Deviation \n $D = \\frac{\\delta_L}{2} \\approx 0.188$") + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  ylab("$\\overline F_L(x)$") + xlab("x")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/tightDistMad_Low.tex", 
                 width=sz, height=sz )
g
dev.off()

#### Medium CAse
#Jump is at 0
pt <- filter(dat, 0 < x, x < .1) %>% select(Med) %>% .[[1,1]]
g <- ggplot(dat, aes(x, Med)) + 
  geom_line(data = filter(dat, x > 0)) +
  annotate("point", x = 0, y = 1, size=1.75) +
  annotate("point", x = 0, y = pt, size=1.75, shape=21) + 
  annotate("text", x = 2.75, y = .85, size = 1.75,
           label="Medium Deviation \n $D = \\frac{\\delta_L + \\delta_M}{2} \\approx 0.478$") + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  ylab("$\\overline F_M(x)$") + xlab("x")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/tightDistMad_Med.tex", 
                 width=sz, height=sz)
g
dev.off()

#High CAse
#Jump is at 0
pt <- filter(dat, 0 < x, x < .1) %>% select(High) %>% .[[1,1]]
g <- ggplot(dat, aes(x, High)) + 
  geom_line(data=filter(dat, x > 0)) +
  annotate("point", x = 0, y = 1, size=1.75) +
  annotate("point", x = 0, y = pt, size=1.75, shape=21) +
  annotate("text", x = 2.75, y = .85, size = 1.75,
           label="High Deviation \n $D = \\frac{\\delta_M + \\delta_H}{2} \\approx 0.665$") + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  ylab("$\\overline F_H(x)$") + xlab("x")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/tightDistMad_High.tex", 
                 width=sz, height=sz )
g
dev.off()

