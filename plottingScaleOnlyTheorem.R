library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)

#simple plots for the upper bound scale only theorem

dat = read_tsv("Results/out_tight_dist_scale_only__4_100.tab")
dat %>% filter(Fbar == 1.) %>% summarise(max(x))

g <- ggplot(data=dat, aes(x, Fbar)) + 
  geom_line() + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  xlab("t") +
  ylab("$P( V \\geq t)$") +
  geom_vline(xintercept=.271, linetype="dotted") + 
  annotate("text", label = "$\\alpha$", 
           x = .271 + .15, y = .125, size=1.75)

g


tikzDevice::tikz(file = "../../TeX/VoPP/Figures/scale_only_tightDistn.tex", 
                 width=1.56, height=1.56 )
g
dev.off()

################
################
################

dat = read_tsv("Results/out_scale_vary_s__1_10.tab")
g<- ggplot(data=dat, aes(x, Bound)) + 
  geom_line() + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  xlab("Scale (S)") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/scale_bound_vary_s.tex", 
                 width=1.56, height=1.56 )
g
dev.off()


################
################

dat = read_tsv("Results/out_scale_vary_M__0.1_1.tab")
g<- ggplot(data=dat, aes(M, Bound)) + 
  geom_line() + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  xlab("Margin (M)") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/scale_bound_vary_M.tex", 
                 width=1.56, height=1.56 )
g
dev.off()


################
################
dat = read_tsv("Results/out_scale_vary_ratio__0.1_1.tab")
g<- ggplot(data=dat, aes(x, Bound)) + 
  geom_line() + 
  theme_bw(base_size = 6) + 
  theme(legend.position="none") +
  xlab("$\\frac{M}{M + S -1}$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/scale_bound_vary_ratio.tex", 
                 width=1.56, height=1.56 )
g
dev.off()
