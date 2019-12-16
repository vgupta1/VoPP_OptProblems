library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)

dat = read_tsv("Results/temp__2_1_1_1.tab")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

dat <- dat %>% mutate(Shape = fct_recode(Shape, `Theorem 3`="UniAnalytic"), 
                      group_type = str_c(Shape, isLB))

g <- dat %>% ggplot(aes(Dev, Bound, color=Shape, group=group_type)) + 
  geom_line(aes(linetype=isLB)) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank()) +
  xlab("Coef. of Deviation $D$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/FullMadPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()

###############
# Similar Crap for CV
###############
dat = read_tsv("Results/tempCV__2_1_1_1.tab")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

dat <- dat %>% mutate(group_type = str_c(Shape, isLB))
g <- 
  dat %>% ggplot(aes(CV, Bound, color=Shape, group=group_type)) + 
  geom_line(aes(linetype=isLB)) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank()) +
  xlab("Coef. of Variation $C$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/FullCVPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()  
