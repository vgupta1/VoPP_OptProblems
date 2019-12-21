library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)

#simple plots for the upper bound theorems

dat = read_tsv("Results/simplePlot__4_1_100.tab")

g <- dat %>% ggplot(aes(D, Bound)) + 
  geom_line(size=1) +
  theme_bw(base_size = 6) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("Coefficient of Deviation $D$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")
g

##Add some dividers to see where it splits
delta_l = 0.37542149115328854 # S = 4
delta_m = 0.5809402158035948  # S = 4
delta_h = .75

g <- g + geom_vline(xintercept=delta_l, linetype="dotted") + 
  geom_vline(xintercept=delta_m, linetype="dotted") + 
  annotate("text", label = "Low", 
           x = 0.15, y = 3.5, size=1.75) + 
  annotate("text", label = "Medium", 
           x = 0.49, y = 3.5, size=1.75) + 
  annotate("text", label = "High", 
           x = 0.7, y = 3.5, size=1.75) + 
  annotate("text", label = "$\\delta_L$", 
           x = delta_l - .025, y = 1, size=1.75) + 
  annotate("text", label = "$\\delta_M$", 
           x = delta_m - .025, y = 1, size=1.75) + 
  annotate("text", label = "$\\delta_H$", 
           x = delta_h + .03, y = 1, size=1.75)

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/simplePlot.tex", width=2, height=4/6 * 2 )
g
dev.off()


#####
#Similar plot but for inverse

g <- dat %>% ggplot(aes(D, 1/Bound)) + 
  geom_line(size=1) +
  theme_bw(base_size = 6) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("Coefficient of Deviation $D$") +
  ylab("Bound on $\\mathcal{R}_{SP} / \\mathcal{R}_{PP}$")
g

g<- g + geom_vline(xintercept=delta_l, linetype="dotted") + 
  geom_vline(xintercept=delta_m, linetype="dotted") + 
  annotate("text", label = "Low", 
           x = 0.15, y = .95, size=1.75) + 
  annotate("text", label = "Medium", 
           x = 0.49, y = .95, size=1.75) + 
  annotate("text", label = "High", 
           x = 0.7, y = .95, size=1.75) + 
  annotate("text", label = "$\\delta_L$", 
           x = delta_l - .025, y = .26, size=1.75) + 
  annotate("text", label = "$\\delta_M$", 
           x = delta_m - .025, y = .26, size=1.75) + 
  annotate("text", label = "$\\delta_H$", 
           x = delta_h + .03, y = .26, size=1.75)

tikzDevice::tikz(file = "../../TeX/VoPP/Figures/singlePriceGuarantee.tex", width=2, height=4/6 * 2 )
g
dev.off()


###Now the fun begins... Try to do it fo ra contour
dat = read_tsv("tempContour__2_11_100.tab")
dat = read_tsv("tempContour__1.1_4_100.tab")

#create a set of custom breaks that includes exp(1) for MHR
bks <- c(1.25, 1.5, 2, exp(1), 3, 3.5)
labs <- c("1.25", "1.50", "2.00", "e", "3.00", "3.50")

g <- ggplot() +
  stat_contour(data=dat, 
               aes(D, S, z = Bound, 
                   linetype= factor(..level.. == exp(1), levels = c(F, T) ), 
                   color= factor(..level..)),  
                    breaks=bks) + 
  scale_color_discrete(labels = labs) + 
    theme_bw(base_size = 6) + 
  theme(legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  scale_linetype_discrete(guide=FALSE) + 
  xlab("Coefficient of Deviation $D$") +
  ylab("Scale $S$")


tikzDevice::tikz(file = "../../TeX/VoPP/Figures/contourPlot.tex", width=2, height=4/6 * 2 )
g
dev.off()


