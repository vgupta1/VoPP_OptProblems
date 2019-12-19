library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)

dat = read_tsv("Results/")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "UniAnalytic", "Unimodal", "Symmetric"),
                      group_type = as.factor(str_c(Shape, isLB))
                      )

dat <- dat %>% mutate(Label = fct_recode(Shape,
               `Symmetric \\& Unimodal`="Symmetric", 
                `Unimodal ($m=1$)` = "Unimodal",
               `Unimodal` = "UniAnalytic") 
               )


#drop the symmetric stuff since we don't discuss in paper anymore
dat <- filter(dat, Shape != "Symmetric", 
                    Bound > 0)

g <- dat %>% ggplot(aes(Dev, Bound, color=Label, group=group_type)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("Coefficient of Deviation $D$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")
g


tikzDevice::tikz(file = "../../TeX/VoPP/Figures/FullMadPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()

###############
# Similar Crap for CV
###############
dat = read_tsv("Results/outCV__2_1_1_1.tab")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "UniAnalytic", "Unimodal", "Symmetric"),
                      group_type = as.factor(str_c(Shape, isLB))
)

dat <- dat %>% mutate(Label = fct_recode(Shape,
                                         `Symmetric \\& Unimodal`="Symmetric", 
                                         `Unimodal ($m=1$)` = "Unimodal",
                                         `Unimodal` = "UniAnalytic") 
)


#drop the symmetric stuff since we don't discuss in paper anymore
dat <- filter(dat, Shape != "Symmetric", 
              Bound > 0)

g <- 
  dat %>% ggplot(aes(CV, Bound, color=Label, group=group_type)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("Coefficient of Variation $C$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

g


tikzDevice::tikz(file = "../../TeX/VoPP/Figures/FullCVPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()  
