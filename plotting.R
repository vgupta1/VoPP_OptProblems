library(tidyverse)
library(ggplot2)
library(stringr)
library(tikzDevice)
library(forcats)

#dat = read_tsv("Results/outMAD__2_1_1_1.tab")
dat = read_tsv("Results/outMAD2__2_0.9_1_1.tab")


dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "UniAnalytic", "Unimodal", "Symmetric"),
                      group_type = as.factor(str_c(Shape, isLB))
                      )

dat <- dat %>% mutate(Label = fct_recode(Shape,
               `Arbitrary Shape` = "None",
               `Symmetric \\& Unimodal`="Symmetric", 
                `Unimodal w/mode $m=1$` = "Unimodal",
               `Unimodal` = "UniAnalytic") 
               )

#drop the symmetric stuff since we don't discuss in paper anymore
dat <- filter(dat, Shape != "Symmetric", 
                    Bound > 0)

g <- dat %>% ggplot(aes(Dev, Bound, color=Label)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  scale_y_continuous(limits=c(NA, 2.8)) +
  theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("(a) Coefficient of Deviation $D$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")
g

dat.labels <- 
  tribble(~Dev, ~Bound, ~Label, 
          .225,   2.77, "Thm. 2", 
          .225, 2.27, "Thm. 7", 
          .225, 1.8, "Thm. 8", 
          .225, 1.4, "Thm. 3", 
          .225, 1.05, "Ex. 1"
  )

g <- g + geom_text(aes(Dev, Bound, label=Label), 
              data=dat.labels, color="black", size=2
              )

tikzDevice::tikz(file = "../../-MS-Minor-Revision-VoPP/Figures/FullMadPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()

###############
# Similar Crap for CV
###############
dat = read_tsv("Results/outCV2__2_0.9_1_1.tab")


dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "UniAnalytic", "Unimodal", "Symmetric"),
                      group_type = as.factor(str_c(Shape, isLB))
)

dat <- dat %>% mutate(Label = fct_recode(Shape,
                                         `Arbitrary Shape` = "None", 
                                         `Symmetric \\& Unimodal`="Symmetric", 
                                         `Unimodal w/mode $m=1$` = "Unimodal",
                                         `Unimodal` = "UniAnalytic") 
)


#drop the symmetric stuff since we don't discuss in paper anymore
dat <- filter(dat, Shape != "Symmetric", 
              Bound > 0)

g <- 
  dat %>% ggplot(aes(CV, Bound, color=Label)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  scale_color_manual(values=c("#F8766D", "#619CFF"))+
    theme(legend.position=c(.235, .72), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("(b) Coefficient of Variation $C$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

g


dat.labels <- 
  tribble(~Dev, ~Bound, ~Label, 
          .5, 2.72, "Thm. 6",
          .5, 2.3, "Thm. 7", 
          .5, 1.5, "Thm. 8", 
          .5, 1.05, "Ex. EC.1"
  )

g<-   g + geom_text(aes(Dev, Bound, label=Label), 
                   data=dat.labels, color="black", size=2)


tikzDevice::tikz(file = "../../-MS-Minor-Revision-VoPP/Figures/FullCVPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()  




###############
# Similar Crap for GM
###############
dat = read_tsv("Results/outGM__2_0.9_1_1.tab")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "Unimodal"),
                      group_type = as.factor(str_c(Shape, isLB))
)

dat <- dat %>% mutate(Label = fct_recode(Shape,
                                         `Arbitrary Shape` = "None", 
                                         `Unimodal w/mode $m=1$` = "Unimodal") 
)


#Debug for now 
names(dat)[1] = "GM"

g <- 
  dat %>% ggplot(aes(GM, Bound, color=Label)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  scale_color_manual(values=c("#F8766D", "#619CFF"))+
  scale_y_continuous(limits=c(NA, 3.2)) +
    theme(legend.position=c(.7, .8), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("(c) Geometric Mean $B$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

g

dat.labels <- 
  tribble(~Dev, ~Bound, ~Label, 
          .85, 2.85, "Thm. 6",
          .85, 2.25, "Thm. 7", 
          .85, 1.85, "Thm. 8", 
          .85, 1.05, "Ex. EC.2"
  )

g <- g + geom_text(aes(Dev, Bound, label=Label), 
                  data=dat.labels, color="black", size=2)


tikzDevice::tikz(file = "../../-MS-Minor-Revision-VoPP/Figures/FullGMPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()  

#######
###############
# Similar Crap for IC
###############
dat = read_tsv("Results/ICBounds__2_0.9_1.0_1.0_0.8.tab")

dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#rename the shapes and reorder
dat <- dat %>% mutate(Shape = as.factor(Shape), 
                      Shape = fct_relevel(Shape, "None", "Unimodal"),
                      group_type = as.factor(str_c(Shape, isLB))
)

dat <- dat %>% mutate(Label = fct_recode(Shape,
                                         `Arbitrary Shape` = "None", 
                                         `Unimodal w/mode $m=1$` = "Unimodal") 
)

g <- 
  dat %>% ggplot(aes(IC, Bound, color=Label)) + 
  geom_line(aes(linetype=isLB), size=1) +
  theme_bw(base_size = 8) + 
  scale_linetype_discrete(guide=FALSE) + 
  scale_color_manual(values=c("#F8766D", "#619CFF"))+
  scale_y_continuous(limits=c(NA, NA)) +
  theme(legend.position=c(.7, .8), 
        legend.title=element_blank(), 
        legend.text=element_text(size=6)) +
  xlab("(d) Fraction of Market $q$") +
  ylab("Bound on $\\mathcal{R}_{PP} / \\mathcal{R}_{SP}$")

g

dat.labels <- 
  tribble(~Dev, ~Bound, ~Label, 
          .8, 1.8, "Thm. 6",
          .8, 1.5, "Thm. 7", 
          .9, 1.25, "Thm. 8", 
          .8, 1.05, "Ex. EC.3"
  )

g <- g + geom_text(aes(Dev, Bound, label=Label), 
                   data=dat.labels, color="black", size=2)

g


tikzDevice::tikz(file = "../../-MS-Minor-Revision-VoPP/Figures/FullICPlot.tex", width=3.125, height=4/6 * 3.125 )
g
dev.off()  

