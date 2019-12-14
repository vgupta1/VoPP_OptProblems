library(tidyverse)
library(ggplot2)
library(stringr)

dat = read_tsv("Results/temp__2_1_1_1.tab")


dat %>% group_by(Shape, isLB) %>%
  summarize(avg = mean(Time))

#Make unique grouping variable
dat <- mutate(dat, group_type = str_c(Shape, isLB))


dat %>% ggplot(aes(Dev, Bound, color=Shape, group=group_type)) + 
  geom_line(aes(linetype=isLB)) + theme_bw() + 
  scale_linetype_discrete(guide=FALSE) + 
  theme(legend.position=c(.25, .75))