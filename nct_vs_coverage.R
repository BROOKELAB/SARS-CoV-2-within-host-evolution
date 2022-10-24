library(tidyverse)
library(rio)
library(here)

nct <- import("nct_vs_coverage.xlsx")
nct.lm <- lm(Coverage~`N gene ct`,nct)
#p-value: 0.02359

ggplot(data = nct, aes(x = `N gene ct`,y = Coverage))+
  geom_point(size=4)+
  geom_errorbar(aes(x=`N gene ct`,ymin=Coverage-SD,ymax=Coverage+SD))+
  xlim(c(20,35)) +
  xlab("N Gene Ct")+
  #scale_y_log10()+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19))
ggsave("figs/ct_vs_coverage.png")
