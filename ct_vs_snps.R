library(tidyverse)
library(rio)
library(here)

naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")

together.ct <-  full_join(naive.ct,vax.ct)
ggplot(data = together.ct, aes(x=ct,y=SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19))+
  xlab("N Gene Ct")+
  ylab("iSNV Count")
ggsave("figs/total_ctSNPs.png")
together.snpct.cor <- cor.test(together.ct$ct,together.ct$SNP_count,method="spearman",
                               exact = F)
#rho = 0.6350996, p-value < 2.2e-16



ct.compare <- t.test(naive.ct$ct,vax.ct$ct)


