library(tidyverse)
library(rio)
library(here)

#overall average iSNV count
load("naive_count_gather.RData")
load("vax_count_gather.RData")
total_count <- rbind(naive_count_gather,vax_count_gather)$SNP_count
total_count <- total_count[which(!is.na(total_count))]
mean(total_count)
#37.56

#naive vs vaccinated iSNV counts
load("naive_count_gather.RData")
load("vax_count_gather.RData")
naive.snpcount <- naive_count_gather$SNP_count[which(!is.na(naive_count_gather$SNP_count))]
vax.snpcount <- vax_count_gather$SNP_count[which(!is.na(vax_count_gather$SNP_count))]
t.test(naive.snpcount,vax.snpcount)

#naive mean = 33.38562
#vax mean = 51.73333
#p-value = 0.05645

#overall average shared iSNV count
load("naive_shared.RData")
load("vaccinated_shared.RData")
total_shared <- rbind(naive_shared,vaccinated_shared)$Shared
mean(total_shared)

#naive vs vaccinated shared iSNV counts
load("naive_shared.RData")
load("vaccinated_shared.RData")
naive_shared <- naive_shared$Shared
vaccinated_shared <- vaccinated_shared$Shared
snpcount.test.shared <- t.test(naive_shared,vaccinated_shared)
#naive mean = 6.65
#vax mean = 5.55
#p-value = 0.6734


#proportion of nucleocapsid iSNVs in naive vs vax participants
naive.ann <- import("naive_annotations.csv")
vax.ann <- import("vaccinated_annotations.csv")

n.dat <- data.frame(
  n = c(sum(naive.n==1),sum(vax.n==1)),
  other =  c(sum(naive.n==0),sum(vax.n==0)),
  row.names = c("naive","vax")
)
n.dat.cut <- data.frame(
  n = c(sum(naive.n==1),sum(vax.n.cut==1)),
  other =  c(sum(naive.n==0),sum(vax.n.cut==0)),
  row.names = c("naive","vax")
)
n.fisher <- fisher.test(n.dat)
#p-value = 1.918e-06
n.fisher.cut <- fisher.test(n.dat.cut)
#p-value = 0.03669

#naive iSNV count over time
load("naive_count_gather.RData")
naive.ct <- import("naive_ct.xlsx")
naive_count_gather <- naive_count_gather%>%
  filter(!is.na(SNP_count))%>%
  mutate(ct = naive.ct$ct)%>%
  filter(ct < 25) #filter for samples unlikely to have high #s of seq artifacts

naive_count_num <- naive_count_gather
naive_count_num$day <- as.numeric(naive_count_num$day)
naive_dailysnp <- lm(SNP_count~day,data = naive_count_num)
summary(naive_dailysnp)
#adjusted r-squared = 0.05007, p-value: 0.02255

#plot iSNV count vs day (FIG 2C)
naive_count_gather$day <- as.character(naive_count_gather$day)
ggplot(data = naive_count_gather, aes(x=day,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:9)))+
  scale_y_continuous(limits = c(0,100))+
  geom_abline(intercept = 2.6549,slope=2.1281,color="grey50")+
  xlab("Time Point") +
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/naive_snps_vs_day.png")

#vaccinated iSNV count over time
load("vax_count_gather.RData")
vax.ct <- import("vaccinated_ct.csv")
vax_count_gather <- vax_count_gather%>%
  filter(!is.na(SNP_count))%>%
  mutate(ct = vax.ct$ct)%>%
  filter(ct < 25)

vax_count_num  <- vax_count_gather
vax_count_num$day <- as.numeric(vax_count_num$day)
vax_dailysnp <- lm(SNP_count~day,data = vax_count_num)
summary(vax_dailysnp)
#adjusted r-squared = 0.2857, p-value = 0.006103  

#plot iSNV count vs day (FIG 2D)
vax_count_gather$day <- as.character(vax_count_gather$day)
ggplot(data = vax_count_gather, aes(x=day,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:7)))+
  scale_y_continuous(limits = c(0,100))+
  geom_abline(intercept = 0.1923,slope=4.3511,color="grey50")+
  xlab("Time Point") +
  ylab("iSNV Count")+
  ggtitle("Immune")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none",plot.title = element_text(size = 22))
ggsave("figs/vax_snps_vs_day.png")









