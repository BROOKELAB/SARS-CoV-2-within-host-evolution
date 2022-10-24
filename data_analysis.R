library(tidyverse)
library(rio)
library(here)

#overall average iSNV count
load("naive_snpcounts.RData")
load("vax_snpcounts.RData")
total.count <- rbind(naive.snpcounts,snpcounts.vax)$SNP_count
mean(total.count)
#37.47

#naive vs vaccinated iSNV counts
load("naive_snpcounts.RData")
load("vax_snpcounts.RData")
t.test(naive.snpcounts$SNP_count,snpcounts.vax$SNP_count)

#naive mean = 33.30065
#vax mean = 51.66667
#p-value = 0.0559

#overall average shared iSNV count
load("naive_shared.RData")
load("vax_shared.RData")
vax.shared <- vax.shared[-c(1,6,12),] #remove dates with only 1 timepoint
total.shared <- rbind(naive.shared,vax.shared)$Shared
mean(total.shared)
#6.31

#naive vs vaccinated shared iSNV counts
load("naive_shared.RData")
load("vax_shared.RData")
vax.shared <- vax.shared[-c(1,6,12),] #remove dates with only 1 timepoint
naive.shared <- naive.shared$Shared
vax.shared <- vax.shared$Shared
t.test(naive.shared,vax.shared)
#naive mean = 6.65
#vax mean = 5.55
#p-value = 0.6734

#naive vs vax ct
naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")
t.test(naive.ct$ct,vax.ct$ct)
#naive = 23.79307 #vax = 25.14644
#p-value = 0.007773

#proportion of nucleocapsid iSNVs in naive vs vax participants
naive.ann <- import("naive_annotations.csv")
vax.ann <- import("vax_annotations.csv")

naive.ann <- naive.ann %>%
  mutate(naive.n = NA)
vax.ann <- vax.ann %>%
  mutate(vax.n = NA)

naive.ann$naive.n[which(naive.ann$GENE == "N")] <- 1
naive.ann$naive.n[which(naive.ann$GENE != "N")] <- 0
vax.ann$vax.n[which(vax.ann$GENE == "N")] <- 1
vax.ann$vax.n[which(vax.ann$GENE != "N")] <- 0

naive.n <- naive.ann$naive.n
vax.n <- vax.ann$vax.n
vax.n.cut <- vax.ann$vax.n[which(vax.ann$ID != "471876")]

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
fisher.test(n.dat)
#p-value = 2.618e-06
fisher.test(n.dat.cut)
#p-value = 0.03907

#naive iSNV count over time
load("naive_snpcounts.RData")
naive.ct <- import("naive_ct.xlsx")
naive.snpcounts <- naive.snpcounts%>%
  mutate(ct = naive.ct$ct)%>%
  filter(ct < 25) #filter for samples unlikely to have high #s of seq artifacts
naive.count.num <- naive.snpcounts
naive.count.num$day_of_infection <- as.numeric(naive.count.num$day_of_infection)
naive.dailysnp <- lm(SNP_count~day_of_infection,data = naive.count.num)
summary(naive.dailysnp)
#adjusted r-squared = 0.05581, p-value: 0.01671

#plot naive iSNV count vs day
naive.snpcounts$day_of_infection <- as.character(naive.snpcounts$day_of_infection)
ggplot(data = naive.snpcounts, aes(x=day_of_infection,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:10)))+
  scale_y_continuous(limits = c(0,100))+
  geom_abline(intercept = 2.4448,slope=1.8561,color="grey50")+
  xlab("Day Post-Enrollment") +
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/naive_snps_vs_day.png")

#vaccinated iSNV count over time
load("vax_snpcounts.RData")
vax.ct <- import("vaccinated_ct.csv")
snpcounts.vax <- snpcounts.vax%>%
  mutate(ct = vax.ct$ct)%>%
  filter(ct < 25)

vax.count.num  <- snpcounts.vax
vax.count.num$day_of_infection <- as.numeric(vax.count.num$day_of_infection)
vax.dailysnp <- lm(SNP_count~day_of_infection,data = vax.count.num)
summary(vax.dailysnp)
#adjusted r-squared = 0.001935 , p-value = 0.3198  

#plot vax iSNV count vs day
snpcounts.vax$day_of_infection <- as.character(snpcounts.vax$day_of_infection)
ggplot(data = snpcounts.vax, aes(x=day_of_infection,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:11)))+
  scale_y_continuous(limits = c(0,100))+
  xlab("Day Post-Enrollment") +
  ylab("iSNV Count")+
  ggtitle("Immune")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none",plot.title = element_text(size = 22))
ggsave("figs/vax_snps_vs_day.png")









