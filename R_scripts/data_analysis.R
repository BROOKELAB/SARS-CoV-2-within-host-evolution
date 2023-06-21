library(tidyverse)
library(rio)
library(here)

#overall average iSNV count
load("naive_snpcounts.RData")
load("vax_snpcounts.RData")
total.count <- rbind(naive.snpcounts,snpcounts.vax)$SNP_count
mean(total.count)
#25.7551

#naive vs vaccinated iSNV counts
load("naive_snpcounts.RData")
load("vax_snpcounts.RData")
t.test(naive.snpcounts$SNP_count,snpcounts.vax$SNP_count)

#naive mean = 14.06604 
#vax mean = 55.97561
#p-value = 3.106e-05
#likely because of difference in Ct thresholds

#overall average shared iSNV count
load("naive_shared.RData")
load("vax_shared.RData")
vax.shared <- vax.shared[-c(1,5,11),] #remove dates with only 1 timepoint
total.shared <- rbind(naive.shared,vax.shared)$Shared
mean(total.shared)
#3.25

#naive vs vaccinated shared iSNV counts
load("naive_shared.RData")
load("vax_shared.RData")
vax.shared <- vax.shared[-c(1,5,11),] #remove dates with only 1 timepoint
naive.shared <- naive.shared$Shared
vax.shared <- vax.shared$Shared
t.test(naive.shared,vax.shared)
#naive mean = 2.05
#vax mean = 6.25
#p-value = 0.1226

#naive vs vax ct 
naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")
t.test(naive.ct$ct,vax.ct$ct)
#naive = 23.79307 #vax = 25.36659
#p-value = 0.003251

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
#p-value = 0.001235
fisher.test(n.dat.cut)
#p-value = 0.2283

#naive iSNV count over time
load("naive_snpcounts.RData")
naive.ct <- import("naive_ct.xlsx")
naive.snpcounts <- naive.snpcounts%>%
  mutate(ct = naive.ct$ct)%>%
  filter(ct < 26) #filter for samples unlikely to have high #s of seq artifacts
naive.count.num <- naive.snpcounts
naive.count.num$day_of_infection <- as.numeric(naive.count.num$day_of_infection)
naive.dailysnp <- lm(SNP_count~day_of_infection,data = naive.count.num)
summary(naive.dailysnp)
#adjusted r-squared = 0.05443, p-value: 0.009203

#plot naive iSNV count vs day
naive.snpcounts$day_of_infection <- as.character(naive.snpcounts$day_of_infection)
ggplot(data = naive.snpcounts, aes(x=day_of_infection,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:10)))+
  scale_y_continuous(limits = c(0,100))+
  geom_abline(intercept = 4.4023,slope=2.1208,color="grey50")+
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
  filter(ct < 26)

vax.count.num  <- snpcounts.vax
vax.count.num$day_of_infection <- as.numeric(vax.count.num$day_of_infection)
vax.dailysnp <- lm(SNP_count~day_of_infection,data = vax.count.num)
summary(vax.dailysnp)
#adjusted r-squared = -0.03176 , p-value = 0.6141  

#plot vax iSNV count vs day
snpcounts.vax$day_of_infection <- as.character(snpcounts.vax$day_of_infection)
ggplot(data = snpcounts.vax, aes(x=day_of_infection,y=SNP_count))+
  geom_point(cex=4)+
  scale_x_discrete(limits = factor(c(1:11)))+
  scale_y_continuous(limits = c(0,100))+
  xlab("Day Post-Enrollment") +
  ylab("iSNV Count")+
  ggtitle("Vaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none",plot.title = element_text(size = 22))
ggsave("figs/vax_snps_vs_day.png")









