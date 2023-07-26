library(tidyverse)
library(rio)
library(here)

#overall average iSNV count
load("naive_snpcounts.RData")
load("vax_snpcounts.RData")
naive.snpcounts <- mutate(naive.snpcounts, "status" = "naive")
vax.snpcounts <- mutate(vax.snpcounts, "status" = "vax")
snpcounts <- rbind(naive.snpcounts, vax.snpcounts)
total.snps <- snpcounts$total_SNPs
mean(total.snps)
#32.74419

#naive vs vaccinated iSNV counts
mean(naive.snpcounts$total_SNPs) #26.83582
mean(vax.snpcounts$total_SNPs) #53.57895
total.test <- glm(total_SNPs ~ status, data = snpcounts, 
                  family = poisson(link = "log"))  
summary(total.test)
#coef = 0.69142 #p < 2e-16

#overall average shared iSNV count
shared <- snpcounts$shared_SNPs
mean(shared)
#3.02907

#naive vs vaccinated shared iSNV counts
mean(naive.snpcounts$shared_SNPs) #1.925373
mean(vax.snpcounts$shared_SNPs) #6.921053
shared.test <- glm(shared_SNPs ~ status, data = snpcounts, 
                  family = poisson(link = "log"))  
summary(shared.test)
#coeff = 1.27945 #p < 2e-16

#naive vs vax ct 
naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")
t.test(naive.ct$ct,vax.ct$ct)
#naive = 23.75664 #vax = 25.14421
#p-value = 0.01113

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
#naive = 5/60 = 0.0833
#vax = 16/71 = 0.225
#p-value = 0.03234

fisher.test(n.dat.cut)
#naive = 5/60 = 0.0833
#vax = 5/32 = 0.156
#p-value = 0.3087

#naive iSNV count over time
naive.snpcounts <- mutate(naive.snpcounts, "ct" = naive.ct$ct)
naive.count.num <- naive.snpcounts
naive.count.num$day_of_infection <- as.numeric(naive.count.num$day_of_infection)
naive.daily.test <- glm(total_SNPs ~ day_of_infection, data = naive.count.num,
            family = poisson(link = "log"))
summary(naive.daily.test)
#coef = 0.156180 #p < 2e-16

#plot naive iSNV count vs day
naive.snpcounts$day_of_infection <- as.character(naive.snpcounts$day_of_infection)
ggplot(data = naive.count.num, aes(x=day_of_infection,y=total_SNPs))+
  geom_point(cex=4)+
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+
  scale_x_continuous(limits = c(0,15), breaks = c(0,5,10,15))+
  xlab("Day post-enrollment") +
  ylab("iSNV count")+
  ggtitle("Unvaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t=12)),
        axis.title.y = element_text(margin = margin(r=15)),
        axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/naive_snps_vs_day.png")

#vax iSNV count over time
vax.snpcounts <- mutate(vax.snpcounts, "ct" = vax.ct$ct)
vax.count.num <- vax.snpcounts
vax.count.num$day_of_infection <- as.numeric(vax.count.num$day_of_infection)
vax.daily.test <- glm(total_SNPs ~ day_of_infection, data = vax.count.num,
                        family = poisson(link = "log"))
summary(vax.daily.test)
#coef =  0.054115 #p < 2e-16

#plot vax iSNV count vs day
vax.snpcounts$day_of_infection <- as.character(vax.snpcounts$day_of_infection)
ggplot(data = vax.count.num, aes(x=day_of_infection,y=total_SNPs))+
  geom_point(cex=4)+
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+
  scale_x_continuous(limits = c(0,15), breaks = c(0,5,10,15))+
  scale_y_continuous(limits = c(0,325), breaks = c(100,200,300))+
  xlab("Day post-enrollment") +
  ylab("iSNV count")+
  ggtitle("Vaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t=12)),
        axis.title.y = element_text(margin = margin(r=15)),
        axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/vax_snps_vs_day.png")

#saliva vs. nasal shared iSNV counts
saliva.ids <- lapply(paste0("sampleIDs_saliva/", list.files("sampleIDs_saliva")),
                     import)
saliva.ids <- bind_rows(saliva.ids)$V1
saliva.snps <- naive.ct[which(naive.ct$vdl_barcode %in% saliva.ids),]  
saliva.snps <- saliva.snps %>%
  mutate("status" = "saliva")%>%
  select(user_id, total_SNPs, shared_SNPs, status)

load("nasal_snpcounts.RData")
nasal.snps <- nasal.snpcounts %>%
  mutate("status" = "nasal") %>%
  dplyr::rename("user_id" = user_ID)%>%
  select(-day_of_infection)

mean(saliva.snps$shared_SNPs) #1.3125
mean(nasal.snps$shared_SNPs) #0.4

env.snps <- rbind(saliva.snps, nasal.snps)
env.test <- glm(shared_SNPs ~ status, data = env.snps, 
                family = poisson(link = "log"))  
summary(env.test)
#coef = 1.1882 #p = 7.28e-07






