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
#29.40116

#naive vs vaccinated iSNV counts
mean(naive.snpcounts$total_SNPs) #23.20149
mean(vax.snpcounts$total_SNPs) #51.26316
total.test <- glm(total_SNPs ~ status, data = snpcounts, 
                  family = poisson(link = "log"))  
summary(total.test)
#coef = 0.79276 #p < 2e-16

#overall average shared iSNV count
shared <- snpcounts$shared_SNPs
mean(shared)
#2.348837

#naive vs vaccinated shared iSNV counts
mean(naive.snpcounts$shared_SNPs) #1.328358
mean(vax.snpcounts$shared_SNPs) #5.947368
shared.test <- glm(shared_SNPs ~ status, data = snpcounts, 
                  family = poisson(link = "log"))  
summary(shared.test)
#coeff = 1.49901 #p < 2e-16

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
#naive = 3/44 = 0.06818182
#vax = 13/64 = 0.203125
#p-value = 0.05915

fisher.test(n.dat.cut)
#naive = 3/44 = 0.06818182
#vax = 3/26 = 0.1153846
#p-value = 0.6635

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

mean(saliva.snps$shared_SNPs) #0.9166667
mean(nasal.snps$shared_SNPs) #0.3166667

env.snps <- rbind(saliva.snps, nasal.snps)
env.test <- glm(shared_SNPs ~ status, data = env.snps, 
                family = poisson(link = "log"))  
summary(env.test)
#coef = 1.0629 #p = 0.000108






