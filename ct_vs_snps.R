library(tidyverse)
library(rio)
library(here)

#ct vs total snps
naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")
vax.ct <- vax.ct[,-1]
together.ct <-  full_join(naive.ct,vax.ct)

ggplot(data = together.ct, aes(x=ct,y=SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19))+
  xlab("N Gene Ct")+
  ylab("iSNV Count")
ggsave("figs/total_ctSNPs.png")
cor.test(together.ct$ct,together.ct$SNP_count,method="spearman",
                               exact = F)
#rho = 0.6350996, p-value < 2.2e-16

#ct vs shared snps
load("naive_daily_shared.RData")
load("vax_daily_shared.RData")

naive.daily.shared <- unlist(naive.daily.shared)
vax.daily.shared <- unlist(vax.daily.shared)

#only 1 timept and therefore no shared snps
vax.ct <- vax.ct[-(which(vax.ct$user_id == "461913" | vax.ct$user_id == "481672" |
                   vax.ct$user_id == "487941")),]
vax.daily.shared <- vax.daily.shared[-c(which(names(vax.daily.shared) == "461913" | 
                                                names(vax.daily.shared) == "481672" | 
                                                names(vax.daily.shared) == "487941"))]

naive.ct.shared <- bind_cols(naive.ct,naive.daily.shared)
colnames(naive.ct.shared)[6] <- "shared_SNP_count"
vax.ct.shared <- bind_cols(vax.ct, vax.daily.shared)
colnames(vax.ct.shared)[6] <- "shared_SNP_count"
total.ct.shared <- bind_rows(naive.ct.shared, vax.ct.shared)

#remove outliers
mean <- mean(total.ct.shared$shared_SNP_count)
std <- sd(total.ct.shared$shared_SNP_count)

t.min <-  mean-(3*std)
t.max <- mean+(3*std)

which(total.ct.shared$shared_SNP_count < t.min | 
        total.ct.shared$shared_SNP_count > t.max)
total.ct.cut <- total.ct.shared[-which(total.ct.shared$shared_SNP_count < t.min | 
                                             total.ct.shared$shared_SNP_count > t.max),]

#plot rho and p-value for different thresholds
cor.matrix <- as.data.frame(matrix(nrow = 6, ncol = 3))
colnames(cor.matrix) <- c("ct","rho","p-value")
cor.matrix$ct <- c(30,29,28,27,26,25)
total.ct.filter <- list()
for(i in 1:nrow(cor.matrix)){
  total.ct.filter[[i]] <- total.ct.cut[total.ct.cut$ct < cor.matrix$ct[i],]
  cor.matrix[i,2] <- (cor.test(total.ct.filter[[i]]$ct,
                               total.ct.filter[[i]]$shared_SNP_count,
                              method="spearman",
                              exact = F))$estimate
  cor.matrix[i,3] <- (cor.test(total.ct.filter[[i]]$ct,
                               total.ct.filter[[i]]$shared_SNP_count,
                               method="spearman",
                               exact = F))$p.value
}

cor.matrix.gather <- cor.matrix %>%
  gather(key = "stat", value = "value", - ct)

ggplot(data = cor.matrix.gather, aes(x = ct, y = value, color = stat))+
  geom_point(cex = 3)+
  geom_line()+
  geom_abline(slope = 0, intercept = 0.05, linetype = "dashed", color = "grey20")+
  ylim(c(-.2,1))+
  scale_x_reverse()+
  scale_color_manual(values = c("red","black"),
                     labels = c("P-value", "Rho"),
                     guide = guide_legend(title = "Statistic", reverse = T))+
  xlab("Ct")+
  ylab("Value")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/correlation_stats.png")
  
  
#plots for different thresholds
total.ct28 <- total.ct.cut[total.ct.cut$ct < 28,]
total.ct27 <- total.ct.cut[total.ct.cut$ct < 27,]
total.ct26 <- total.ct.cut[total.ct.cut$ct < 26,]
total.ct25 <- total.ct.cut[total.ct.cut$ct < 25,]

(length(total.ct28$user_id) - length(total.ct25$user_id))/(length(total.ct28$user_id))
#0.4216216 (ie 42% reduction between ct28 and ct25 datasets)

ggplot(data = total.ct28, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,28),breaks = c(14,16,18,20,22,24,26,28))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 28")
ggsave("figs/plot28.png")

ggplot(data = total.ct27, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,28),breaks = c(14,16,18,20,22,24,26,28))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 27")
ggsave("figs/plot27.png")

ggplot(data = total.ct26, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,26),breaks = c(14,16,18,20,22,24,26))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 26")
ggsave("figs/plot26.png")

ggplot(data = total.ct25, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,26),breaks = c(14,16,18,20,22,24,26))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 25")
ggsave("figs/plot25.png")

#monte carlo simulation of downsampled ct = 28 dataset
set.seed(100)
randoms <- list()
for(i in 1:100){
  randoms[[i]] <- sample(1:length(total.ct28$shared_SNP_count), 
                         length(total.ct26$shared_SNP_count))
}

ct28.random <- list()
rhos <- list()
pvals <- list()
for(i in seq_along(randoms)){
  ct28.random[[i]] <- total.ct28[randoms[[i]],]
  rhos[[i]] <- cor.test(ct28.random[[i]]$ct, 
                        ct28.random[[i]]$shared_SNP_count,
                        method="spearman",
                        exact = F)$estimate
  pvals[[i]] <- cor.test(ct28.random[[i]]$ct, 
                         ct28.random[[i]]$shared_SNP_count,
                         method="spearman",
                         exact = F)$p.value
}
rhos <- unlist(rhos)
pvals <- unlist(pvals)
mean(rhos) #0.1503527
mean(pvals) #0.1150407


