library(tidyverse)
library(rio)
library(here)

#ct vs total snps
naive.ct <- import("naive_ct.xlsx")
vax.ct <- import("vaccinated_ct.csv")
vax.ct <- vax.ct[,-1]
together.ct <-  full_join(naive.ct,vax.ct)
#together.ct.cut <- together.ct[which(together.ct$ct < 25),]

ggplot(data = together.ct, aes(x=ct,y=SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        axis.title.y = element_text(vjust = +1.8))+
  xlab("N Gene Ct")+
  ylab("iSNV Count")
ggsave("figs/total_ctSNPs.png")
cor.test(together.ct$ct,together.ct$SNP_count,method="spearman",
                               exact = F)
#rho = 0.6448757, p-value < 2.2e-16

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

#plot rho and p-value for different thresholds
cor.matrix <- as.data.frame(matrix(nrow = 9, ncol = 3))
colnames(cor.matrix) <- c("ct","rho","p-value")
cor.matrix$ct <- c(30,29,28,27,26,25,24,23,22)
total.ct.filter <- list()
for(i in 1:nrow(cor.matrix)){
  total.ct.filter[[i]] <- total.ct.shared[total.ct.shared$ct < cor.matrix$ct[i],]
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
total.ct30 <- total.ct.shared[total.ct.shared$ct < 30,]
total.ct28 <- total.ct.shared[total.ct.shared$ct < 28,]
total.ct26 <- total.ct.shared[total.ct.shared$ct < 26,]
total.ct24 <- total.ct.shared[total.ct.shared$ct < 24,]

ggplot(data = total.ct30, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,30),breaks = c(14,16,18,20,22,24,26,28,30))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 30")
ggsave("figs/plot30.png")

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

ggplot(data = total.ct24, aes(x = ct, y = shared_SNP_count))+
  geom_point(size = 4)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
        plot.title = element_text(size = 22))+
  scale_x_continuous(limits = c(14,24),breaks = c(14,16,18,20,22,24))+
  scale_y_continuous(limits = c(0,10), breaks = c(0,2,4,6,8,10))+
  xlab("Ct")+
  ylab("Shared iSNV Count")+
  ggtitle("Threshold = 24")
ggsave("figs/plot24.png")

