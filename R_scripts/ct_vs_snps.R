library(tidyverse)
library(rio)
library(here)

#correlation between Ct and total iSNVs
naive.ct <- import("naive_ct.xlsx")
summary(glm(total_SNPs ~ ct, naive.ct, family = poisson(link = "log")))
#coef = 0.4491 #p < 2e-16
ggplot(data = naive.ct, aes(x=ct,y=total_SNPs))+
  geom_point(size = 4)+
  xlab("Ct")+
  ylab("iSNV count")+
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+
  scale_x_continuous(limits = c(15,30),breaks = c(15,20,25,30))+
  theme_bw()+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=19),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/saliva_ctSNPs.png")

#correlation between Ct and shared iSNVs at different Ct thresholds
cor.matrix <- as.data.frame(matrix(nrow = 8, ncol = 3))
colnames(cor.matrix) <- c("ct","estimate","p-value")
cor.matrix$ct <- c(28,27,26,25,24,23,22,20)
filter.ct <- list()
for(i in 1:nrow(cor.matrix)){
  filter.ct[[i]] <- naive.ct[naive.ct$ct < cor.matrix$ct[i],]
  cor.matrix[i,2] <- glm(shared_SNPs ~ ct, filter.ct[[i]], 
                         family = poisson(link = "log"))$coefficients[2]
  cor.matrix[i,3] <- coef(summary(glm(shared_SNPs ~ ct, filter.ct[[i]], 
                         family = poisson(link = "log"))))[,4][2]
}

cor.matrix.gather <- cor.matrix %>%
  gather(key = "stat", value = "value", - ct)

ggplot(data = cor.matrix.gather, aes(x = ct, y = value, shape = stat))+
  geom_point(cex = 3)+
  geom_line()+
  geom_abline(slope = 0, intercept = 0.05, linetype = "dashed", color = "grey50")+
  ylim(c(-.2,1))+
  scale_x_reverse()+
  scale_shape_manual(values = c(19,21),
                     labels = c("Coefficient", "P-value"),
                     guide = guide_legend(title = "Statistic"))+
  xlab("Ct")+
  ylab("Value")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19))
ggsave("figs/correlation_stats.png")

#plots for different thresholds
ct28 <- naive.ct[naive.ct$ct < 28,]
ct26 <- naive.ct[naive.ct$ct < 26,]
ct24 <- naive.ct[naive.ct$ct < 24,]
ct22 <- naive.ct[naive.ct$ct < 22,]

ct.list <- list("28"= ct28,"26"=ct26,"24"=ct24,"22"=ct22)

ct.plots <- list()
for(i in seq_along(ct.list)){
  ct.plots[[i]] <- ggplot(data = ct.list[[i]], aes(x = ct, y = shared_SNPs))+
    geom_point(size = 4)+
    #geom_smooth(method = "glm", method.args = list(family = "poisson"))+
    scale_x_continuous(limits = c(15,29),breaks = c(16,18,20,22,24,26,28))+
    scale_y_continuous(limits = c(0,12), breaks = c(0,2,4,6,8,10,12))+
    xlab("Ct")+
    ylab("Shared iSNV count")+
    ggtitle(paste("Threshold =",names(ct.list)[[i]]))+
    theme_bw()+
    theme(axis.title = element_text(size=22),axis.text = element_text(size=19),
          plot.title = element_text(size = 22))
}

#correlations for vaccinees
vax.ct <- import("vaccinated_ct.csv")
summary(glm(shared_SNPs ~ ct, vax.ct, family = poisson(link = "log")))
#coef = -0.05731 #p = 0.00806

ggplot(data = vax.ct, aes(x=ct,y=shared_SNPs))+
  geom_point(size = 4)+
  xlab("Ct")+
  ylab("iSNV count")+
  ggtitle("Vaccinee samples")+
  geom_smooth(method = "glm", method.args = list(family = "poisson"))+
  scale_x_continuous(limits = c(15,30),breaks = c(15,20,25,30))+
  theme_bw()+
  theme(axis.title = element_text(size=22),
        axis.text = element_text(size=19),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))
ggsave("figs/vax_ctSNPs.png")


#correlations for nasal samples
load("nasal_snpcounts.RData")
nasal.ct <- import("all_nasal_user_info.xlsx")%>%
  select(user_ID, day_of_infection, cn)%>%
  mutate("shared_SNPs" = nasal.snpcounts$shared_SNPs)%>%
  mutate("total_SNPs" = nasal.snpcounts$total_SNPs)

summary(glm(shared_SNPs ~ cn, nasal.ct, family = poisson(link = "log")))
#coef = 0.005064 #p = 0.927

ggplot(data = nasal.ct, aes(x=cn,y=shared_SNPs))+
  geom_point(size = 4)+
  scale_x_continuous(limits = c(13,30),breaks = c(15,20,25,30))+
  scale_y_continuous(limits = c(0,3), breaks = c(0,1,2,3))+
  xlab("CN")+
  ylab("Shared iSNV count")+
  ggtitle("Nasal samples")+
  theme_bw()+
  theme(axis.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t =10)),
        axis.title.y = element_text(margin = margin(r =15)),
        axis.text = element_text(size=19),
        plot.title = element_text(size = 22))
ggsave("figs/nasal_ctSNPs.png")



