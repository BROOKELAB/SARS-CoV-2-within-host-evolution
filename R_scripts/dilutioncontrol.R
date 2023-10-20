library(tidyverse)
library(rio)
library(Nmisc)
library(here)

##Ct = 23.6##
ct23.dir <- "dilution_output_tables/"
ct23.files <- list.files(ct23.dir,pattern = "ivar")
ct23.files <- ct23.files[c(1:10)]
ct23.full <- paste0(ct23.dir,ct23.files)
ct23 <- list()
for(i in 1:length(ct23.full)){
  ct23[[i]] <- import(ct23.full[[i]])
}
names(ct23) <- c("50%A","10%A","5%A","2%A","1%A","50%B","10%B","5%B","2%B","1%B")

ct23.all <- ct23

snp.list <- ct23$`50%A` %>%
  mutate(ct23$`50%A`, "SNP" = paste0(REF,POS,ALT)) %>%
  filter(SNP == "C3267T" | SNP == "C5388A" | SNP == "T6954C" | SNP == "G11287-TCTGGTTTT"
         | SNP == "A21764-TACATG" | SNP == "T21990-TTA" | SNP == "A23063T" | SNP == "C23271A"
         | SNP == "C23604A" | SNP == "C23709T" | SNP == "T24506G" | SNP == "G24914C"
         | SNP == "C27972T" | SNP == "G28048T" | SNP == "A28111G" | SNP == "G28280C"
         | SNP == "A28281T" | SNP == "T28282A" | SNP == "C28977T") %>%
  select(SNP)%>%
  distinct()
save(snp.list, file = "alpha_SNPs.RData")

alpha.filter <- function(control){
  control <- control %>% 
    mutate("SNP" = paste0(REF,POS,ALT)) %>%
    select(SNP, ALT_FREQ) %>%
    distinct()
  
  table <- left_join(snp.list, control)
  table[is.na(table)] <- 0
  return(table)
}

ct23 <- lapply(ct23, alpha.filter)

ct23.R1 <- keep_at(ct23, ends_with("A"))
names(ct23.R1) <- gsub("A","",names(ct23.R1))
ct23.R2 <- keep_at(ct23, ends_with("B"))
names(ct23.R2) <- gsub("B","",names(ct23.R2))

for(i in seq_along(ct23.R1)){
  ct23.R1[[i]] <- ct23.R1[[i]] %>%
    mutate("DIL" = names(ct23.R1)[[i]]) %>%
    relocate(DIL, .before = SNP)
  ct23.R2[[i]] <- ct23.R2[[i]] %>%
    mutate("DIL" = names(ct23.R2)[[i]]) %>%
    relocate(DIL, .before = SNP)
}

ct23.table <- map2(ct23.R1, ct23.R2, full_join, by = c("DIL","SNP"))
ct23.table <- bind_rows(ct23.table)
colnames(ct23.table)[c(3:4)] <- c("R1", "R2")
write.csv(ct23.table, file = "dilution_controls/ct23_dilutions.csv")

cor.test(ct23.table$R1, ct23.table$R2, method = "pearson") #cor = 0.9958906 #p < 2.2e-16

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")

#plot
ggplot(data=ct23.table, aes(x=R1,y=R2,col = DIL)) +
  geom_point(size=4) +
  labs(x="R1 allele frequency",y="R2 allele frequency")+
  scale_x_continuous(limits = c(0,0.51),
                     breaks = c(0,0.1, 0.2, 0.3,0.4, 0.5))+
  scale_y_continuous(limits = c(0,0.51),
                     breaks = c(0,0.1, 0.2, 0.3,0.4, 0.5))+
  ggtitle("Ct = 23.6")+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  guides(color = guide_legend(title = "% Spike-in"))
ggsave("figs/ct23_dilution.png")


##Ct = 26##
ct26.dir <- "dilution_output_tables/"
ct26.files <- list.files(ct26.dir,pattern = "ivar")
ct26.files <- ct26.files[c(11:20)]
ct26.full <- paste0(ct26.dir,ct26.files)
ct26 <- list()
for(i in 1:length(ct26.full)){
  ct26[[i]] <- import(ct26.full[[i]])
}
names(ct26) <- c("50%A","10%A","5%A","2%A","1%A","50%B","10%B","5%B","2%B","1%B")

ct26.all <- ct26

load("alpha_SNPs.RData")
alpha.filter <- function(control){
  control <- control %>% 
    mutate("SNP" = paste0(REF,POS,ALT)) %>%
    select(SNP, ALT_FREQ) %>%
    distinct()
  
  table <- left_join(snp.list, control)
  table[is.na(table)] <- 0
  return(table)
}

ct26 <- lapply(ct26, alpha.filter)

ct26.R1 <- keep_at(ct26, ends_with("A"))
names(ct26.R1) <- gsub("A","",names(ct26.R1))
ct26.R2 <- keep_at(ct26, ends_with("B"))
names(ct26.R2) <- gsub("B","",names(ct26.R2))

for(i in seq_along(ct26.R1)){
  ct26.R1[[i]] <- ct26.R1[[i]] %>%
    mutate("DIL" = names(ct26.R1)[[i]]) %>%
    relocate(DIL, .before = SNP)
  ct26.R2[[i]] <- ct26.R2[[i]] %>%
    mutate("DIL" = names(ct26.R2)[[i]]) %>%
    relocate(DIL, .before = SNP)
}

ct26.table <- map2(ct26.R1, ct26.R2, full_join, by = c("DIL","SNP"))
ct26.table <- bind_rows(ct26.table)
colnames(ct26.table)[c(3:4)] <- c("R1", "R2")
write.csv(ct26.table, file = "dilution_controls/ct26_dilutions.csv")

cor.test(ct26.table$R1, ct26.table$R2, method = "pearson") #cor = 0.9924263 #p < 2.2e-16

#plot
five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data=ct26.table, aes(x=R1,y=R2,col = DIL)) +
  geom_point(size=4) +
  labs(x="R1 allele frequency",y="R2 allele frequency")+
  scale_x_continuous(limits = c(0,0.8),
                     breaks = c(0, 0.2, 0.4, 0.6,0.8))+
  scale_y_continuous(limits = c(0,0.8),
                     breaks = c(0, 0.2, 0.4, 0.6,0.8))+
  ggtitle("Ct = 26")+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  guides(color = guide_legend(title = "% Spike-in"))
ggsave("figs/ct26_dilution.png")

##Ct = 28##
ct28.dir <- "dilution_output_tables/"
ct28.files <- list.files(ct28.dir,pattern = "ivar")
ct28.files <- ct28.files[c(21:30)]
ct28.full <- paste0(ct28.dir,ct28.files)
ct28 <- list()
for(i in 1:length(ct28.full)){
  ct28[[i]] <- import(ct28.full[[i]])
}
names(ct28) <- c("50%A","10%A","5%A","2%A","1%A","50%B","10%B","5%B","2%B","1%B")

ct28.all <- ct28

load("alpha_SNPs.RData")
alpha.filter <- function(control){
  control <- control %>% 
    mutate("SNP" = paste0(REF,POS,ALT)) %>%
    select(SNP, ALT_FREQ) %>%
    distinct()
  
  table <- left_join(snp.list, control)
  table[is.na(table)] <- 0
  return(table)
}

ct28 <- lapply(ct28, alpha.filter)

ct28.R1 <- keep_at(ct28, ends_with("A"))
names(ct28.R1) <- gsub("A","",names(ct28.R1))
ct28.R2 <- keep_at(ct28, ends_with("B"))
names(ct28.R2) <- gsub("B","",names(ct28.R2))

for(i in seq_along(ct28.R1)){
  ct28.R1[[i]] <- ct28.R1[[i]] %>%
    mutate("DIL" = names(ct28.R1)[[i]]) %>%
    relocate(DIL, .before = SNP)
  ct28.R2[[i]] <- ct28.R2[[i]] %>%
    mutate("DIL" = names(ct28.R2)[[i]]) %>%
    relocate(DIL, .before = SNP)
}

ct28.table <- map2(ct28.R1, ct28.R2, full_join, by = c("DIL","SNP"))
ct28.table <- bind_rows(ct28.table)
colnames(ct28.table)[c(3:4)] <- c("R1", "R2")
write.csv(ct28.table, file = "dilution_controls/ct28_dilutions.csv")

cor.test(ct28.table$R1, ct28.table$R2, method = "pearson") #cor = 0.9667631 #p < 2.2e-16

#plot
five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data=ct28.table, aes(x=R1,y=R2,col = DIL)) +
  geom_point(size=4) +
  labs(x="R1 allele frequency",y="R2 allele frequency")+
  scale_x_continuous(limits = c(0,0.51),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  scale_y_continuous(limits = c(0,0.51),
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))+
  ggtitle("Ct = 28")+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  guides(color = guide_legend(title = "% Spike-in"))
ggsave("figs/ct28_dilution.png")

#spike-in vs correlation
alpha.groups <- bind_rows(ct23.table,ct26.table,ct28.table)
alpha.groups <- alpha.groups %>%
  group_by(DIL, .add = T)%>%
  group_split()
names(alpha.groups) <- c("1%","10%","2%","5%","50%")

alpha.cor <- as.data.frame(matrix(data = NA, nrow = 5, ncol = 2))
colnames(alpha.cor) <- c("% Spike-in","Correlation")
alpha.cor[,1] <- c(50,10,5,2,1)
alpha.cor[1,2] <- cor.test(alpha.groups$`50%`$R1, alpha.groups$`50%`$R2, method = "pearson")$estimate
alpha.cor[2,2] <- cor.test(alpha.groups$`10%`$R1, alpha.groups$`10%`$R2, method = "pearson")$estimate
alpha.cor[3,2] <- cor.test(alpha.groups$`5%`$R1, alpha.groups$`5%`$R2, method = "pearson")$estimate
alpha.cor[4,2] <- cor.test(alpha.groups$`2%`$R1, alpha.groups$`2%`$R2, method = "pearson")$estimate
alpha.cor[5,2] <- cor.test(alpha.groups$`1%`$R1, alpha.groups$`1%`$R2, method = "pearson")$estimate

ggplot(data = alpha.cor, aes(x = `% Spike-in` , y = Correlation))+
  geom_point(cex = 4)+
  geom_line()+
  scale_x_continuous(breaks = c(1,2,5,10,50))+
  scale_y_continuous(limits = c(0,1), breaks = c(0,0.25,0.50,0.75,1))+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/dilution_correlation.png")

#correlations for non-alpha iSNVs
other.filter <- function(ctall){
  ctall <- mutate(ctall,"SNP" = paste0(REF,POS,ALT))
  ctall <- ctall[-which(ctall$SNP %in% snp.list$SNP),]
  ctall <- ctall[-which(ctall$POS %in% c(6696,11074,15965,29051,78,187,635,1059,2094,3037,
                                         3130,6696,6990,8022,10323,10741,11074,11567,13408,
                                         14786,19684,20148,21137,24034,24378,25563,26144,
                                         26461,26681,27964,28077,28253,28262,28281,28472,
                                         28826,28854,29051,29700,29760)),]
  ctall <- ctall %>%
    select(SNP, ALT_FREQ, TOTAL_DP) %>%
    distinct()
  return(ctall)
}
ct23.other <- lapply(ct23.all, other.filter)
names(ct23.other) <- names(ct23.all)
ct26.other <- lapply(ct26.all, other.filter)
names(ct26.other) <- names(ct26.all)
ct28.other <- lapply(ct28.all, other.filter)
names(ct28.other) <- names(ct28.all)

other.table <- function(ct.other){
  R1 <- keep_at(ct.other, ends_with("A"))
  names(R1) <- gsub("A","",names(R1))
  R2 <- keep_at(ct.other, ends_with("B"))
  names(R2) <- gsub("B","",names(R2))
  
  for(i in seq_along(R1)){
    R1[[i]] <- R1[[i]] %>%
      mutate("DIL" = names(R1)[[i]]) %>%
      relocate(DIL, .before = SNP)
    R2[[i]] <- R2[[i]] %>%
      mutate("DIL" = names(R2)[[i]]) %>%
      relocate(DIL, .before = SNP)
  }
  
  table <- map2_dfr(R1, R2, full_join, by = c("DIL","SNP"))
  colnames(table)[c(3:6)] <- c("R1_FREQ", "R1_DP", "R2_FREQ","R2_DP")
  table[is.na(table)] <- 0
  
  table <- table %>%
    filter(R1_DP >= 1000 | R2_DP >= 1000)%>%
    select(-R1_DP, - R2_DP, -DIL)
  
  return(table)
}

ct.other.list <- list(ct23.other,ct26.other,ct28.other)

other.tables <- lapply(ct.other.list, other.table)
other.table <- bind_rows(other.tables)
cor.test(other.table$R1_FREQ, other.table$R2_FREQ, method = "p")
#cor = 0.996902  #p-value < 2.2e-16

other.10 <- other.table %>%
  filter(R1_FREQ >= 0.10 | R2_FREQ >= 0.10)
other.5 <- other.table %>%
  filter(R1_FREQ >= 0.05 | R2_FREQ >= 0.05)
other.3 <- other.table %>%
  filter(R1_FREQ >= 0.03 | R2_FREQ >= 0.03)
other.1 <- other.table %>%
  filter(R1_FREQ >= 0.01 | R2_FREQ >= 0.01)

#correlations of non-alpha iSNVs
ggplot(data = other.1, aes(x=R1_FREQ, y = R2_FREQ))+
  geom_point()+
  xlab("R1 frequency")+
  ylab("R2 freqeuncy")+
  ggtitle("Cutoff = 0.01")+
  theme_bw()+
  theme(plot.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/nonalpha_01.png")
cor.test(other.1$R1_FREQ, other.1$R2_FREQ, method = "p")
#cor = 0.996902 #pval < 2.2e-16

ggplot(data = other.3, aes(x=R1_FREQ, y = R2_FREQ))+
  geom_point()+
  xlab("R1 frequency")+
  ylab("R2 freqeuncy")+
  ggtitle("Cutoff = 0.03")+
  theme_bw()+
  theme(plot.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/nonalpha_03.png")
cor.test(other.3$R1_FREQ, other.3$R2_FREQ, method = "p") 
#cor = 0.9974495 #p < 2.2e-16

ggplot(data = other.5, aes(x=R1_FREQ, y = R2_FREQ))+
  geom_point()+
  xlab("R1 frequency")+
  ylab("R2 freqeuncy")+
  ggtitle("Cutoff = 0.05")+
  theme_bw()+
  theme(plot.title = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/nonalpha_05.png")
cor.test(other.5$R1_FREQ, other.5$R2_FREQ, method = "p")
#cor = 0.9967089 #p < 2.2e-16

#missing iSNVS at various frequencies
freq.cutoff <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 2))
colnames(freq.cutoff) <- c("threshold","missing")
freq.cutoff[,1] <- c(.01,.03,.05,.1)
freq.cutoff[1,2] <- length(which(other.1$R1_FREQ == 0 | other.1$R2_FREQ == 0))
freq.cutoff[2,2] <- length(which(other.3$R1_FREQ == 0 | other.3$R2_FREQ == 0))
freq.cutoff[3,2] <- length(which(other.5$R1_FREQ == 0 | other.5$R2_FREQ == 0)) 
freq.cutoff[4,2] <- length(which(other.10$R1_FREQ == 0 | other.10$R2_FREQ == 0))

ggplot(data = freq.cutoff, aes(x = threshold, y = missing))+
  geom_point(cex = 4)+
  geom_line()+
  scale_x_continuous(breaks = c(0.01,0.03,0.05,0.1))+
  scale_y_continuous(limits = c(0,3000), breaks = c(0,1000,2000,3000))+
  xlab("iSNV frequency threshold")+
  ylab("Missing iSNVs")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)))
ggsave("figs/freq_threshold.png")






