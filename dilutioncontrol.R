library(tidyverse)
library(rio)
library(here)

#ct = 23.6
ct23.dir <- "dilution_output_tables/"
ct23.files <- list.files(ct23.dir,pattern = "ivar")
ct23.files <- ct23.files[c(1:10)]
ct23.full <- paste0(ct23.dir,ct23.files)
ct23 <- list()
for(i in 1:length(ct23.full)){
  ct23[[i]] <- import(ct23.full[[i]])
}

ct23.all <- ct23

for(i in 1:length(ct23)){
  ct23[[i]] <- ct23[[i]]%>%
    select(POS,REF,ALT,ALT_FREQ,REF_DP,ALT_DP)%>%
    filter(POS == "3267" | POS == "5388" | POS == "6954" | POS == "11287"
           | POS == "21764" | POS == "21990" | POS == "23063" | POS == "23271"
           | POS == "23604" | POS == "23709" | POS == "24506" | POS == "24914"
           | POS == "27972" | POS == "28048" | POS == "28111" | POS == "28280"
           | POS == "28281" | POS == "28282" | POS == "28977") %>%
    distinct()
}
ct23.filenames <- paste0("dilution_controls/ct23/",
                         gsub(".ivar.tsv","_dilution",ct23.files),".csv")
for(i in 1:length(ct23)){
  write.csv(ct23[[i]],ct23.filenames[[i]])
}

for(i in seq_along(ct23.all)){
  ct23.all[[i]] <- ct23.all[[i]] %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP)%>%
    filter(ALT_FREQ>=.03,ALT_FREQ<=.97)%>%
    filter(TOTAL_DP >= 1000)%>%
    filter(POS != 6696)%>% #these are likely sequencing artifacts
    filter(POS!= 11074)%>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)%>%
    filter(POS != 187) %>%
    filter(POS != 1059) %>%
    filter(POS != 2094) %>%
    filter(POS != 3037) %>%
    filter(POS != 3130) %>%
    filter(POS != 6696) %>%
    filter(POS != 6990) %>%
    filter(POS != 8022) %>%
    filter(POS != 10323) %>%
    filter(POS != 10741) %>%
    filter(POS != 11074) %>%
    filter(POS != 13408) %>%
    filter(POS != 14786) %>%
    filter(POS != 19684) %>%
    filter(POS != 20148) %>%
    filter(POS != 21137) %>%
    filter(POS != 24034) %>%
    filter(POS != 24378) %>%
    filter(POS != 25563) %>%
    filter(POS != 26144) %>%
    filter(POS != 26461) %>%
    filter(POS != 26681) %>%
    filter(POS != 28077) %>%
    filter(POS != 28826) %>%
    filter(POS != 28854) %>%
    filter(POS != 29051) %>%
    filter(POS != 29700)
}

ct23.lengths <- vector(length = length(ct23.all))
for(i in seq_along(ct23.all)){
  ct23.lengths[[i]] <- length(ct23.all[[i]]$POS)
}
names(ct23.lengths) <- c("50_R1_23","10_R1_23","5_R1_23","2_R1_23","1_R1_23",
                         "50_R2_23","10_R2_23","5_R2_23","2_R2_23","1_R2_23")
ct23.r1 <- ct23.lengths[c(1:5)]
ct23.r2 <- ct23.lengths[c(6:10)]

plot(ct23.r1,ct23.r2)
cor.test(ct23.r1,ct23.r2, method = "pearson") #cor = 0.9994602 # p = 1.505e-05

dilution23 <- read.csv("dilution_controls/ct23_dilutions.csv")
colnames(dilution23)[[4]] <- "Spike-In %"
cor.test(dilution23$R1,dilution23$R2, method = "pearson") #cor = 0.994992 #p < 2.2e-16

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data=dilutiongroup, aes(x=R1,y=R2,col=`Spike-In %`)) +
  geom_point(size=4) +
  labs(x="R1 Allele Frequency",y="R2 Allele Frequency")+
  scale_color_discrete(breaks = c("1%","2%","5%","10%","50%"))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  ggtitle("Ct = 23.6")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))
ggsave("figs/ct23_dilution.png")

#Ct = 26
ct26.dir <- "dilution_output_tables/"
ct26.files <- list.files(ct26.dir,pattern = "ivar")
ct26.files <- ct26.files[c(11:20)]
ct26.full <- paste0(ct26.dir,ct26.files)
ct26 <- list()
for(i in 1:length(ct26.full)){
  ct26[[i]] <- import(ct26.full[[i]])
}
ct26.all <- ct26

for(i in 1:length(ct26)){
  ct26[[i]] <- ct26[[i]]%>%
    select(POS,REF,ALT,ALT_FREQ,REF_DP,ALT_DP)%>%
    filter(POS == "3267" | POS == "5388" | POS == "6954" | POS == "11287"
           | POS == "21764" | POS == "21990" | POS == "23063" | POS == "23271"
           | POS == "23604" | POS == "23709" | POS == "24506" | POS == "24914"
           | POS == "27972" | POS == "28048" | POS == "28111" | POS == "28280"
           | POS == "28281" | POS == "28282" | POS == "28977") %>%
    distinct()
}
ct26.filenames <- paste0("dilution_controls/ct26/",
                         gsub(".ivar.tsv","_dilution",ct26.files),".csv")
for(i in 1:length(ct26)){
  write.csv(ct26[[i]],ct26.filenames[[i]])
}

for(i in seq_along(ct26.all)){
  ct26.all[[i]] <- ct26.all[[i]] %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP)%>%
    filter(ALT_FREQ>=.03,ALT_FREQ<=.97)%>%
    filter(TOTAL_DP >= 1000)%>%
    filter(POS != 6696)%>% #these are likely sequencing artifacts
    filter(POS!= 11074)%>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)%>%
    filter(POS != 187) %>%
    filter(POS != 1059) %>%
    filter(POS != 2094) %>%
    filter(POS != 3037) %>%
    filter(POS != 3130) %>%
    filter(POS != 6696) %>%
    filter(POS != 6990) %>%
    filter(POS != 8022) %>%
    filter(POS != 10323) %>%
    filter(POS != 10741) %>%
    filter(POS != 11074) %>%
    filter(POS != 13408) %>%
    filter(POS != 14786) %>%
    filter(POS != 19684) %>%
    filter(POS != 20148) %>%
    filter(POS != 21137) %>%
    filter(POS != 24034) %>%
    filter(POS != 24378) %>%
    filter(POS != 25563) %>%
    filter(POS != 26144) %>%
    filter(POS != 26461) %>%
    filter(POS != 26681) %>%
    filter(POS != 28077) %>%
    filter(POS != 28826) %>%
    filter(POS != 28854) %>%
    filter(POS != 29051) %>%
    filter(POS != 29700)
}

ct26.lengths <- vector(length = length(ct26.all))
for(i in seq_along(ct26.all)){
  ct26.lengths[[i]] <- length(ct26.all[[i]]$POS)
}
names(ct26.lengths) <- c("50_R1_26","10_R1_26","5_R1_26","2_R1_26","1_R1_26",
                         "50_R2_26","10_R2_26","5_R2_26","2_R2_26","1_R2_26")
ct26.r1 <- ct26.lengths[c(1:5)]
ct26.r2 <- ct26.lengths[c(6:10)]

plot(ct26.r1,ct26.r2)
cor.test(ct26.r1,ct26.r2, method = "pearson") #cor = 0.9944083 #p = 0.0005015

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
dilution26 <- read.csv("dilution_controls/ct26_dilutions.csv")
colnames(dilution26)[[4]] <- "Spike-In %"
cor.test(dilution26$R1,dilution26$R2, method = "pearson") #cor = 0.9924263 #p < 2.2e-16


ggplot(data=dilution26, aes(x=R1,y=R2,col=`Spike-In %`)) +
  geom_point(size=4) +
  #geom_smooth(method = "lm")+
  labs(x="R1 Allele Frequency",y="R2 Allele Frequency")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  ggtitle("Ct = 26")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))
ggsave("figs/ct26_dilution.png")
  
  
#Ct = 28
ct28.dir <- "dilution_output_tables/"
ct28.files <- list.files(ct28.dir,pattern = "ivar")
ct28.files <- ct28.files[c(21:30)]
ct28.full <- paste0(ct28.dir,ct28.files) 
ct28 <- list()
for(i in 1:length(ct28.full)){
  ct28[[i]] <- import(ct28.full[[i]])
}
ct28.all <- ct28
for(i in 1:length(ct28)){
  ct28[[i]] <- ct28[[i]]%>%
    select(POS,REF,ALT,ALT_FREQ,REF_DP,ALT_DP)%>%
    filter(POS == "3267" | POS == "5388" | POS == "6954" | POS == "11287"
           | POS == "21764" | POS == "21990" | POS == "23063" | POS == "23271"
           | POS == "23604" | POS == "23709" | POS == "24506" | POS == "24914"
           | POS == "27972" | POS == "28048" | POS == "28111" | POS == "28280"
           | POS == "28281" | POS == "28282" | POS == "28977") %>%
    distinct()
}
ct28.filenames <- paste0("dilution_controls/ct28/",
                         gsub(".ivar.tsv","_dilution",ct28.files),".csv")
for(i in 1:length(ct28)){
  write.csv(ct28[[i]],ct28.filenames[[i]])
}

for(i in seq_along(ct28.all)){
  ct28.all[[i]] <- ct28.all[[i]] %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP)%>%
    filter(ALT_FREQ>=.03,ALT_FREQ<=.97)%>%
    filter(TOTAL_DP >= 1000)%>%
    filter(POS != 6696)%>% #these are likely sequencing artifacts
    filter(POS!= 11074)%>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)%>%
    filter(POS != 187) %>%
    filter(POS != 1059) %>%
    filter(POS != 2094) %>%
    filter(POS != 3037) %>%
    filter(POS != 3130) %>%
    filter(POS != 6696) %>%
    filter(POS != 6990) %>%
    filter(POS != 8022) %>%
    filter(POS != 10323) %>%
    filter(POS != 10741) %>%
    filter(POS != 11074) %>%
    filter(POS != 13408) %>%
    filter(POS != 14786) %>%
    filter(POS != 19684) %>%
    filter(POS != 20148) %>%
    filter(POS != 21137) %>%
    filter(POS != 24034) %>%
    filter(POS != 24378) %>%
    filter(POS != 25563) %>%
    filter(POS != 26144) %>%
    filter(POS != 26461) %>%
    filter(POS != 26681) %>%
    filter(POS != 28077) %>%
    filter(POS != 28826) %>%
    filter(POS != 28854) %>%
    filter(POS != 29051) %>%
    filter(POS != 29700)
}

ct28.lengths <- vector(length = length(ct28.all))
for(i in seq_along(ct28.all)){
  ct28.lengths[[i]] <- length(ct28.all[[i]]$POS)
}
names(ct28.lengths) <- c("50_R1_28","10_R1_28","5_R1_28","2_R1_28","1_R1_28",
                         "50_R2_28","10_R2_28","5_R2_28","2_R2_28","1_R2_28")
ct28.r1 <- ct28.lengths[c(1:5)]
ct28.r2 <- ct28.lengths[c(6:10)]

plot(ct28.r1,ct28.r2)
cor.test(ct28.r1,ct28.r2, method = "pearson") #0.946761

dilution28 <- read.csv("dilution_controls/ct28_dilutions.csv")
colnames(dilution28)[[4]] <- "Spike-In %"
cor.test(dilution28$R1,dilution28$R2, method = "pearson") #cor = 0.9667631 #p < 2.2e-16

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data=dilution28, aes(x=R1,y=R2,col=`Spike-In %`)) +
  geom_point(size=4) +
  #geom_smooth(method = "lm")+
  labs(x="R1 Allele Frequency",y="R2 Allele Frequency")+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")+
  ggtitle("Ct = 28")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.text = element_text(size=19),legend.title = element_text(size=22),
        axis.title.x = element_text(margin = margin(t=15)),
        axis.title.y = element_text(margin = margin(r=15)),
        plot.title = element_text(size = 22))+
  scale_color_manual(values = five.palette, breaks = c("1%","2%","5%","10%","50%"))
ggsave("figs/ct28_dilution.png")

#freq threshold estimates
dilutiongroups <- bind_rows(dilution23,dilution26,dilution28)
dilutiongroups <- dilutiongroups %>%
  group_split(`Spike-In %`)

cor.test(dilutiongroups$`50%`$R1, dilutiongroups$`50%`$R2, method = "pearson") 
#cor = 0.9340078 #p < 2.2e-16
cor.test(dilutiongroups$`10%`$R1, dilutiongroups$`10%`$R2, method = "pearson")
#cor = 0.8861263 #p < 2.2e-16
cor.test(dilutiongroups$`5%`$R1, dilutiongroups$`5%`$R2, method = "pearson")
#cor = 0.8389804 #p = 2.358e-15
cor.test(dilutiongroups$`2%`$R1, dilutiongroups$`2%`$R2, method = "pearson")
#cor = 0.519296 #p = 5.732e-05
cor.test(dilutiongroups$`1%`$R1, dilutiongroups$`1%`$R2, method = "pearson")
#cor = 0.5627796 #p = 9.443e-06




