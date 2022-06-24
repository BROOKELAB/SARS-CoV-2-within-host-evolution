library(tidyverse)
library(rio)
library(here)
#Ct = 26
ct26.dir <- "07_09_21_output_tables/"
ct26.files <- list.files(ct26.dir,pattern = "ivar")
ct26.files <- ct26.files[c(39:44,46:49)]
ct26.full <- paste0(ct26.dir,ct26.files)
ct26 <- list()
for(i in 1:length(ct26.full)){
  ct26[[i]] <- import(ct26.full[[i]])
}
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

#this .csv is a compilation of all csvs in ct26 labeled by dilution group
dilution26 <- read.csv("dilution_controls/ct26_dilutions.csv")

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
colnames(dilution26)[[4]] <- "Spike-In %"
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
ct28.dir <- "07_09_21_output_tables/"
ct28.files <- list.files(ct28.dir,pattern = "ivar")
ct28.files <- ct28.files[c(50:55,57:60)] #54 through 63
ct28.full <- paste0(ct28.dir,ct28.files) 
ct28 <- list()
for(i in 1:length(ct28.full)){
  ct28[[i]] <- import(ct28.full[[i]])
}
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

dilution28 <- read.csv("dilution_controls/ct28_dilutions.csv")
colnames(dilution28)[[4]] <- "Spike-In %"

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


#Ct = 23.6
dilutiongroup <- read.csv("dilution_controls/B117_dilutiongroups.csv")
colnames(dilutiongroup)[[4]] <- "Spike-In %"
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
 
coroverall <- cor.test(dilutiongroup$R1,dilutiongroup$R2)
#cor = 0.994992 #p < 2.2e-16

basedilutiongroup <- read.csv("COVID_data/dilution_controls/B117_dilutiongroups_base.csv")
cor1 <- cor.test(basedilutiongroup$R1_1,basedilutiongroup$R2_1)
#cor = 0.4986644 #p = 0.04928
cor2 <- cor.test(basedilutiongroup$R1_2,basedilutiongroup$R2_2)
#cor = -0.2405076 #p = 0.3696
cor5 <- cor.test(basedilutiongroup$R1_5,basedilutiongroup$R2_5)
#cor = 0.7210975 #p = 0.00162
cor10 <- cor.test(basedilutiongroup$R1_10,basedilutiongroup$R2_10)
#cor = 0.5348757 #p = 0.03278
cor50 <- cor.test(basedilutiongroup$R1_50,basedilutiongroup$R2_50)
#cor = 0.9260392 # p = 2.665e-07
sd1 <- sd(basedilutiongroup$R1_1,basedilutiongroup$R2_1)
#0.002518875
sd2 <- sd(basedilutiongroup$R1_2,basedilutiongroup$R2_2)
#0.003081759
sd5 <- sd(basedilutiongroup$R1_5,basedilutiongroup$R2_5)
#0.01128875
sd10 <- sd(basedilutiongroup$R1_10,basedilutiongroup$R2_10)
#0.01726524
sd50 <- sd(basedilutiongroup$R1_50,basedilutiongroup$R2_50)
#0.06513716
sdoverall <-sd(dilutiongroup$R1,dilutiongroup$R2)

plot(basedilutiongroup$R1_2,basedilutiongroup$R2_2)
group2 <- dilutiongroup[17:32,]
group1 <- dilutiongroup[1:16,]
plot1 <- plot(group1$R1,group1$R2)
abline(lm(group1$R2~group1$R1))
plot2 <- plot(group2$R1,group2$R2)

##ct = 26##
dilution26 <- read.csv("dilution_controls/ct26_dilutions.csv")
ggplot(data=dilution26, aes(x=R1,y=R2,col=Group)) +
  geom_point() +
  #geom_smooth(method = "lm")+
  labs(x="log10(R1 Allele Frequency)",y="log10(R2 Allele Frequency)")+
  scale_color_discrete(breaks = c("1%","2%","5%","10%","50%"))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")

##ct = 28##
dilution28 <- read.csv("dilution_controls/ct28_dilutions.csv")
ggplot(data=dilution28, aes(x=R1,y=R2,col=Group)) +
  geom_point() +
  #geom_smooth(method = "lm")+
  labs(x="log10(R1 Allele Frequency)",y="log10(R2 Allele Frequency)")+
  scale_color_discrete(breaks = c("1%","2%","5%","10%","50%"))+
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(trans = "log10")

