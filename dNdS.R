library(tidyverse)
library(rio)
library(gdata)
library(here)

#naives
naive.dirs <- list.dirs("naive_output_tables")[-c(1:3)]
freq.read <- function(dir){
  freq.full <- paste0(dir,"/",
                            list.files(dir,pattern = ".ivar"))
  freq.files <- list()
  for(i in seq_along(freq.full)){
    freq.files[[i]] <- import(freq.full[[i]])
  }
  return(freq.files)
}
naive.freq.files <- lapply(naive.dirs,freq.read)
names(naive.freq.files) <- gsub("naive_output_tables/","",naive.dirs)
for(i in seq_along(naive.freq.files)){
  naive.freq.files[[i]] <- lapply(naive.freq.files[[i]],select,
                             POS,REF,ALT,ALT_DP,ALT_FREQ,REF_DP)
  naive.freq.files[[i]] <- lapply(naive.freq.files[[i]],distinct)
}

snp.read <- function(dir){
  snp.full <- paste0(dir,"/",list.files(dir,pattern = "snp"))
  snp.files <- list()
  for(i in seq_along(snp.full)){
    snp.files[[i]] <- import(snp.full[[i]])
  }
  return(snp.files)
}
naive.snp.files <- lapply(naive.dirs,snp.read)
names(naive.snp.files) <- gsub("naive_output_tables/","",naive.dirs)
for(i in seq_along(naive.snp.files)){
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],dplyr::rename,"GENE" = `ANN[*].GENE`)
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],dplyr::rename,"IMP" = `ANN[*].IMPACT`)
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],dplyr::rename,"EFF" =  `EFF[*].EFFECT`)
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],dplyr::rename,"AA" = `EFF[*].AA`)
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],filter,IMP != "MODIFIER")
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],filter,IMP != "MODIFIER,MODIFIER")
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],distinct)
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],mutate,"EFF.CUT" = strsplit(EFF,","))
}
for(i in seq_along(naive.snp.files)){
  for(j in seq_along(naive.snp.files[[i]])){
    for(k in seq_along(naive.snp.files[[i]][[j]]$EFF)){
      naive.snp.files[[i]][[j]]$EFF[[k]] <- naive.snp.files[[i]][[j]]$EFF.CUT[[k]][[1]]
    }
  }
}
for(i in seq_along(naive.snp.files)){
  naive.snp.files[[i]] <- lapply(naive.snp.files[[i]],select,POS,REF,ALT,GENE,IMP,EFF,AA)
}

#need to standardize notation for indels between freqs and snps
for(i in seq_along(naive.freq.files)){
  for(j in seq_along(naive.freq.files[[i]])){
    for(k in seq_along(naive.freq.files[[i]][[j]]$ALT)){
      if(substr(naive.freq.files[[i]][[j]]$ALT[[k]],0,1)=="+"){
        naive.freq.files[[i]][[j]]$ALT[[k]] <- paste0(naive.freq.files[[i]][[j]]$REF[[k]], 
                                                 substr(naive.freq.files[[i]][[j]]$ALT[[k]],2,nchar(naive.freq.files[[i]][[j]]$ALT[[k]])))
      }
      if(substr(naive.freq.files[[i]][[j]]$ALT[[k]],0,1)=="-"){
        naive.freq.files[[i]][[j]]$REF[[k]] <- paste0(naive.freq.files[[i]][[j]]$REF[[k]],
                                                 substr(naive.freq.files[[i]][[j]]$ALT[[k]],2,nchar(naive.freq.files[[i]][[j]]$ALT[[k]])))
        naive.freq.files[[i]][[j]]$ALT[[k]] <- substr(naive.freq.files[[i]][[j]]$REF[[k]],0,1)
      }
    }
  }
}

#merge freq info with snp annotations
merge <- function(freq,snp){
  merged <- list()
  for(i in seq_along(freq)){
    merged[[i]] <- inner_join(freq[[i]],snp[[i]])
    merged[[i]] <- merged[[i]]%>%
      filter(ALT_FREQ >= .03)%>%
      filter(ALT_DP + REF_DP >=1000)%>%
      filter(POS != 6696)%>%
      filter(POS != 11074)%>%
      filter(POS != 29051)%>%
      distinct()
  }
  return(merged)
}
naives <- list()
for(i in seq_along(naive.freq.files)){
  naives[[i]] <- merge(naive.freq.files[[i]],naive.snp.files[[i]])
}
names(naives) <- names(naive.freq.files)
save(naives,file = "naives.RData")

#calculate dN/dS for each user
load("naives.RData")
calc.dNdS <- function(user){
  eff <- list()
  for(i in 1:length(user)){
    eff[[i]] <- user[[i]]$EFF
  }
  syn <- list()
  nonsyn <- list()
  dNdS <- list()
  for(j in 1:length(eff)){
    syn[[j]] <- (length(grep(pattern = "synonymous_variant",eff[[j]])))/9803 #9803 estimates number of synonymous sites
    nonsyn[[j]] <- (length(eff[[j]]) - syn[[j]])/19606 #19606 estimates number of nonsynonymous sites
    dNdS[[j]] <- (nonsyn[[j]])/(syn[[j]])
  }
  return(dNdS)
}
naive.dnds <- lapply(naives,calc.dNdS)
naive.names <- gsub("user_","",names(naives))

#combine all dNdS info into a dataframe
for(i in seq_along(naive.dnds)){
  length(naive.dnds[[i]]) <- 10
}
naive.dnds.table <- as.data.frame(rlist::list.cbind(naive.dnds))
colnames(naive.dnds.table) <- naive.names

#plot dNdS vs day
naive.dnds.table <- naive.dnds.table %>%
  mutate("day"=as.numeric(c(1:10)))%>%
  gather(key = "Participant ID", value = "dN/dS",-day)

#remove null/na/inf values
for(i in 1:length(naive.dnds.table$`dN/dS`)){
  if(is.null(naive.dnds.table$`dN/dS`[[i]])){
    naive.dnds.table$`dN/dS`[[i]] <- NA
  }
  if(is.nan(naive.dnds.table$`dN/dS`[[i]])){
    naive.dnds.table$`dN/dS`[[i]] <- NA
  }
  if(is.infinite(naive.dnds.table$`dN/dS`[[i]])){
    naive.dnds.table$`dN/dS`[[i]] <- NA
  }
}

naive.dnds.table$`dN/dS` <- as.numeric(naive.dnds.table$`dN/dS`)
naive.dnds.table$day <- as.numeric(naive.dnds.table$day)
naive.dnds.lm <- lm(`dN/dS`~day,data=naive.dnds.table)
summary(naive.dnds.lm)
#p-value: 0.03871, Adjusted R-squared: 0.0226 
naive.dnds.table$day <- as.factor(naive.dnds.table$day)
ggplot(data = naive.dnds.table,aes(x=day,y=`dN/dS`)) +
  geom_point(cex=4)+
  xlab("Day")+
  ggtitle("Unvaccinated")+
  geom_abline(intercept=1.8711177,slope=0.0633897, color= "grey50")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/naive_dnds_lm.png")

#vax
vax.dirs <- list.dirs("vaccinated_output_tables")[-c(1,10)]
freq.read <- function(dir){
  freq.full <- paste0(dir,"/",
                      list.files(dir,pattern = ".ivar"))
  freq.files <- list()
  for(i in seq_along(freq.full)){
    freq.files[[i]] <- import(freq.full[[i]])
  }
  return(freq.files)
}
vax.freq.files <- lapply(vax.dirs,freq.read)
names(vax.freq.files) <- gsub("vaccinated_output_tables/","",vax.dirs)
for(i in seq_along(vax.freq.files)){
  vax.freq.files[[i]] <- lapply(vax.freq.files[[i]],select,
                                  POS,REF,ALT,ALT_DP,ALT_FREQ,REF_DP)
  vax.freq.files[[i]] <- lapply(vax.freq.files[[i]],distinct)
}

snp.read <- function(dir){
  snp.full <- paste0(dir,"/",list.files(dir,pattern = "snp"))
  snp.files <- list()
  for(i in seq_along(snp.full)){
    snp.files[[i]] <- import(snp.full[[i]])
  }
  return(snp.files)
}
vax.snp.files <- lapply(vax.dirs,snp.read)
names(vax.snp.files) <- gsub("vaccinated_output_tables/","",vax.dirs)
for(i in seq_along(vax.snp.files)){
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],dplyr::rename,"GENE" = `ANN[*].GENE`)
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],dplyr::rename,"IMP" = `ANN[*].IMPACT`)
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],dplyr::rename,"EFF" =  `EFF[*].EFFECT`)
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],dplyr::rename,"AA" = `EFF[*].AA`)
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],filter,IMP != "MODIFIER")
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],filter,IMP != "MODIFIER,MODIFIER")
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],distinct)
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],mutate,"EFF.CUT" = strsplit(EFF,","))
}
for(i in seq_along(vax.snp.files)){
  for(j in seq_along(vax.snp.files[[i]])){
    for(k in seq_along(vax.snp.files[[i]][[j]]$EFF)){
      vax.snp.files[[i]][[j]]$EFF[[k]] <- vax.snp.files[[i]][[j]]$EFF.CUT[[k]][[1]]
    }
  }
}
for(i in seq_along(vax.snp.files)){
  vax.snp.files[[i]] <- lapply(vax.snp.files[[i]],select,POS,REF,ALT,GENE,IMP,EFF,AA)
}

#need to standardize notation for indels between freqs and snps
for(i in seq_along(vax.freq.files)){
  for(j in seq_along(vax.freq.files[[i]])){
    for(k in seq_along(vax.freq.files[[i]][[j]]$ALT)){
      if(substr(vax.freq.files[[i]][[j]]$ALT[[k]],0,1)=="+"){
        vax.freq.files[[i]][[j]]$ALT[[k]] <- paste0(vax.freq.files[[i]][[j]]$REF[[k]], 
                                                      substr(vax.freq.files[[i]][[j]]$ALT[[k]],2,nchar(vax.freq.files[[i]][[j]]$ALT[[k]])))
      }
      if(substr(vax.freq.files[[i]][[j]]$ALT[[k]],0,1)=="-"){
        vax.freq.files[[i]][[j]]$REF[[k]] <- paste0(vax.freq.files[[i]][[j]]$REF[[k]],
                                                      substr(vax.freq.files[[i]][[j]]$ALT[[k]],2,nchar(vax.freq.files[[i]][[j]]$ALT[[k]])))
        vax.freq.files[[i]][[j]]$ALT[[k]] <- substr(vax.freq.files[[i]][[j]]$REF[[k]],0,1)
      }
    }
  }
}

#merge freq info with snp annotations
merge <- function(freq,snp){
  merged <- list()
  for(i in seq_along(freq)){
    merged[[i]] <- inner_join(freq[[i]],snp[[i]])
    merged[[i]] <- merged[[i]]%>%
      filter(ALT_FREQ >= .03)%>%
      filter(ALT_DP + REF_DP >=1000)%>%
      filter(POS != 6696)%>%
      filter(POS != 11074)%>%
      filter(POS != 29051)%>%
      distinct()
  }
  return(merged)
}
vaccinateds <- list()
for(i in seq_along(vax.freq.files)){
  vaccinateds[[i]] <- merge(vax.freq.files[[i]],vax.snp.files[[i]])
}
names(vaccinateds) <- names(vax.freq.files)
save(vaccinateds,file = "vaccinateds.RData")

#calculate dN/dS for each user
load("vaccinateds.RData")
calc.dNdS <- function(user){
  eff <- list()
  for(i in 1:length(user)){
    eff[[i]] <- user[[i]]$EFF
  }
  syn <- list()
  nonsyn <- list()
  dNdS <- list()
  for(j in 1:length(eff)){
    syn[[j]] <- (length(grep(pattern = "synonymous_variant",eff[[j]])))/9803 #9803 estimates number of synonymous sites
    nonsyn[[j]] <- (length(eff[[j]]) - syn[[j]])/19606 #19606 estimates number of nonsynonymous sites
    dNdS[[j]] <- (nonsyn[[j]])/(syn[[j]])
  }
  return(dNdS)
}
vax.dnds <- lapply(vaccinateds,calc.dNdS)
vax.names <- gsub("dir_","",names(vaccinateds))

#combine all dNdS info into a dataframe
for(i in seq_along(vax.dnds)){
  length(vax.dnds[[i]]) <- 8
}
vax.dnds.table <- as.data.frame(rlist::list.cbind(vax.dnds))
colnames(vax.dnds.table) <- vax.names

#plot dNdS vs day
vax.dnds.table <- vax.dnds.table %>%
  mutate("day"=as.numeric(c(1:8)))%>%
  gather(key = "Participant ID", value = "dN/dS",-day)

#remove null/na/inf values
for(i in 1:length(vax.dnds.table$`dN/dS`)){
  if(is.null(vax.dnds.table$`dN/dS`[[i]])){
    vax.dnds.table$`dN/dS`[[i]] <- NA
  }
  if(is.nan(vax.dnds.table$`dN/dS`[[i]])){
    vax.dnds.table$`dN/dS`[[i]] <- NA
  }
  if(is.infinite(vax.dnds.table$`dN/dS`[[i]])){
    vax.dnds.table$`dN/dS`[[i]] <- NA
  }
}

vax.dnds.table$`dN/dS` <- as.numeric(vax.dnds.table$`dN/dS`)
vax.dnds.table$day <- as.numeric(vax.dnds.table$day)
vax.dnds.lm <- lm(`dN/dS`~day,data=vax.dnds.table)
summary(vax.dnds.lm)
#p-value: 0.5019, Adjusted R-squared: -0.01246 
vax.dnds.table$day <- as.factor(vax.dnds.table$day)
ggplot(data = vax.dnds.table,aes(x=day,y=`dN/dS`)) +
  geom_point(cex=4)+
  xlab("Day")+
  ggtitle("Vaccinated")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19),
        legend.position = "none", plot.title = element_text(size = 22))
ggsave("figs/vax_dnds_lm.png")

#compare naive and vax dnds values
all.dnds.table <- rbind(naive.dnds.table,vax.dnds.table)
all.dnds.table <- mutate(all.dnds.table,"Status"=NA)

for(i in 1:length(all.dnds.table$`Participant ID`)){
  if(all.dnds.table$`Participant ID`[[i]] %in% naive.names){
    all.dnds.table$Status[[i]] <- "Naive"
  } else {
    all.dnds.table$Status[[i]] <- "Vaccinated"
  }
}
all.dnds.table$`Participant ID` <- as.character(all.dnds.table$`Participant ID`)
all.dnds.table$`dN/dS` <- as.numeric(all.dnds.table$`dN/dS`)

only.naive <- all.dnds.table[c(1:200),]
only.vax <- all.dnds.table[c(201:296),]
mean.only.naive <- mean(na.exclude(only.naive$`dN/dS`)) #2.150727
mean.only.vax <- mean(na.exclude(only.vax$`dN/dS`)) #2.811289
t.test(only.naive$`dN/dS`,only.vax$`dN/dS`)
#p-value = 0.0002031

ggplot(data=all.dnds.table,aes(x=`Participant ID`,y=`dN/dS`,fill=Status,color=Status))+
  geom_dotplot(binaxis = "y",binwidth = .19,stackdir = "center")+
  xlab("Participant")+
  ylab("dN/dS")+
  geom_segment(aes(x=0,y=mean.only.naive,xend=20,yend=mean.only.naive),color="grey30")+
  geom_segment(aes(x=21,y=mean.only.vax,xend=33,yend=mean.only.vax),color="grey30")+
  scale_x_discrete(limits = c("432686","432870","433227", "435786","435805","438577",
                              "442978","444332","444446","444633","445602","449650",
                              "450241","450348","451152","451709","453058","459597",
                              "471467","471588","461913","471876","475670","475951",
                              "481242","481672","482828","484249","486422","487250",
                              "487297","487941"))+
  scale_fill_manual(values = c("white","black"))+
  scale_color_manual(values = c("black","black"))+
  theme_bw()+
  theme(axis.title = element_text(size = 22), 
        axis.text = element_text(size = 19),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15))
ggsave("figs/dnds_comparison.png")




