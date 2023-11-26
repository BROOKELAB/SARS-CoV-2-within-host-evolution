library(tidyverse)
library(rio)
library(gdata)
library(Nmisc)
library(dunn.test)
library(here)

#function to calculate FST
calc.FST <- function(sample1,sample2){
  sample1 <- unlist(lapply(sample1, as.numeric))
  sample2 <- unlist(lapply(sample2, as.numeric))
  sample.tab <- as_tibble(data.frame(matrix(data = NA, ncol = 2)))
  sample.tab <- na.omit(as_tibble(cbind(sample1,sample2)))
  
  avg.freq <- list()
  s1.haplotype <- list()
  s2.haplotype <- list()
  pop.haplotype <- list()
  
  if((dim(sample.tab)[1]) > 0){ 
    sample1 <- sample.tab$sample1
    sample2 <- sample.tab$sample2
    for(i in seq_along(sample1)) {
      avg.freq[[i]] <- mean(c(sample1[[i]],sample2[[i]]))
      s1.haplotype[[i]] <- 1 - (sample1[[i]]^2 + 
                                  ((1 - sample1[[i]])^2))
      s2.haplotype[[i]] <- 1 - (sample2[[i]]^2 + 
                                  ((1 - sample2[[i]])^2))
      pop.haplotype[[i]] <- 1 - (avg.freq[[i]]^2 + 
                                   ((1 - avg.freq[[i]])^2))
    } 
    s1.haplotype <- unlist(s1.haplotype)
    s2.haplotype <- unlist(s2.haplotype)
    pop.haplotype <- unlist(pop.haplotype)
    s1.mean <- mean(s1.haplotype)
    s2.mean <- mean(s2.haplotype)
    HS <- mean(c(s1.mean,s2.mean))
    HT <- mean(pop.haplotype)
    FST <- (HT-HS)/HT
    if(FST == "NaN"){ 
      FST <- 0
    }
  } else {
    FST <- NA
  }
  return(FST)
}

#calculate threshold FST based on between-run variation
load("runcontrol.RData")
threshold <- calc.FST(runcontrol$RUN1,runcontrol$RUN2) #this will be the FST threshold of detection
#threshold: 0.001276611

#calculate within-nasal, within-saliva, between-environment FSTs for each user
load("comparison_list.RData")

calc.nasal.FST <- function(user){
  nasal <- user[,grep("NASAL",colnames(user))]
  cj <- expand.grid(seq_along(nasal), seq_along(nasal))
  cj <- cj[-which(cj$Var1 >= cj$Var2),]
  fst <- list()
  for(i in seq_along(cj$Var1)){
    fst[[i]] <- calc.FST(nasal[[cj$Var1[[i]]]],nasal[[cj$Var2[[i]]]])
    if(!is.na(fst[[i]])){
      if(fst[[i]] <= threshold){
        fst[[i]] <- 0
      }
    }
  }
  return(fst)
}
nasal.fst.list <- lapply(comp.list,calc.nasal.FST)
nasal.fst <- lapply(nasal.fst.list, unlist)
all.nasal.fst <- unlist(nasal.fst)

calc.saliva.FST <- function(user){
  saliva <- user[,grep("SALIVA",colnames(user))]
  cj <- expand.grid(seq_along(saliva), seq_along(saliva))
  cj <- cj[-which(cj$Var1 >= cj$Var2),]
  fst <- list()
  for(i in seq_along(cj$Var1)){
    fst[[i]] <- calc.FST(saliva[[cj$Var1[[i]]]],saliva[[cj$Var2[[i]]]])
    if(!is.na(fst[[i]])){
      if(fst[[i]] <= threshold){
        fst[[i]] <- 0
      }
    }
  }
  return(fst)
}
saliva.fst.list <- lapply(comp.list, calc.saliva.FST)
saliva.fst <- lapply(saliva.fst.list, unlist)
all.saliva.fst <- unlist(saliva.fst)

calc.between.FST <- function(user){
  nasal <- user[,grep("NASAL",colnames(user))]
  saliva <- user[,grep("SALIVA",colnames(user))]
  cj <- expand.grid(seq_along(nasal), seq_along(saliva))
  fst <- list()
  for(i in seq_along(cj$Var1)){
    fst[[i]] <- calc.FST(nasal[[cj$Var1[[i]]]],saliva[[cj$Var2[[i]]]])
    if(!is.na(fst[[i]])){
      if(fst[[i]] <= threshold){
        fst[[i]] <- 0
      }
    }
  }
  return(fst)
}
between.fst.list <- lapply(comp.list, calc.between.FST)
between.fst <- lapply(between.fst.list, unlist)
all.between.fst <- unlist(between.fst)

user.fst <- list()
for(i in seq_along(comp.list)){
  user.fst[[i]] <- matrix(nrow = length(between.fst[[i]]),ncol = 3, data = NA)
  colnames(user.fst[[i]]) <- c("nasal_FST","saliva_FST","between_FST")
  length(nasal.fst[[i]]) <- length(between.fst[[i]])
  length(saliva.fst[[i]]) <- length(between.fst[[i]])
  user.fst[[i]][,1] <- nasal.fst[[i]]
  user.fst[[i]][,2] <- saliva.fst[[i]]
  user.fst[[i]][,3] <- between.fst[[i]]
  user.fst[[i]] <- as.data.frame(user.fst[[i]])
}
names(user.fst) <- names(comp.list)
save(user.fst, file = "user_fst.RData")

for(i in seq_along(user.fst)){
  filenames <- paste0("FST_tables/",names(user.fst),".csv")
  write.csv(user.fst[[i]],filenames[[i]])
}

#heatmaps
heatmap.matrix <- function(user,fstuser){
  nasal <- user[,grep("NASAL",colnames(user))]
  saliva <- user[,grep("SALIVA",colnames(user))]
  fullgrid <- expand.grid(seq_along(nasal), seq_along(saliva))
  grid <- fullgrid[-which(fullgrid$Var1 >= fullgrid$Var2),]
  usermat <- as.data.frame(matrix(ncol = 3, 
                                  nrow = length(grid$Var1) + 
                                    length(grid$Var1)+
                                    length(fullgrid$Var1) + 
                                    (2*length(nasal))))
  colnames(usermat) <- c("Env1","Env2","FST")
  
  usermat$Env1 <- c(paste0("N",grid$Var1), paste0("S",grid$Var1), 
                    paste0("N",fullgrid$Var1),paste0("N",seq_along(nasal)),
                    paste0("S",seq_along(saliva)))
  usermat$Env2 <- c(paste0("N",grid$Var2), paste0("S",grid$Var2), 
                    paste0("S",fullgrid$Var2),paste0("N",seq_along(nasal)),
                    paste0("S",seq_along(saliva)))
  usermat$FST <- c(fstuser$nasal_FST[1:length(grid$Var1)],
                   fstuser$saliva_FST[1:length(grid$Var1)],
                   fstuser$between_FST,
                   rep(0,2*length(nasal)))
  return(usermat)
}

comp.list.cut <- comp.list %>%
  discard(~ dim(.x)[1] == 0)

user.fst.cut <- user.fst %>%
  discard(~ length(na.omit(unlist(.x))) == 0) #5

FST.matrix <- map2(comp.list.cut, user.fst.cut, heatmap.matrix)

FST.maps <- list()
for(i in seq_along(FST.matrix)){
  FST.maps[[i]] <- ggplot(data = FST.matrix[[i]], aes(x=Env1, y=Env2,fill=FST))+
    geom_tile() +
    scale_fill_viridis_c(limits = c(0,1),
                         breaks = c(0,.25,.5,.75,1))+
    ggtitle(gsub("user_","",names(FST.matrix)[[i]]))+
    theme_bw()+
    theme(axis.title = element_blank(),
          axis.text = element_text(size = 19,
                                   face = "bold"),
          legend.title = element_text(size = 22,
                                      face = "bold"),
          legend.text = element_text(size = 19,
                                     face = "bold"),
          plot.title = element_text(size = 22,
                                    face = "bold"))
  
}

names(FST.maps) <- gsub("user_","",names(FST.matrix))
plot(FST.maps$`433227`)
ggsave("figs/FST_433227.png")
plot(FST.maps$`471467`)
ggsave("figs/FST_471467.png")
plot(FST.maps$`444633`)
ggsave("figs/FST_444633.png")


