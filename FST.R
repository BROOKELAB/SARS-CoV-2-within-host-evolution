library(tidyverse)
library(rio)
library(gdata)
library(data.table)
library(here)

load("comparison_list.RData")

calc.FST <- function(sample1,sample2){
  sample1 <- lapply(sample1, as.numeric)
  sample2 <- lapply(sample2, as.numeric)
  
  avg.freq <- list()
  s1.haplotype <- list()
  s2.haplotype <- list()
  pop.haplotype <- list()
  
  for(i in seq_along(sample1)) {
    if((!is.na(sample1[[i]])) && (!is.na(sample2[[i]]))){
      avg.freq[[i]] <- mean(c(sample1[[i]],sample2[[i]]))
      s1.haplotype[[i]] <- 1 - (sample1[[i]]^2 + 
                                  ((1 - sample1[[i]])^2))
      s2.haplotype[[i]] <- 1 - (sample2[[i]]^2 + 
                                  ((1 - sample2[[i]])^2))
      pop.haplotype[[i]] <- 1 - (avg.freq[[i]]^2 + 
                                   ((1 - avg.freq[[i]])^2))
    }
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
  return(FST)
}

calc.nasal.FST <- function(user){
  nasal <- user[,grep("NASAL",colnames(user))]
  cj <- expand.grid(seq_along(nasal), seq_along(nasal))
  cj <- cj[-which(cj$Var1 >= cj$Var2),]
  fst <- list()
  for(i in seq_along(cj$Var1)){
    fst[[i]] <- calc.FST(nasal[[cj$Var1[[i]]]],nasal[[cj$Var2[[i]]]])
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

stat.matrix <- matrix(nrow = length(comp.list), ncol = 3)
colnames(stat.matrix) <- c("nasal_vs_saliva","nasal_vs_between","saliva_vs_between")
rownames(stat.matrix) <- c(names(comp.list))

for(i in seq_along(rownames(stat.matrix))){
  if(length(which(!is.na(nasal.fst[[i]]))) > 1 &&
     length(which(!is.na(saliva.fst[[i]]))) > 1){
    stat.matrix[i,1] <- (t.test(nasal.fst[[i]],saliva.fst[[i]],var.equal = T))$p.value
  } else {
    stat.matrix[i,1] <- NA
  }
  stat.matrix[i,2] <- (t.test(nasal.fst[[i]],between.fst[[i]],var.equal = T))$p.value
  stat.matrix[i,3] <- (t.test(saliva.fst[[i]],between.fst[[i]],var.equal = T))$p.value
}

total.matrix <- matrix(nrow = 1, ncol = 3)
colnames(total.matrix) <- c("nasal_vs_saliva","nasal_vs_between","saliva_vs_between")
rownames(total.matrix) <- "total"
total.matrix[1,1] <- (t.test(all.nasal.fst,all.saliva.fst,var.equal = T))$p.value
total.matrix[1,2] <- (t.test(all.nasal.fst,all.between.fst,var.equal = T))$p.value
total.matrix[1,3] <- (t.test(all.saliva.fst,all.between.fst,var.equal = T))$p.value

stat.matrix <- rbind(stat.matrix,total.matrix)

write.csv(stat.matrix,"FST_significance.csv")


#heatmaps
heatmap.matrix <- function(user,fstuser){
  nasal <- user[,grep("NASAL",colnames(user))]
  saliva <- user[,grep("SALIVA",colnames(user))]
  fullgrid <- expand.grid(seq_along(nasal), seq_along(saliva))
  grid <- fullgrid[-which(fullgrid$Var1 >= fullgrid$Var2),]
  usermat <- as.data.frame(matrix(ncol = 3, 
                                  nrow = length(which(!is.na(fstuser$nasal_FST))) + 
                                    length(which(!is.na(fstuser$saliva_FST)))+
                                    length(fstuser$between_FST) + 
                                    (2*length(nasal))))
  colnames(usermat) <- c("Env1","Env2","FST")
  
  usermat$Env1 <- c(paste0("N",grid$Var1), paste0("S",grid$Var1), 
                    paste0("N",fullgrid$Var1),paste0("N",seq_along(nasal)),
                    paste0("S",seq_along(saliva)))
  usermat$Env2 <- c(paste0("N",grid$Var2), paste0("S",grid$Var2), 
                    paste0("S",fullgrid$Var2),paste0("N",seq_along(nasal)),
                    paste0("S",seq_along(saliva)))
  usermat$FST <- c(fstuser$nasal_FST[which(!is.na(fstuser$nasal_FST))],
                   fstuser$saliva_FST[which(!is.na(fstuser$saliva_FST))],
                   fstuser$between_FST,
                   rep(0,2*length(nasal)))
  return(usermat)
}

FST.matrix <- list()
for(i in seq_along(comp.list)){
  FST.matrix[[i]] <- heatmap.matrix(comp.list[[i]],user.fst[[i]])
}
names(FST.matrix) <- gsub("user_","",names(comp.list))

FST.maps <- list()
for(i in seq_along(FST.matrix)){
  FST.maps[[i]] <- ggplot(data = FST.matrix[[i]], aes(x=Env1, y=Env2,fill=FST))+
    geom_tile() +
    scale_fill_viridis_c(limits = c(0,1),
                         breaks = c(0,.25,.5,.75,1))+
    ggtitle(names(FST.matrix)[[i]])+
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
names(FST.maps) <- names(FST.matrix)

ggpubr::ggarrange(FST.maps$`459597`,FST.maps$`433227`,FST.maps$`432686`,
                  nrow = 1, ncol = 3)



