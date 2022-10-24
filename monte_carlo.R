library(tidyverse)
library(rio)
library(here)

load("comparison_list.RData")
comp.list <- lapply(comp.list,select,-c(POS,REF,ALT))

#group tables by day of sampling
split.time1 <- function(comp.user){
  time1 <- comp.user[,c(1,2,3)]
  return(time1)
}
split.time2 <- function(comp.user){
  time2 <- comp.user[,c(1,4,5)]
  return(time2)
}
split.time3 <- function(comp.user){
  time3 <- comp.user[,c(1,6,7)]
  return(time3)
}
split.time4 <- function(comp.user){
  time4 <- comp.user[,c(1,8,9)]
  return(time4)
}
split.time5 <- function(comp.user){
  time5 <- comp.user[,c(1,10,11)]
  return(time5)
}
split.time6 <- function(comp.user){
  time6 <- comp.user[,c(1,12,13)]
  return(time6)
}
split.times <- list(split.time1,split.time2,split.time3,split.time4,
                    split.time5,split.time6)

final.times <- list()
for(i in seq_along(comp.list)){
  final.times[[i]] <- as.numeric(gsub("SALIVA_FREQ_","",last(colnames(comp.list[[i]]))))
}
final.times <- unlist(final.times)

time.splitter <- function(final.time,user){
  time.splits <- list()
  for(i in 1:final.time){
    time.splits[[i]] <- split.times[[i]](user)
  }
  return(time.splits)
}

time.splits <- list()
for(i in seq_along(comp.list)){
  time.splits[[i]] <- time.splitter(final.times[[i]],comp.list[[i]])
}
names(time.splits) <- names(comp.list)

#randomization
daily.scrambler <- function(user.day){
  snp <- select(user.day,SNP)
  no.snp <- select(user.day,-SNP)
  row.shuffle <- no.snp[sample(nrow(no.snp)),]
  
  number <- sample(length(row.shuffle[,1]),1) #this is the number of rows to randomly swap 
  swap.list <- sample(length(row.shuffle[,1]),number)
  to.swap <- row.shuffle[swap.list,]
  no.swap <- row.shuffle[-swap.list,]
  swapped <- to.swap #initialize
  for(i in seq_along(to.swap[,1])){
    swapped[,1][[i]] <- to.swap[,2][[i]]
    swapped[,2][[i]] <- to.swap[,1][[i]]
  }
  
  final.shuffle <- bind_rows(swapped,no.swap)
  final.shuffle <- bind_cols(snp,final.shuffle)
  return(final.shuffle)
}

scrambler100 <- function(user.day){
  daily.scrambled100 <- list()
  for(i in 1:100){
    daily.scrambled100[[i]] <- daily.scrambler(user.day)
  }
  return(daily.scrambled100)
}

scrambled <- list()
for(i in seq_along(time.splits)){
  scrambled[[i]] <- lapply(time.splits[[i]],scrambler100)
}
names(scrambled) <- names(time.splits)

#put the tables back in comp.table format (ie all sample dates in one table)
trials <- c(1:100)
daygroup <- function(user,trials){
  group <- list()
  for(i in seq_along(user)){
    group[[i]] <- user[[i]][[trials]]
  }
  return(group)
}

user.daygroups <- function(user,trials){
  daily.groups <- list()
  for(i in seq_along(trials)){
    daily.groups[[i]] <- daygroup(user,trials[[i]])
  }
  return(daily.groups)
}

scrambled.daily <- list()
for(i in seq_along(scrambled)){
  scrambled.daily[[i]] <- user.daygroups(scrambled[[i]],trials)
}
names(scrambled.daily) <- names(scrambled)

humpty.dumpty <- function(scrambled.daily){ #put the split tables together again :)
  scrambled.join <- list()
  for(i in seq_along(scrambled.daily)){
    scrambled.join[[i]] <- scrambled.daily[[i]] %>%
      reduce(full_join,by = "SNP")
  }
  return(scrambled.join)
}
scrambled.comp <- lapply(scrambled.daily, humpty.dumpty)

#FST calculations
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

FST100 <- function(user){
  nasal.fst <- list()
  for(i in seq_along(user)){
    nasal.fst[[i]] <- calc.nasal.FST(user[[i]])
  }
  nasal.fst <- lapply(nasal.fst,unlist)
  saliva.fst <- list()
  for(i in seq_along(user)){
    saliva.fst[[i]] <- calc.saliva.FST(user[[i]])
  }
  saliva.fst <- lapply(saliva.fst,unlist)
  between.fst <- list()
  for(i in seq_along(user)){
    between.fst[[i]] <- calc.between.FST(user[[i]])
  }
  between.fst <- lapply(between.fst,unlist)
  
  user.fst <- list()
  for(i in seq_along(user)){
    user.fst[[i]] <- matrix(nrow = length(between.fst[[i]]),ncol = 3, data = NA)
    colnames(user.fst[[i]]) <- c("nasal_FST","saliva_FST","between_FST")
    length(nasal.fst[[i]]) <- length(between.fst[[i]])
    length(saliva.fst[[i]]) <- length(between.fst[[i]])
    user.fst[[i]][,1] <- nasal.fst[[i]]
    user.fst[[i]][,2] <- saliva.fst[[i]]
    user.fst[[i]][,3] <- between.fst[[i]]
    user.fst[[i]] <- as.data.frame(user.fst[[i]])
  } 
  #same day comparison FSTs (ie n1 vs s1, n2 vs s2, etc) are always the same
  #(for all 100 replicates)
  #this makes sense since its the same pool every time
  return(user.fst)
}

fst.simulations <- lapply(scrambled.comp,FST100)
combined.fst.sims <- lapply(fst.simulations, bind_rows)

sim.stat.matrix <- function(sim){
  stat.matrix <- matrix(nrow = length(sim), ncol = 3)
  colnames(stat.matrix) <- c("nasal_vs_saliva","nasal_vs_between","saliva_vs_between")
  
  for(i in seq_along(sim)){
    stat.matrix[i,1] <- (t.test(sim[[i]]$nasal_FST,sim[[i]]$saliva_FST,var.equal = T))$p.value
    stat.matrix[i,2] <- (t.test(sim[[i]]$nasal_FST,sim[[i]]$between_FST,var.equal = T))$p.value
    stat.matrix[i,3] <- (t.test(sim[[i]]$saliva_FST,sim[[i]]$between_FST,var.equal = T))$p.value
  }
  stat.matrix <- as.data.frame(stat.matrix)
  return(stat.matrix)
}
sim.stat.matrices <- lapply(fst.simulations, sim.stat.matrix)

p.count <- list()
for(i in seq_along(sim.stat.matrices)){
  p.count[[i]] <- list("nasal_vs_saliva" = NA,"nasal_vs_between"=NA,"saliva_vs_between"=NA)
  p.count[[i]][1] <- length(which(sim.stat.matrices[[i]]$nasal_vs_saliva < 0.05))
  p.count[[i]][2] <- length(which(sim.stat.matrices[[i]]$nasal_vs_between < 0.05))
  p.count[[i]][3] <- length(which(sim.stat.matrices[[i]]$saliva_vs_between < 0.05))
  p.count[[i]] <- unlist(p.count[[i]])
}
p.count <- bind_rows(p.count)
rownames(p.count) <- names(sim.stat.matrices)

#sim vs real
load("user_fst.RData")
real.comparisons <- function(real,fake){
  compare.nasal <- list()
  compare.saliva <- list()
  compare.between <- list()
  for(i in seq_along(fake)){
    compare.nasal[[i]] <- (t.test(real$nasal_FST, fake[[i]]$nasal_FST,var.equal = T))$p.value
    compare.saliva[[i]] <- (t.test(real$saliva_FST, fake[[i]]$saliva_FST,var.equal = T))$p.value
    compare.between[[i]] <- (t.test(real$between_FST, fake[[i]]$between_FST,var.equal = T))$p.value
  }
  compare.nasal <- unlist(compare.nasal)
  compare.saliva <- unlist(compare.saliva)
  compare.between <- unlist(compare.between)
  compare <- bind_cols(compare.nasal,compare.saliva,compare.between)
  colnames(compare) <- c("p_nasal","p_saliva","p_between")
  return(compare)
}

p.comparisons <- list()
for(i in seq_along(user.fst)){
  p.comparisons[[i]] <- real.comparisons(user.fst[[i]], fst.simulations[[i]])
}
names(p.comparisons) <- names(user.fst)
average.p.comp <- as.data.frame(matrix(ncol = 3,nrow = length(p.comparisons)))
colnames(average.p.comp) <- c("p_nasal","p_saliva","p_between")
rownames(average.p.comp) <- names(p.comparisons)
for(i in seq_along(p.comparisons)){
  average.p.comp[i,1] <- mean(p.comparisons[[i]]$p_nasal)
  average.p.comp[i,2] <- mean(p.comparisons[[i]]$p_saliva)
  average.p.comp[i,3] <- mean(p.comparisons[[i]]$p_between)
}

#histogram distributions
real <- bind_rows(user.fst)
fake <- lapply(fst.simulations,bind_rows)
fake <- bind_rows(fake)

t.test(real$nasal_FST, fake$nasal_FST) #p-value < 2.2e-16
t.test(real$saliva_FST, fake$saliva_FST) #p-value < 2.2e-16
t.test(real$between_FST, fake$between_FST) #p-value < 2.2e-16

nasal.real <- real %>% select(nasal_FST)
nasal.fake <- fake %>% select(nasal_FST)

ggplot(data = nasal.real, aes(x = nasal_FST))+
  geom_histogram(aes(y=stat(density), fill = "Real"), alpha = 0.3, 
                 bins = 70, color = "red")+
  geom_histogram(data = nasal.fake, aes(y = stat(density),x = nasal_FST,
                                         fill = "Simulated"),
                 alpha = 0.3, bins = 70, color = "blue")+
  geom_density(data = nasal.real, aes(x = nasal_FST), color = "red")+
  geom_density(data = nasal.fake, aes(x = nasal_FST), color = "blue")+
  scale_color_manual(values = c("red","blue"))+
  scale_fill_manual(values = c("red","blue"), name = "Data")+
  xlab("FST")+
  ylab("Density")+
  ggtitle("Nasal FST")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        plot.title = element_text(size = 22))+
  guides(fill = guide_legend(override.aes = list(color = "white")))
ggsave("figs/nasal_simulation.png")

saliva.real <- real %>% select(saliva_FST)
saliva.fake <- fake %>% select(saliva_FST)

ggplot(data = saliva.real, aes(x = saliva_FST))+
  geom_histogram(aes(y=stat(density), fill = "Real"), alpha = 0.3, 
                 bins = 70, color = "red")+
  geom_histogram(data = saliva.fake, aes(y = stat(density),x = saliva_FST,
                                          fill = "Simulated"),
                 alpha = 0.3, bins = 70, color = "blue")+
  geom_density(data = saliva.real, aes(x = saliva_FST), color = "red")+
  geom_density(data = saliva.fake, aes(x = saliva_FST), color = "blue")+
  scale_color_manual(values = c("red","blue"))+
  scale_fill_manual(values = c("red","blue"), name = "Data")+
  xlab("FST")+
  ylab("Density")+
  ggtitle("Saliva FST")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        plot.title = element_text(size = 22))+
  guides(fill = guide_legend(override.aes = list(color = "white")))
ggsave("figs/saliva_simulation.png")

between.real <- real %>% select(between_FST)
between.fake <- fake %>% select(between_FST)

ggplot(data = between.real, aes(x = between_FST))+
  geom_histogram(aes(y=stat(density), fill = "Real"), alpha = 0.3, 
                 bins = 70, color = "red")+
  geom_histogram(data = between.fake, aes(y = stat(density),x = between_FST,
                                          fill = "Simulated"),
                 alpha = 0.3, bins = 70, color = "blue")+
  geom_density(data = between.real, aes(x = between_FST), color = "red")+
  geom_density(data = between.fake, aes(x = between_FST), color = "blue")+
  scale_color_manual(values = c("red","blue"))+
  scale_fill_manual(values = c("red","blue"), name = "Data")+
  xlab("FST")+
  ylab("Density")+
  ggtitle("Nasal-Saliva FST")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 19),
        plot.title = element_text(size = 22))+
  guides(fill = guide_legend(override.aes = list(color = "white")))
ggsave("figs/between_simulation.png")



