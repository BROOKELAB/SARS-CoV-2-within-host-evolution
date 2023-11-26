library(tidyverse)
library(rio)
library(compiler)
library(here)

load("user_fst.RData")

set.seed(7)

#randomize Fst values across groups for each user
randomize.tables <- function(user){
  all <- unname(unlist(user))
  all <- all[which(!is.na(all))]
  all <- sample(all, length(all))
  nasal <- all[1:length(na.omit(user$nasal_FST))] 
  saliva <- all[(length(nasal)+1) : ((length(nasal)) + (length(na.omit(user$saliva_FST))))]
  between <- all[(length(nasal)+ (length(saliva)+1)) : ((length(nasal)+ (length(saliva))) + (length(na.omit(user$between_FST))))]
  
  length(nasal) <- length(between)
  length(saliva) <- length(between)
  
  nasal <- as.data.frame(nasal)
  saliva <- as.data.frame(saliva)
  between <- as.data.frame(between)
  
  new <- as.data.frame(bind_cols(nasal, saliva, between))
  colnames(new) <- c("nasal_FST","saliva_FST","between_FST")
  return(new)
}

#repeat randomization x10000 per user
# + pull out differences in group means for each randomized dataset
multi.random <- function(user){
  scrambled.tabs <- list()
  scrambled.ns <- list()
  scrambled.nb <- list()
  scrambled.sb <- list()
  for(i in 1:10000){
    scrambled.tabs[[i]] <- randomize.tables(user)
  }
  for(i in seq_along(scrambled.tabs)){
    scrambled.ns[[i]] <- abs(mean(na.omit((scrambled.tabs[[i]]$nasal_FST))) - 
                               mean(na.omit(scrambled.tabs[[i]]$saliva_FST)))
    scrambled.nb[[i]] <- abs(mean(na.omit((scrambled.tabs[[i]]$nasal_FST))) - 
                          mean(na.omit(scrambled.tabs[[i]]$between_FST)))
    scrambled.sb[[i]] <- abs(mean(na.omit((scrambled.tabs[[i]]$saliva_FST))) - 
                          mean(na.omit(scrambled.tabs[[i]]$between_FST)))
  }
  scrambled.stats <- list(scrambled.ns, scrambled.nb, scrambled.sb)
  names(scrambled.stats) <- c("nasal_saliva", "nasal_between", "saliva_between")
  return(scrambled.stats)
}
multi.random <- cmpfun(multi.random)

random.stats <- lapply(user.fst, multi.random) #this will take a bit

#randomize x10000 for combined cross-user dataset + pull out stats
all.users <- bind_rows(user.fst)
total.random.stats <- multi.random(all.users)

#compare randomized tables to real data (per user)
p.calc <- function(user, random){
  trials <- 10000
  
  user.ns <- abs(mean(na.omit(user$nasal_FST)) - mean(na.omit(user$saliva_FST)))
  user.nb <- abs(mean(na.omit(user$nasal_FST)) - mean(na.omit(user$between_FST)))
  user.sb <- abs(mean(na.omit(user$saliva_FST)) - mean(na.omit(user$between_FST)))
  
  if(!is.nan(user.ns)){
    p.ns <- length(which(random$nasal_saliva >= user.ns)) / trials
    p.nb <- length(which(random$nasal_between >= user.nb)) / trials
    p.sb <- length(which(random$saliva_between >= user.sb)) / trials
  } else {
    p.ns <- NA
    p.nb <- NA
    p.sb <- NA
  }
  p.list <- c(p.ns, p.nb, p.sb)
  p.list[which(p.list == 0)] <- "<0.0001"
  names(p.list) <- c("nasal_vs_saliva", "nasal_vs_between", "saliva_vs_between")
  return(p.list)
}
p.mat <- as.data.frame(bind_rows(map2(user.fst, random.stats, p.calc)))

#compare randomized table to real data (overall)
p.overall <- p.calc(all.users, total.random.stats)
fst.sig.matrix <- bind_rows(p.mat, p.overall)
rownames(fst.sig.matrix) <- c(names(user.fst),"total")

write.csv(fst.sig.matrix, "FST_significance_mc.csv")




