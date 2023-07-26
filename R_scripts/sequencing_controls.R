library(tidyverse)
library(rio)
library(gtools)
library(here)

#2264562, 2390233, 2415377, 2415388, 2421070, 2425945, 2436517, 2447784, 2447896

bar_2264562 <- import("naive_output_tables/all_ivar/2264562.ivar.tsv")
bar_2390233 <- import("naive_output_tables/all_ivar/2390233.ivar.tsv")
bar_2415377 <- import("naive_output_tables/all_ivar/2415377.ivar.tsv")
bar_2415388 <- import("naive_output_tables/all_ivar/2415388.ivar.tsv")
bar_2421070 <- import("naive_output_tables/all_ivar/2421070.ivar.tsv")
bar_2425945 <- import("naive_output_tables/all_ivar/2425945.ivar.tsv")
bar_2436517 <- import("naive_output_tables/all_ivar/2436517.ivar.tsv")
bar_2447784 <- import("naive_output_tables/all_ivar/2447784.ivar.tsv")
bar_2447896 <- import("naive_output_tables/all_ivar/2447896.ivar.tsv")


run1 <- list(bar_2264562, bar_2390233, bar_2415377, bar_2415388, bar_2421070,
          bar_2425945, bar_2436517, bar_2447784, bar_2447896)

names(run1) <- c("2264562", "2390233", "2415377", "2415388", "2421070", 
                   "2425945", "2436517", "2447784", "2447896")

run2 <- lapply((paste0("sequence_controls/",list.files("sequence_controls/", 
                                                           pattern = "ivar"))),
                  import)
names(run2) <- names(run1)

#test correlation between iSNVs >0.03 freq and >1000 reads
run1.cut <- lapply(run1, mutate, "SNP" = paste0(REF,POS,ALT))
run1.cut <- lapply(run1.cut, filter, ALT_FREQ >= 0.03)
run1.cut <- lapply(run1.cut, filter, TOTAL_DP >= 1000)
run1.cut <- lapply(run1.cut, select, SNP, ALT_FREQ)
run1.cut <- lapply(run1.cut, distinct)

run2.cut <- lapply(run2, mutate, "SNP" = paste0(REF,POS,ALT))
run2.cut <- lapply(run2.cut, filter, ALT_FREQ >= 0.03)
run2.cut <- lapply(run2.cut, filter, TOTAL_DP >= 1000)
run2.cut <- lapply(run2.cut, select, SNP, ALT_FREQ)
run2.cut <- lapply(run2.cut, distinct)

run.compare <- map2(run1.cut,run2.cut,full_join, by = "SNP")
run.compare <- lapply(run.compare, distinct)
run.cors <- list()
for(i in seq_along(run.compare)){
  run.compare[[i]][is.na(run.compare[[i]])] <- 0
  colnames(run.compare[[i]])[c(2:3)] <- c("R1","R2")
  run.cors[[i]] <- cor.test(run.compare[[i]]$R1, run.compare[[i]]$R2, method = "p")$estimate
}
run.cors <- unlist(run.cors)
save(run.cors, file = "run_cors.RData") #will use for quality control in future analyses

#test just shared from 444332
#2415377 #2425945 #2436517 #2447784
user.run1 <- list(run1$`2415377`, run1$`2421070`, run1$`2425945`,run1$`2436517`,run1$`2447784`)
user.run1 <- lapply(user.run1, select, POS,REF,ALT,ALT_FREQ,TOTAL_DP)
names(user.run1) <- c("2415377","2421070","2425945","2436517","2447784")

user.run2 <- list(run2$`2415377`, run2$`2421070`, run2$`2425945`,run2$`2436517`,run2$`2447784`)
user.run2 <- lapply(user.run2, select, POS,REF,ALT,ALT_FREQ,TOTAL_DP)
names(user.run2) <- names(user.run1)

#filter
run.filter <- function(run){
  run <- run %>%
    filter(ALT_FREQ >= 0.03)%>%
    filter(TOTAL_DP >= 1000)
  run <- run[-which(run$POS %in% c(6696,11074,15965,29051,187,1059,2094,3037,
                                   3130,6696,6990,8022,10323,10741,11074,13408,
                                   14786,19684,20148,21137,24034,24378,25563,26144,
                                   26461,26681,28077,28826,28854,29051,29700,29760)),]
  run <- run %>%
    select(POS,REF,ALT,ALT_FREQ)
  return(run)
}
user.run1.filter <- lapply(user.run1,run.filter)
user.run2.filter <- lapply(user.run2, run.filter)

#intersect
allpos.run1 <- lapply(user.run1.filter,select, POS,REF,ALT)
allpos.run2 <- lapply(user.run2.filter,select, POS,REF,ALT)

snp.intersect <- function(allpos){
  expand <- expand.grid(seq_along(allpos), seq_along(allpos))
  expand <- expand[,c("Var2","Var1")]
  expand <- expand[which(expand$Var1 != expand$Var2),]
  intersecting <- vector(mode = "list", length = length(expand))
  for(i in seq_along(expand$Var2)){
    intersecting[[i]] <- dplyr::intersect(allpos[[expand$Var2[[i]]]],
                                          allpos[[expand$Var1[[i]]]])
  }
  flat <- bind_rows(intersecting)
  common <- distinct(flat)
  return(common)
}

intersect.run1 <- snp.intersect(allpos.run1)
intersect.run2 <- snp.intersect(allpos.run2)

varianttable <- function(intersectuser,everythinguser){
  freq <- list()
  for(i in 2:(length(everythinguser)+1)){
    freq[[1]] <- intersectuser
    freq[[i]] <- everythinguser[[i-1]]
  }
  joiner <- freq %>%
    purrr::reduce(left_join, by=c("POS","REF","ALT"))%>%
    distinct()
  numvec <- 1:((dim(joiner)[[2]]-3) - ((dim(joiner)[[2]]-3)/2))
  
  numrep <- list()
  for(i in 1:length(numvec)){
    numrep[[i]]<- rep(numvec[[i]],2)
  }
  numrep <- unlist(numrep)
  colrep <- rep(c("freq_","depth_"),length(numvec))
  colnames(joiner) <- c("POS","REF","ALT",paste0(colrep,numrep))
  
  for(i in 1:dim(joiner)[[1]]){
    for(j in 1:dim(joiner)[[2]]){
      if(dim(joiner)[[1]] != 0){
        if(strsplit(colnames(joiner)[[j]],"_")[[1]][1] == "depth"){
          if((!is.na(joiner[i,j])) && (joiner[i,j]<1000)){
            joiner[i,(j-1)] <- paste(joiner[i,(j-1)],"(low dp)")
          }
        }
      }
    }
  }
  
  colvec <- strsplit(colnames(joiner),"_")
  keepvec <- list()
  for(i in 1:length(colvec)){
    if(colvec[[i]][1]!= "depth"){
      keepvec[[i]] <-i
    }
  }
  keepvec <- keepvec[!sapply(keepvec,is.null)]
  keepvec <- unlist(keepvec)
  joiner <- joiner[,keepvec]
  return(joiner)
} 

run1.table <- varianttable(intersect.run1, user.run1)
run2.table <- varianttable(intersect.run2, user.run2)

#replace NA with 0.01 
NA_01 <- function(vartable){
  if(dim(vartable)[[1]] != 0){
    for(i in 1:length(vartable$POS)){
      for(j in 1:length(colnames(vartable))){
        if(is.na(vartable[i,j])){
          vartable[i,j] <- 0.01
        }
      }
    }
  }
  return(vartable)
}
run1.table <- NA_01(run1.table)
run2.table <- NA_01(run2.table)

#remove consensus variants
remove.consensus <- function(user){
  info <- user[,1:3]
  freqs <- user[,grep("freq",colnames(user))]
  keeps <- freqs
  for(i in seq_along(rownames(keeps))){
    for(j in seq_along(colnames(keeps))){
      if(!identical(grep("(low dp)", keeps[i,j]),integer(0))){
        keeps[i,j] <- 1
      }
      if(length(which(keeps[i,] < 0.97)) < 2){
        keeps[i,] <- "X"
      }
    }
  }
  info <- info[which(keeps$freq_1 != "X"),]
  freqs <- freqs[which(keeps$freq_1 != "X"),]
  final <- cbind(info,freqs)
  return(final)
}

run1.table <- remove.consensus(run1.table)
run2.table <- remove.consensus(run2.table)

#generate list of shared variants
intersect.run1 <- run1.table[,c(1:3)] %>%
  mutate("SNP" = paste0(REF,POS,ALT)) %>%
  relocate(SNP, .before = POS)
intersect.run2 <- run2.table[,c(1:3)] %>%
  mutate("SNP" = paste0(REF,POS,ALT)) %>%
  relocate(SNP, .before = POS)

user.run1 <- lapply(user.run1, mutate, "SNP" = paste0(REF,POS,ALT))
user.run2 <- lapply(user.run2, mutate, "SNP" = paste0(REF,POS,ALT))

for(i in seq_along(user.run1)){
  user.run1[[i]] <- user.run1[[i]][which(user.run1[[i]]$SNP %in% intersect.run1$SNP),]
  user.run2[[i]] <- user.run2[[i]][which(user.run2[[i]]$SNP %in% intersect.run2$SNP),]
}

user.run1 <- lapply(user.run1, select, SNP,ALT_FREQ)
user.run2 <- lapply(user.run2, select, SNP,ALT_FREQ)

date.label <- function(user.run){
  for(i in seq_along(user.run)){
    if(names(user.run)[[i]] == "2415377"){
      user.run[[i]] <- user.run[[i]] %>%
        mutate("DAY" = 1)
    }
    if(names(user.run)[[i]] == "2421070"){
      user.run[[i]] <- user.run[[i]] %>%
        mutate("DAY" = 2)
    }
    if(names(user.run)[[i]] == "2425945"){
      user.run[[i]] <- user.run[[i]] %>%
        mutate("DAY" = 3)
    }
    if(names(user.run)[[i]] == "2436517"){
      user.run[[i]] <- user.run[[i]] %>%
        mutate("DAY" = 5)
    }
    if(names(user.run)[[i]] == "2447784"){
      user.run[[i]] <- user.run[[i]] %>%
        mutate("DAY" = 6)
    }
  }
  return(user.run)
}
user.run1 <- date.label(user.run1)
user.run2 <- date.label(user.run2)

user.table <- list()
for(i in seq_along(user.run1)){
  user.table[[i]] <- full_join(user.run1[[i]],user.run2[[i]], by = c("SNP","DAY"))
  user.table[[i]] <- user.table[[i]] %>%
    relocate("DAY", .after = "SNP")
  colnames(user.table[[i]]) <- c("SNP","DAY","RUN1","RUN2")
  user.table[[i]] <- distinct(user.table[[i]])
  user.table[[i]][is.na(user.table[[i]])] <- 0
}

user.table <- bind_rows(user.table)
cor.test(user.table$RUN1, user.table$RUN2, method = "p") #cor = 0.9754242 #p-value = 0.02458

export(user.table, "runcontrol_user.csv")

#test for Ct effect
#2264562, 2390233, 2415388, 2421070, 2447896
#432870  #435786  #435805  #444332  #449650
ct.list <- as_tibble(c(21.34, 25.8, 27.94,23.93,24.9))

ct.run1 <- list(run1$`2264562`,run1$`2390233`,run1$`2415388`,run1$`2421070`,run1$`2447896`)
ct.run1 <- lapply(ct.run1, mutate,"SNP" = paste0(REF,POS,ALT))
names(ct.run1) <- c("2264562","2390233","2415388","2421070","2447896")
ct.run2 <- list(run2$`2264562`,run2$`2390233`,run2$`2415388`,run2$`2421070`,run2$`2447896`)
ct.run2 <- lapply(ct.run2, mutate,"SNP" = paste0(REF,POS,ALT))
names(ct.run2) <- names(ct.run1)


load("naive_intersecting.RData")
naive.intersecting <- lapply(naive.intersecting, mutate, "SNP" = paste0(REF,POS,ALT))

ct.run1$`2264562` <- ct.run1$`2264562`[which(ct.run1$`2264562`$SNP %in% 
                                               naive.intersecting$user_432870$SNP),]
ct.run1$`2390233` <- ct.run1$`2390233` [which(ct.run1$`2390233`$SNP %in% 
                                               naive.intersecting$user_435786$SNP),]
ct.run1$`2415388` <- ct.run1$`2415388` [which(ct.run1$`2415388`$SNP %in% 
                                                naive.intersecting$user_435805$SNP),]
ct.run1$`2421070` <- ct.run1$`2421070` [which(ct.run1$`2421070`$SNP %in% 
                                                naive.intersecting$user_444332$SNP),]
ct.run1$`2447896` <- ct.run1$`2447896` [which(ct.run1$`2447896`$SNP %in% 
                                                naive.intersecting$user_449650$SNP),]

ct.run2$`2264562` <- ct.run2$`2264562`[which(ct.run2$`2264562`$SNP %in% 
                                               naive.intersecting$user_432870$SNP),]
ct.run2$`2390233` <- ct.run2$`2390233` [which(ct.run2$`2390233`$SNP %in% 
                                                naive.intersecting$user_435786$SNP),]
ct.run2$`2415388` <- ct.run2$`2415388` [which(ct.run2$`2415388`$SNP %in% 
                                                naive.intersecting$user_435805$SNP),]
ct.run2$`2421070` <- ct.run2$`2421070` [which(ct.run2$`2421070`$SNP %in% 
                                                naive.intersecting$user_444332$SNP),]
ct.run2$`2447896` <- ct.run2$`2447896` [which(ct.run2$`2447896`$SNP %in% 
                                                naive.intersecting$user_449650$SNP),]

ct.run1 <- lapply(ct.run1,select,SNP,ALT_FREQ)
ct.run2 <- lapply(ct.run2, select,SNP,ALT_FREQ)

for(i in seq_along(ct.run1)){
  ct.run1[[i]] <- mutate(ct.run1[[i]], "Ct" = rep(ct.list[i,], length(ct.run1[[i]]$SNP)))
  ct.run2[[i]] <- mutate(ct.run2[[i]], "Ct" = rep(ct.list[i,], length(ct.run2[[i]]$SNP)))
}

ct.table <- list()
for(i in seq_along(ct.run1)){
  ct.table[[i]] <- full_join(ct.run1[[i]],ct.run2[[i]], by = c("SNP","Ct"))
  colnames(ct.table[[i]]) <- c("SNP","RUN1","Ct","RUN2")
  ct.table[[i]] <- relocate(ct.table[[i]],Ct, .before = SNP)
  ct.table[[i]] <- distinct(ct.table[[i]])
  ct.table[[i]][is.na(ct.table[[i]])] <- 0
}
ct.table[[1]] <- ct.table[[1]] %>%
  mutate("PARTICIPANT" = 432870)%>%
  mutate("DAY" = 1)%>%
  relocate(DAY, .before = Ct)%>%
  relocate(PARTICIPANT, .before = DAY)
ct.table[[2]] <- ct.table[[2]] %>%
  mutate("PARTICIPANT" = 435786)%>%
  mutate("DAY" = 3)%>%
  relocate(DAY, .before = Ct)%>%
  relocate(PARTICIPANT, .before = DAY)
ct.table[[3]] <- ct.table[[3]] %>%
  mutate("PARTICIPANT" = 435805)%>%
  mutate("DAY" = 9) %>%
  relocate(DAY, .before = Ct)%>%
  relocate(PARTICIPANT, .before = DAY)
ct.table[[4]] <- ct.table[[4]] %>%
  mutate("PARTICIPANT" = 444332) %>%
  mutate("DAY" = 2)%>%
  relocate(DAY, .before = Ct)%>%
  relocate(PARTICIPANT, .before = DAY)
ct.table[[5]] <- ct.table[[5]] %>%
  mutate("PARTICIPANT" = 449650) %>%
  mutate("DAY" = 1)%>%
  relocate(DAY, .before = Ct)%>%
  relocate(PARTICIPANT, .before = DAY)

ct.all <- as.data.frame(bind_rows(ct.table))
ct.all$Ct <- as.double(ct.all$Ct)
ct.all <- arrange(ct.all, Ct, desc = T)
cor.test(ct.all$RUN1, ct.all$RUN2, method = "p") 
#cor = 0.9966904, p-value < 2.2e-16
export(ct.all, file = "runcontrol_ct.csv")

#filter out high freq and low freq
ct.all.filter <- ct.all %>%
  filter(RUN1 > 0.03 | RUN2 > 0.03) %>%
  filter(RUN1 < 0.97 | RUN2 < 0.97)
cor.test(ct.all.filter$RUN1, ct.all.filter$RUN2, method = "p")
#cor = 0.8287419 #p = 0.04148

#combine
ct.all.filter <- ct.all.filter[,-c(1:3)]
user.table <- user.table[,-2]
runcontrol <- bind_rows(ct.all.filter, user.table)
cor.test(runcontrol$RUN1, runcontrol$RUN2, method = "p")
#cor = 0.833268 #p = 0.00275
save(runcontrol, file = "runcontrol.RData")






