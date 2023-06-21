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

#test just shared from 444332
#2415377 #2425945 #2436517 #2447784
user.run1 <- list(run1$`2415377`, run1$`2421070`, run1$`2425945`,run1$`2436517`,run1$`2447784`)
user.run1 <- lapply(user.run1, mutate,"SNP" = paste0(REF,POS,ALT))
names(user.run1) <- c("2415377","2421070","2425945","2436517","2447784")

user.run2 <- list(run2$`2415377`, run2$`2421070`, run2$`2425945`,run2$`2436517`,run2$`2447784`)
user.run2 <- lapply(user.run2, mutate,"SNP" = paste0(REF,POS,ALT))
names(user.run2) <- names(user.run1)

#filter
run.filter <- function(run){
  run <- run %>%
    filter(ALT_FREQ >= 0.03)%>%
    filter(ALT_FREQ <= 0.97)%>%
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
    filter(POS != 29700) %>%
    filter(POS != 29760) %>%
    #mutate("SNP" = paste0(REF,POS,ALT))%>%
    select(SNP,ALT_FREQ)
  return(run)
}
user.run1.filter <- lapply(user.run1,run.filter)
user.run2.filter <- lapply(user.run2, run.filter)

#intersect
allpos.run1 <- lapply(user.run1.filter,select, SNP)
allpos.run2 <- lapply(user.run2.filter,select, SNP)

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
                                               naive.intersecting$`432870`$SNP),]
ct.run1$`2390233` <- ct.run1$`2390233` [which(ct.run1$`2390233`$SNP %in% 
                                               naive.intersecting$`435786`$SNP),]
ct.run1$`2415388` <- ct.run1$`2415388` [which(ct.run1$`2415388`$SNP %in% 
                                                naive.intersecting$`435805`$SNP),]
ct.run1$`2421070` <- ct.run1$`2421070` [which(ct.run1$`2421070`$SNP %in% 
                                                naive.intersecting$`444332`$SNP),]
ct.run1$`2447896` <- ct.run1$`2447896` [which(ct.run1$`2447896`$SNP %in% 
                                                naive.intersecting$`449650`$SNP),]

ct.run2$`2264562` <- ct.run2$`2264562`[which(ct.run2$`2264562`$SNP %in% 
                                               naive.intersecting$`432870`$SNP),]
ct.run2$`2390233` <- ct.run2$`2390233` [which(ct.run2$`2390233`$SNP %in% 
                                                naive.intersecting$`435786`$SNP),]
ct.run2$`2415388` <- ct.run2$`2415388` [which(ct.run2$`2415388`$SNP %in% 
                                                naive.intersecting$`435805`$SNP),]
ct.run2$`2421070` <- ct.run2$`2421070` [which(ct.run2$`2421070`$SNP %in% 
                                                naive.intersecting$`444332`$SNP),]
ct.run2$`2447896` <- ct.run2$`2447896` [which(ct.run2$`2447896`$SNP %in% 
                                                naive.intersecting$`449650`$SNP),]

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

ct.all <- as.data.frame(bind_rows(ct.table))
ct.all$Ct <- as.double(ct.all$Ct)
ct.all <- ct.all %>%
  arrange(Ct, desc = T) %>%
  mutate("PARTICIPANT" = c(432870,432870,449650,449650,435786,435805,435805)) %>%
  relocate(PARTICIPANT, .before = Ct) %>%
  mutate("DAY" = c(1,1,1,1,3,9,9)) %>%
  relocate(DAY, .after = PARTICIPANT)
cor.test(ct.all$RUN1, ct.all$RUN2, method = "p") #cor = 0.8419865, p = 0.01747
export(ct.all, file = "runcontrol_ct.csv")

user.table.cut <- user.table[,-2]
ct.all.cut <- ct.all[,-c(1:3)]
runcontrol <- bind_rows(user.table.cut,ct.all.cut)
cor.test(runcontrol$RUN1, runcontrol$RUN2) #cor = 0.847303 #p = 0.0009929
save(runcontrol, file = "runcontrol.RData")




