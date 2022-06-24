library(tidyverse)
library(rio)
library(rlist)
library(gdata)
library(gtools)
library(here)

dirlist <- list.dirs("naive_output_tables")[-c(1:3)]

#pare down data to necessary elements
sleek <- function(dir){
  files <- list.files(dir,pattern = "ivar")
  files <- gtools::mixedsort(sort(files))
  full <- paste0(dir,"/",files)
  user <- list()
  for(i in 1:length(full)){
    user[[i]] <- import(full[i])
  }
  sleekuser <- list()
  for(i in 1:length(user)){
    sleekuser[[i]] <- user[[i]] %>%
      filter(ALT_FREQ>=.03,ALT_FREQ<=.97)%>%
      filter(ALT_DP + REF_DP >=1000)%>%
      select(POS,REF,ALT,ALT_DP,REF_DP,ALT_FREQ)%>%
      distinct()
  }
  return(sleekuser)
} #POS,REF,ALT,ALT_DP,REF_DP,ALT_FREQ + data thresholds
everything <- function(dir){
  files <- list.files(dir,pattern = "ivar")
  files <- gtools::mixedsort(sort(files))
  full <- paste0(dir,"/",files)
  user <- list()
  for(i in 1:length(full)){
    user[[i]] <- import(full[i])
  }
  everythinguser <- list()
  for(i in 1:length(user)){
    everythinguser[[i]] <- user[[i]] %>%
      select(POS,REF,ALT,ALT_FREQ,TOTAL_DP)%>%
      distinct()
  }
  return(everythinguser)
} #POS,REF,ALT,ALT_FREQ,TOTAL_DP + no thresholds
eff <- function(dir){
  files <- list.files(dir,pattern = "snp")
  files <- gtools::mixedsort(sort(files))
  full <- paste0(dir,"/",files)
  user <- lapply(full,import)
  effuser <- list()
  for(i in seq_along(user)){
    effuser[[i]] <- as_tibble(user[[i]]) %>%
      distinct()%>%
      dplyr::rename("GENE" = `ANN[*].GENE`)%>%
      dplyr::rename("IMP" = `ANN[*].IMPACT`)%>%
      dplyr::rename("EFF" =  `EFF[*].EFFECT`)%>%
      dplyr::rename("AA" = `EFF[*].AA`)%>%
      filter(IMP != "MODIFIER")%>%
      filter(IMP != "MODIFIER,MODIFIER")%>%
      mutate("EFF.CUT" = strsplit(EFF,","))%>%
      mutate("GENE.CUT"=strsplit(GENE,","))
  }
  for(i in seq_along(effuser)){
    for(j in seq_along(effuser[[i]]$POS)){
      effuser[[i]]$EFF[[j]] <- effuser[[i]]$EFF.CUT[[j]][1]
    }
  }
  
  for(i in seq_along(effuser)){
    for(j in seq_along(effuser[[i]]$POS)){
      effuser[[i]]$GENE[[j]] <- effuser[[i]]$GENE.CUT[[j]][1]
    }
  }
  
  for(i in seq_along(effuser)){
    for(j in seq_along(effuser[[i]]$EFF)){
      if(effuser[[i]]$EFF[[j]] == "downstream_gene_variant"){
        effuser[[i]]$GENE[[j]] <- "UTR"
      }
      if(effuser[[i]]$EFF[[j]] == "upstream_gene_variant"){
        effuser[[i]]$GENE[[j]] <- "UTR"
      }
      if(is.na(effuser[[i]]$EFF[[j]])){
        effuser[[i]]$GENE[[j]] <- "UTR"
      }
    }
  }
  
  for(i in seq_along(effuser)){
    effuser[[i]] <- effuser[[i]] %>%
      select(POS,REF,ALT,GENE,EFF,AA)
  }
  return(effuser)
} #gene annotation files

sleekuser <- lapply(dirlist,sleek)
everythinguser <- lapply(dirlist,everything)
effuser <- lapply(dirlist,eff)

names(sleekuser) <- c("432686","432870","433227","435786","435805","438577",
                      "442978","444332","444446","444633","445602","449650",
                      "450241","450348","451152","451709","453058","459597",
                      "471467","471588")
names(everythinguser) <- names(sleekuser)
names(effuser) <- names(sleekuser)
names(dirlist) <- names(sleekuser)

#count number of SNPs for each user on each day
snpcount <- function(user){
  snplength <-list()
  for (i in seq_along(user)){
    user[[i]] <- user[[i]] %>%
      filter(POS != 6696)%>% #these are likely sequencing artifacts
      filter(POS!= 11074)%>%
      filter(POS != 15965)%>%
      filter(POS!= 29051)
    snplength[[i]] <- length(user[[i]]$POS)
  }
  return(snplength)
}
snpcounts <- lapply(sleekuser,snpcount)

snpcounts_dat <- data.frame(matrix(ncol = 20,nrow = 10))
rownames(snpcounts_dat) <- paste("time point",1:10)
colnames(snpcounts_dat) <- names(dirlist)
for(i in seq_along(snpcounts)){
  for(j in seq_along(snpcounts[[i]])){
    snpcounts_dat[j,i] <- snpcounts[[i]][[j]]
  }
}

save(snpcounts_dat, file = "snpcounts_dat.RData")
write.csv(snpcounts_dat,"naive_snpcounts.csv")

#pull out SNP positions
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
allpos <- lapply(sleekuser,positions)

#pull out variants that do NOT appear more than once within users
keepunique <- function(pos){
  splat <- list.rbind(pos)
  unique  <- distinct(splat)
  unique <- unique %>%
    filter(POS != 6696)%>%
    filter(POS != 11074)%>%
    filter(POS != 15965)%>%
    filter(POS != 29051)
  return(unique)
}
uniqueSNPs <- list()
uniquecounts <- list()
for(i in seq_along(allpos)){
  uniqueSNPs[[i]] <- keepunique(allpos[[i]])
  uniquecounts[[i]] <- length(uniqueSNPs[[i]]$POS)
}

#pull out intersecting iSNVs
var1 <- function(allpos){
  time1 <- list()
  for(i in 1:length(allpos)){
    time1[[i]] <- dplyr::intersect(allpos[[1]],allpos[[i]])
  }
  return(time1)
}
var2 <- function(allpos){
  time2 <- list()
  for(i in 1:length(allpos)){
    time2[[i]] <- dplyr::intersect(allpos[[2]],allpos[[i]])
  }
  return(time2)
}
var3 <- function(allpos){
  time3 <- list()
  for(i in 1:length(allpos)){
    time3[[i]] <- dplyr::intersect(allpos[[3]],allpos[[i]])
  }
  return(time3)
}
var4 <- function(allpos){
  time4 <- list()
  for(i in 1:length(allpos)){
    time4[[i]] <- dplyr::intersect(allpos[[4]],allpos[[i]])
  }
  return(time4)
}
var5 <- function(allpos){
  time5 <- list()
  for(i in 1:length(allpos)){
    time5[[i]] <- dplyr::intersect(allpos[[5]],allpos[[i]])
  }
  return(time5)
}
var6 <- function(allpos){
  time6 <- list()
  for(i in 1:length(allpos)){
    time6[[i]] <- dplyr::intersect(allpos[[6]],allpos[[i]])
  }
  return(time6)
}
var7 <- function(allpos){
  time7 <- list()
  for(i in 1:length(allpos)){
    time7[[i]] <- dplyr::intersect(allpos[[7]],allpos[[i]])
  }
  return(time7)
}
var8 <- function(allpos){
  time8 <- list()
  for(i in 1:length(allpos)){
    time8[[i]] <- dplyr::intersect(allpos[[8]],allpos[[i]])
  }
  return(time8)
}
var9 <- function(allpos){
  time9 <- list()
  for(i in 1:length(allpos)){
    time9[[i]] <- dplyr::intersect(allpos[[9]],allpos[[i]])
  }
  return(time9)
}
var10 <- function(allpos){
  time10 <- list()
  for(i in 1:length(allpos)){
    time10[[i]] <- dplyr::intersect(allpos[[10]],allpos[[i]])
  }
  return(time10)
}


var_432686 <- list("time1"=var1(allpos[[1]]),"time2"=var2(allpos[[1]]),
                   "time3"=var3(allpos[[1]]),"time4"=var4(allpos[[1]]),
                   "time5"=var5(allpos[[1]]),"time6"=var6(allpos[[1]]))
var_432870 <- list("time1"=var1(allpos[[2]]),"time2"=var2(allpos[[2]]),
                   "time3"=var3(allpos[[2]]),"time4"=var4(allpos[[2]]),
                   "time5"=var5(allpos[[2]]),"time6"=var6(allpos[[2]]),
                   "time7"=var7(allpos[[2]]),"time8"=var8(allpos[[2]]),
                   "time9"=var9(allpos[[2]]))
var_433227 <- list("time1"=var1(allpos[[3]]),"time2"=var2(allpos[[3]]),
                   "time3"=var3(allpos[[3]]),"time4"=var4(allpos[[3]]),
                   "time5"=var5(allpos[[3]]),"time6"=var6(allpos[[3]]),
                   "time7"=var7(allpos[[3]]))
var_435786 <- list("time1"=var1(allpos[[4]]),"time2"=var2(allpos[[4]]),
                   "time3"=var3(allpos[[4]]),"time4"=var4(allpos[[4]]),
                   "time5"=var5(allpos[[4]]),"time6"=var6(allpos[[4]]))
var_435805 <- list("time1"=var1(allpos[[5]]),"time2"=var2(allpos[[5]]),
                   "time3"=var3(allpos[[5]]),"time4"=var4(allpos[[5]]),
                   "time5"=var5(allpos[[5]]),"time6"=var6(allpos[[5]]),
                   "time7"=var7(allpos[[5]]),"time8"=var8(allpos[[5]]),
                   "time9"=var9(allpos[[5]]))
var_438577 <- list("time1"=var1(allpos[[6]]),"time2"=var2(allpos[[6]]),
                   "time3"=var3(allpos[[6]]),"time4"=var4(allpos[[6]]),
                   "time5"=var5(allpos[[6]]),"time6"=var6(allpos[[6]]))
var_442978 <- list("time1"=var1(allpos[[7]]),"time2"=var2(allpos[[7]]),
                   "time3"=var3(allpos[[7]]),"time4"=var4(allpos[[7]]),
                   "time5"=var5(allpos[[7]]))
var_444332 <- list("time1"=var1(allpos[[8]]),"time2"=var2(allpos[[8]]),
                   "time3"=var3(allpos[[8]]),"time4"=var4(allpos[[8]]),
                   "time5"=var5(allpos[[8]]),"time6"=var6(allpos[[8]]),
                   "time7"=var7(allpos[[8]]),"time8"=var8(allpos[[8]]))
var_444446 <- list("time1"=var1(allpos[[9]]),"time2"=var2(allpos[[9]]),
                   "time3"=var3(allpos[[9]]),"time4"=var4(allpos[[9]]),
                   "time5"=var5(allpos[[9]]),"time6"=var6(allpos[[9]]),
                   "time7"=var7(allpos[[9]]),"time8"=var8(allpos[[9]]),
                   "time9"=var9(allpos[[9]]))
var_444633 <- list("time1"=var1(allpos[[10]]),"time2"=var2(allpos[[10]]),
                   "time3"=var3(allpos[[10]]),"time4"=var4(allpos[[10]]),
                   "time5"=var5(allpos[[10]]))
var_445602 <- list("time1"=var1(allpos[[11]]),"time2"=var2(allpos[[11]]),
                   "time3"=var3(allpos[[11]]),"time4"=var4(allpos[[11]]),
                   "time5"=var5(allpos[[11]]),"time6"=var6(allpos[[11]]),
                   "time7"=var7(allpos[[11]]),"time8"=var8(allpos[[11]]),
                   "time9"=var9(allpos[[11]]),"time10"=var10(allpos[[11]]))
var_449650 <- list("time1"=var1(allpos[[12]]),"time2"=var2(allpos[[12]]),
                   "time3"=var3(allpos[[12]]),"time4"=var4(allpos[[12]]),
                   "time5"=var5(allpos[[12]]),"time6"=var6(allpos[[12]]),
                   "time7"=var7(allpos[[12]]),"time8"=var8(allpos[[12]]))
var_450241 <- list("time1"=var1(allpos[[13]]),"time2"=var2(allpos[[13]]),
                   "time3"=var3(allpos[[13]]),"time4"=var4(allpos[[13]]),
                   "time5"=var5(allpos[[13]]),"time6"=var6(allpos[[13]]))
var_450348 <- list("time1"=var1(allpos[[14]]),"time2"=var2(allpos[[14]]),
                   "time3"=var3(allpos[[14]]),"time4"=var4(allpos[[14]]),
                   "time5"=var5(allpos[[14]]),"time6"=var6(allpos[[14]]),
                   "time7"=var7(allpos[[14]]),"time8"=var8(allpos[[14]]),
                   "time9"=var9(allpos[[14]]))
var_451152 <- list("time1"=var1(allpos[[15]]),"time2"=var2(allpos[[15]]),
                   "time3"=var3(allpos[[15]]),"time4"=var4(allpos[[15]]),
                   "time5"=var5(allpos[[15]]),"time6"=var6(allpos[[15]]),
                   "time7"=var7(allpos[[15]]),"time8"=var8(allpos[[15]]))
var_451709 <- list("time1"=var1(allpos[[16]]),"time2"=var2(allpos[[16]]),
                   "time3"=var3(allpos[[16]]),"time4"=var4(allpos[[16]]),
                   "time5"=var5(allpos[[16]]),"time6"=var6(allpos[[16]]),
                   "time7"=var7(allpos[[16]]),"time8"=var8(allpos[[16]]),
                   "time9"=var9(allpos[[16]]))
var_453058 <- list("time1"=var1(allpos[[17]]),"time2"=var2(allpos[[17]]),
                   "time3"=var3(allpos[[17]]),"time4"=var4(allpos[[17]]),
                   "time5"=var5(allpos[[17]]),"time6"=var6(allpos[[17]]),
                   "time7"=var7(allpos[[17]]))
var_459597 <- list("time1"=var1(allpos[[18]]),"time2"=var2(allpos[[18]]),
                   "time3"=var3(allpos[[18]]),"time4"=var4(allpos[[18]]),
                   "time5"=var5(allpos[[18]]),"time6"=var6(allpos[[18]]),
                   "time7"=var7(allpos[[18]]))
var_471467 <- list("time1"=var1(allpos[[19]]),"time2"=var2(allpos[[19]]),
                   "time3"=var3(allpos[[19]]),"time4"=var4(allpos[[19]]),
                   "time5"=var5(allpos[[19]]),"time6"=var6(allpos[[19]]),
                   "time7"=var7(allpos[[19]]),"time8"=var8(allpos[[19]]),
                   "time9"=var9(allpos[[19]]),"time10"=var10(allpos[[19]]))
var_471588 <- list("time1"=var1(allpos[[20]]),"time2"=var2(allpos[[20]]),
                   "time3"=var3(allpos[[20]]),"time4"=var4(allpos[[20]]),
                   "time5"=var5(allpos[[20]]),"time6"=var6(allpos[[20]]),
                   "time7"=var7(allpos[[20]]),"time8"=var8(allpos[[20]]),
                   "time9"=var9(allpos[[20]]))

var_list <- list("var_432686"=var_432686,"var_432870"=var_432870,
                 "var_433227"=var_433227,"var_435786"=var_435786,
                 "var_435805"=var_435805,"var_438577"=var_438577,
                 "var_442978"=var_442978,"var_444332"=var_444332,
                 "var_444446"=var_444446,"var_444633"=var_444633,
                 "var_445602"=var_445602,"var_449650"=var_449650,
                 "var_450241"=var_450241,"var_450348"=var_450348,
                 "var_451152"=var_451152,"var_451709"=var_451709,
                 "var_453058"=var_453058,"var_459597"=var_459597,
                 "var_471467"=var_471467,"var_471588"=var_471588)

#remove self-comparison
selfevict <- function(var){
  for(i in 1:length(var)){
    var[[i]][i] <- NA
  }
  return(var)
}
var_list <- lapply(var_list,selfevict)

#write file for each user to a .csv
var_filenames <- paste0("naive_intersect/",names(var_list),".csv")
for(i in seq_along(var_list)){
  capture.output(var_list[[i]],file=var_filenames[[i]])
}

#pull out a list of mutations that appear > once for each user (shared iSNVs)
keepcommon <- function(var){
  flat <- flatten(var)
  splat <- na.omit(rlist::list.rbind(flat))
  common <- distinct(splat)
  return(common)
}
naive_intersecting <- lapply(var_list,keepcommon)
names(naive_intersecting) <- gsub("var","int",names(var_list))
naive_intersecting <- lapply(naive_intersecting,arrange,POS)
save(naive_intersecting,file = "naive_intersecting.RData")

#count shared iSNVs present on each day
shared_lengths <- function(all,intersecting){
  all_cut <- list()
  for(i in seq_along(all)){
    all_cut[[i]] <- all[[i]] %>%
      filter(POS != 6696) %>%
      filter(POS != 11074) %>%
      filter(POS != 15965) %>%
      filter(POS!= 29051)
  }
  subset <- list()
  subset_lengths <- list()
  for(i in seq_along(all_cut)){
    subset[[i]] <- which(all_cut[[i]]$POS %in% intersecting$POS)
    subset_lengths[[i]] <- length(subset[[i]])
  }
  subset_lengths <- unlist(subset_lengths)
  return(subset_lengths)
}
naive_daily_shared <- list()
for(i in seq_along(naive_intersecting)){
  naive_daily_shared[[i]] <- shared_lengths(allpos[[i]],naive_intersecting[[i]])
}
save(naive_daily_shared, file = "naive_daily_shared.RData")

#create variant tables to track shared iSNV frequencies over time
varianttable <-function(intersectuser,everythinguser){
  freq <- list()
  for(i in 2:(length(everythinguser)+1)){
    freq[[1]] <- intersectuser
    freq[[i]] <- everythinguser[[i-1]]
  }
  joiner <- freq %>%
    purrr::reduce(left_join, by=c("POS","REF","ALT"))
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
      if(strsplit(colnames(joiner)[[j]],"_")[[1]][1] == "depth"){
        if((!is.na(joiner[i,j])) && (joiner[i,j]<1000)){
          joiner[i,(j-1)] <- paste(joiner[i,(j-1)],"(low dp)")
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
naive.vartables <- list()
for(i in seq_along(naive_intersecting)){
  naive.vartables[[i]] <- varianttable(naive_intersecting[[i]],everythinguser[[i]])
}
names(naive.vartables) <- names(sleekuser)

#replace NA with 0.01 (for plotting on log scale)
NA_01 <- function(vartable){
  for(i in 1:length(vartable$POS)){
    for(j in 1:length(colnames(vartable))){
      if(is.na(vartable[i,j])){
        vartable[i,j] <- 0.01
      }
    }
  }
  return(vartable)
}
naive.vartables <- lapply(naive.vartables,NA_01)

#reformat SNP notation with no +/-'s (eg. REF = C, ALT= +T --> REF = C, ALT = CT)
for(i in seq_along(naive.vartables)){
  for(j in seq_along(naive.vartables[[i]]$ALT)){
    if(substr(naive.vartables[[i]]$ALT[[j]],0,1)=="+"){
      naive.vartables[[i]]$ALT[[j]] <- paste0(naive.vartables[[i]]$REF[[j]], 
                                               substr(naive.vartables[[i]]$ALT[[j]],2,nchar(naive.vartables[[i]]$ALT[[j]])))
    }
    if(substr(naive.vartables[[i]]$ALT[[j]],0,1)=="-"){
      naive.vartables[[i]]$REF[[j]] <- paste0(naive.vartables[[i]]$REF[[j]],
                                               substr(naive.vartables[[i]]$ALT[[j]],2,nchar(naive.vartables[[i]]$ALT[[j]])))
      naive.vartables[[i]]$ALT[[j]] <- substr(naive.vartables[[i]]$REF[[j]],0,1)
    }
  }
}

#annotate intersecting iSNVs with gene functions
annotable <-function(vartables,effuser){
  annotable <- vartables%>%
    mutate("GENE"=NA)%>%
    mutate("EFF"=NA)%>%
    mutate("AA"=NA)%>%
    select(POS,REF,ALT,GENE,EFF,AA)
  for(x in 1:length(annotable$POS)){
    for(i in 1:length(effuser)){
      for(j in 1:length(effuser[[i]]$POS)){
        if((effuser[[i]]$POS[[j]] == annotable$POS[[x]]) &&
           (effuser[[i]]$REF[[j]] == annotable$REF[[x]]) &&
           (effuser[[i]]$ALT[[j]] == annotable$ALT[[x]])){
          annotable$GENE[[x]] <- effuser[[i]]$GENE[[j]]
          annotable$EFF[[x]] <- effuser[[i]]$EFF[[j]]
          annotable$AA[[x]] <- effuser[[i]]$AA[[j]]
        }
      }
    }
  }
  return(annotable)
} 
naive.eff <- list()
for(i in seq_along(naive.vartables)){
  naive.eff[[i]] <- annotable(naive.vartables[[i]],effuser[[i]])
}
names(naive.eff) <- names(effuser)

#annotate UTR as 5' or 3' UTR
for(i in seq_along(naive.eff)){
  for(j in seq_along(naive.eff[[i]]$POS)){
    if(naive.eff[[i]]$POS[[j]] <= 265){
      naive.eff[[i]]$GENE[[j]] <- "5'UTR"
    }
    if(naive.eff[[i]]$POS[[j]] >= 29675){
      naive.eff[[i]]$GENE[[j]] <- "3'UTR"
    }
  }
}

#add columns for ID and reformat to REF-POS-ALT SNP notation
for(i in seq_along(naive.eff)){
  naive.eff[[i]] <- naive.eff[[i]] %>%
    mutate("ID" = names(naive.eff)[[i]]) %>%
    relocate(ID,.before = POS) %>%
    mutate("SNP"=paste0(REF,POS,ALT))%>%
    relocate(SNP,.after = ID)%>%
    select(ID,SNP,POS,REF,ALT,GENE,EFF,AA)
}

#reformat EFF column to read synonymous/nonsynonymous/utr
for(i in seq_along(naive.eff)){
  for(j in seq_along(naive.eff[[i]]$GENE)){
    if((naive.eff[[i]]$GENE[[j]] == "5'UTR") |
       (naive.eff[[i]]$GENE[[j]] == "3'UTR")){
      naive.eff[[i]]$EFF[[j]] <- "UTR"
    }
    if(naive.eff[[i]]$EFF[[j]] == "synonymous_variant"){
      naive.eff[[i]]$EFF[[j]] <- "synonymous"
    }
  }
}
for(i in seq_along(naive.eff)){
  for(j in seq_along(naive.eff[[i]]$EFF)){
    if((naive.eff[[i]]$EFF[[j]] != "UTR") &&
       (naive.eff[[i]]$EFF[[j]] != "synonymous")){
      naive.eff[[i]]$EFF[[j]] <- "nonsynonymous"
    }
  }
}

#get rid of extraneous punctuation in AA column
for(i in seq_along(naive.eff)){
  for(j in seq_along(naive.eff[[i]]$AA)){
    naive.eff[[i]]$AA[[j]] <- strsplit(naive.eff[[i]]$AA[[j]],",")[[1]][1]
  }
}

#filter out seq artefacts
for(i in seq_along(naive.eff)){
  naive.eff[[i]] <- naive.eff[[i]] %>%
    filter(POS != 6696) %>%
    filter(POS != 11074) %>%
    filter(POS != 15965) %>%
    filter(POS!= 29051)
}

#write individual annotation csv files for each user
for(i in seq_along(naive.eff)){
  filenames <- paste0("naive_annotations/annotated_",names(naive.eff),".csv")
  write.csv(naive.eff[[i]],filenames[[i]])
}

#bind annotations for all users into one large table
naive.annotable <- bind_rows(naive.eff)
write.csv(naive.annotable,"naive_annotations.csv")

#reformat variant plots to REF-POS-ALT SNP notation
for(i in seq_along(naive.vartables)){
  for(j in seq_along(naive.vartables[[i]]$POS)){
    naive.vartables[[i]]$POS[[j]] <- paste0(naive.vartables[[i]]$REF[[j]], 
      naive.vartables[[i]]$POS[[j]], naive.vartables[[i]]$ALT[[j]])
  }
}

for(i in seq_along(naive.vartables)){
  naive.vartables[[i]] <- naive.vartables[[i]] %>%
    rename("SNP"=POS)%>%
    select(!REF)%>%
    select(!ALT)
}

#transpose data (makes things easier to plot)
for(i in seq_along(naive.vartables)){
  naive.vartables[[i]] <- as.data.frame(t(naive.vartables[[i]]))
  colnames(naive.vartables[[i]]) <- naive.vartables[[i]][1,]
  naive.vartables[[i]] <- naive.vartables[[i]][-1,]
  rownames(naive.vartables[[i]]) <- paste0("freq_",1:dim(naive.vartables[[i]])[[1]])
}

#write csvs for variant frequency tracking tables
for(i in seq_along(naive.vartables)){
  filenames <- paste0("naive_variant_tables/user_",names(naive.vartables),".csv")
  write.csv(naive.vartables[[i]],filenames[[i]])
}

#count the total number of intersecting mutations for each user
load("naive_intersecting.RData")
naive_intersecting_cut <- list()
naive_intersect_lengths <- list()
for(i in seq_along(naive_intersecting)){
  naive_intersecting_cut[[i]] <- naive_intersecting[[i]] %>%
    filter(POS!= 6696)%>%
    filter(POS!= 11074)%>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)
  naive_intersect_lengths[[i]] <- length(naive_intersecting_cut[[i]]$POS)
}

#calculate percentage of intersecting (shared) mutants for each user
naive_intersect_percent <- list()
for(i in seq_along(uniquecounts)){
  naive_intersect_percent[[i]] <- naive_intersect_lengths[[i]]/(naive_intersect_lengths[[i]]+uniquecounts[[i]])
}

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
naive_shared <- as.data.frame(matrix(nrow = 20,ncol = 4))
colnames(naive_shared) <- c("Participant ID", "Shared", "Unique", "Total")
naive_shared[,1] <- names(dirlist)
naive_shared[,2] <- as.numeric(unlist(naive_intersect_lengths))
naive_shared[,3] <- as.numeric(unlist(uniquecounts))
naive_shared[,4] <- naive_shared[,2] + naive_shared[,3]
naive_shared$`Participant ID` <- as.character(naive_shared$`Participant ID`)
save(naive_shared, file = "naive_shared.RData")

#plot daily iSNV counts (FIG 2A)
load("snpcounts_dat.Rdata")
naive_count_gather <- snpcounts_dat%>%
  mutate("day"=1:10)%>%
  gather(key = "Participant ID",value="SNP_count",-day)

twenty.palette <-c("#9D6A90","#94719A","#8978A2","#7C7FA9","#6D86AD","#5D8DAF",
                   "#4C93AF","#3C99AC","#2E9EA7","#26A39F","#2AA796","#36AB8C",
                   "#46AE81","#58B075","#6AB269","#7CB25E","#8FB355","#A2B24D",
                   "#B5B148","#C8AF46")

ggplot(data=naive_count_gather,aes(x=`Participant ID`,y=SNP_count,fill = `Participant ID`))+
  geom_dotplot(binaxis = "y", binwidth = .08, stackdir = "center", color = NA)+ #binwidth .08
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  scale_y_log10()+
  geom_col(data = naive_shared,aes(x=`Participant ID`,y=Shared), 
           fill=NA,color = "black")+
  geom_col(data=naive_shared,aes(x=`Participant ID`,y=Total),
           fill=NA,color="grey40")+
  scale_fill_manual(values = twenty.palette)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size = 19),
        legend.title = element_text(size = 18),legend.text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle=50, margin = margin(t=3), 
                                   hjust = 1),
        axis.title.x = element_text(margin = margin(t=10)),
        plot.title = element_text(size = 22))
ggsave("figs/naive_SNPs_shared_total_box.png")
save(naive_count_gather,file="naive_count_gather.Rdata")






