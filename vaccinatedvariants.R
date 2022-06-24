library(tidyverse)
library(rio)
library(rlist)
library(gtools)
library(here)

dirlist_vax <- list.dirs("vaccinated_output_tables")[-c(1,10)]

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
} #POS,REF,ALT,ALT_DP,REF_DP,ALT_FREQ
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
}  #gene annotation files

sleekuser_vax <- lapply(dirlist_vax,sleek)
names(sleekuser_vax) <- c("461913","471876","475670","475951","481242","481672",
                          "482828","484249","486422","487250","487297",
                          "487941")
names(dirlist_vax) <- names(sleekuser_vax)
everythinguser_vax <- lapply(dirlist_vax,everything)
names(everythinguser_vax) <- names(sleekuser_vax)
effuser_vax <- lapply(dirlist_vax,eff)
names(effuser_vax) <- names(sleekuser_vax)

#count snps
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
snpcounts_vax <- lapply(sleekuser_vax,snpcount)

snpcounts_dat_vax <- data.frame(matrix(ncol = 12,nrow = 8))
rownames(snpcounts_dat_vax) <- paste("time point",1:8)
colnames(snpcounts_dat_vax) <- names(dirlist_vax)
for(i in seq_along(snpcounts_vax)){
  for(j in seq_along(snpcounts_vax[[i]])){
    snpcounts_dat_vax[j,i] <- snpcounts_vax[[i]][[j]]
  }
}

save(snpcounts_dat_vax, file = "snpcounts_dat_vax.RData")
write.csv(snpcounts_dat_vax,"vaccinated_snpcounts.csv")

#pull out SNP positions
positions <- function(sleekuser){
  allpos <- list()
  for (i in 1:length(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
allpos_vax <- lapply(sleekuser_vax,positions)
names(allpos_vax) <- names(sleekuser_vax)

#pull out SNPs that do NOT appear more than once in a single user
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
uniqueSNPs_vax <- list()
uniquecounts_vax <- list()
for(i in 1:length(allpos_vax)){
  uniqueSNPs_vax[[i]] <- keepunique(allpos_vax[[i]])
  uniquecounts_vax[[i]] <- length(uniqueSNPs_vax[[i]]$POS)
}
names(uniquecounts_vax) <- names(sleekuser_vax)
uniquecounts_vax <- uniquecounts_vax[-c(1,6,12)] #take out users w/ only 1 timept

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

#only includes users with more than one time point
var_471876 <- list("time1"=var1(allpos_vax[[2]]),"time2"=var2(allpos_vax[[2]]),
                   "time3"=var3(allpos_vax[[2]]),"time4"=var4(allpos_vax[[2]]),
                   "time5"=var5(allpos_vax[[2]]),"time6"=var6(allpos_vax[[2]]),
                   "time7"=var7(allpos_vax[[2]]),"time8"=var8(allpos_vax[[2]]))
var_475670 <- list("time1"=var1(allpos_vax[[3]]),"time2"=var2(allpos_vax[[3]]),
                   "time3"=var3(allpos_vax[[3]]),"time4"=var4(allpos_vax[[3]]),
                   "time5"=var5(allpos_vax[[3]]),"time6"=var6(allpos_vax[[3]]),
                   "time7"=var7(allpos_vax[[3]]))
var_475951 <- list("time1"=var1(allpos_vax[[4]]),"time2"=var2(allpos_vax[[4]]),
                   "time3"=var3(allpos_vax[[4]]),"time4"=var4(allpos_vax[[4]]))
var_481242 <- list("time1"=var1(allpos_vax[[5]]),"time2"=var2(allpos_vax[[5]]),
                   "time3"=var3(allpos_vax[[5]]),"time4"=var4(allpos_vax[[5]]),
                   "time5"=var5(allpos_vax[[5]]),"time6"=var6(allpos_vax[[5]]),
                   "time7"=var7(allpos_vax[[5]]))
var_482828 <- list("time1"=var1(allpos_vax[[7]]),"time2"=var2(allpos_vax[[7]]),
                   "time3"=var3(allpos_vax[[7]]),"time4"=var4(allpos_vax[[7]]))
var_484249 <- list("time1"=var1(allpos_vax[[8]]),"time2"=var2(allpos_vax[[8]]))
var_486422 <- list("time1"=var1(allpos_vax[[9]]),"time2"=var2(allpos_vax[[9]]),
                   "time3"=var3(allpos_vax[[9]]))
var_487250 <- list("time1"=var1(allpos_vax[[10]]),"time2"=var2(allpos_vax[[10]]),
                   "time3"=var3(allpos_vax[[10]]),"time4"=var4(allpos_vax[[10]]))
var_487297 <- list("time1"=var1(allpos_vax[[11]]),"time2"=var2(allpos_vax[[11]]),
                   "time3"=var3(allpos_vax[[11]]))
var_list_vax <- list("var_471876"=var_471876,"var_475670"=var_475670,
                 "var_475951"=var_475951,"var_481242"=var_481242,
                 "var_482828"=var_482828,"var_484249"=var_484249,
                 "var_486422"=var_486422,"var_487250"=var_487250,
                 "var_487297"=var_487297)
            
#remove self-comparison
selfevict <- function(var){
  for(i in 1:length(var)){
    var[[i]][i] <- NA
  }
  return(var)
}
for(i in 1:length(var_list_vax)){
  var_list_vax[[i]] <- selfevict(var_list_vax[[i]])
}

#write file for each user to a .csv
var_filenames_vax <- paste0("vaccinated_intersect/",names(var_list_vax),".csv")
for(i in 1:length(var_list_vax)){
  capture.output(var_list_vax[[i]],file=var_filenames_vax[[i]])
}

#pull out a list of mutations that appear > once for each user (shared iSNVs)
keepcommon <- function(var){
  flat <- flatten(var)
  splat <- na.omit(list.rbind(flat))
  common <- distinct(splat)
  return(common)
}
vax_intersecting <- lapply(var_list_vax,keepcommon)
names(vax_intersecting) <- gsub("var","int",names(var_list_vax))
vax_intersecting <- lapply(vax_intersecting,arrange,POS)
save(vax_intersecting,file = "vax_intersecting.RData")

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
allpos_vax_cut <- allpos_vax[-c(1,6,12)]
vax_daily_shared <- list()
for(i in seq_along(vax_intersecting)){
  vax_daily_shared[[i]] <- shared_lengths(allpos_vax_cut[[i]],vax_intersecting[[i]])
}
save(vax_daily_shared, file = "vax_daily_shared.RData")

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
everythinguser_vax_cut <- everythinguser_vax[-c(1,6,12)]
vax.vartables <- list()
for(i in 1:length(vax_intersecting)){
  vax.vartables[[i]] <- varianttable(vax_intersecting[[i]],everythinguser_vax_cut[[i]])
}
names(vax.vartables) <- names(vax_intersecting)

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
for(i in 1:length(vax.vartables)){
  vax.vartables[[i]] <- NA_01(vax.vartables[[i]])
}

#reformat SNP notation with no +/-'s (eg. REF = C, ALT= +T --> REF = C, ALT = CT)
for(i in 1:length(vax.vartables)){
  for(j in 1:length(vax.vartables[[i]]$ALT)){
    if(substr(vax.vartables[[i]]$ALT[[j]],0,1)=="+"){
      vax.vartables[[i]]$ALT[[j]] <- paste0(vax.vartables[[i]]$REF[[j]], 
                                              substr(vax.vartables[[i]]$ALT[[j]],2,nchar(vax.vartables[[i]]$ALT[[j]])))
    }
    if(substr(vax.vartables[[i]]$ALT[[j]],0,1)=="-"){
      vax.vartables[[i]]$REF[[j]] <- paste0(vax.vartables[[i]]$REF[[j]],
                                              substr(vax.vartables[[i]]$ALT[[j]],2,nchar(vax.vartables[[i]]$ALT[[j]])))
      vax.vartables[[i]]$ALT[[j]] <- substr(vax.vartables[[i]]$REF[[j]],0,1)
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

effuser_vax_cut <- effuser_vax[-c(1,6,12)] #remove users w/ 1 timept
vax.eff <- list()
for(i in 1:length(vax.vartables)){
  vax.eff[[i]] <- annotable(vax.vartables[[i]],effuser_vax_cut[[i]])
}

#annotate UTR as 5'/3'/general UTR
for(i in 1:length(vax.eff)){
  for(j in 1:length(vax.eff[[i]]$POS)){
    if(vax.eff[[i]]$POS[[j]] <= 265){
      vax.eff[[i]]$GENE[[j]] <- "5'UTR"
    }
    if((vax.eff[[i]]$POS[[j]]>28259) &&
       (vax.eff[[i]]$POS[[j]]<28274)){
      vax.eff[[i]]$GENE[[j]] <- "UTR"
    }
    if(vax.eff[[i]]$POS[[j]] >= 29675){
      vax.eff[[i]]$GENE[[j]] <- "3'UTR"
    }
  }
}

#add columns for ID and reformat to REF-POS-ALT SNP notation
names(vax.eff) <- gsub("int_","",names(vax.vartables))
for(i in 1:length(vax.eff)){
  vax.eff[[i]] <- vax.eff[[i]] %>%
    mutate("ID" = names(vax.eff)[[i]]) %>%
    relocate(ID,.before = POS) %>%
    mutate("SNP"=paste0(REF,POS,ALT))%>%
    relocate(SNP,.after = ID)%>%
    select(ID,SNP,POS,REF,ALT,GENE,EFF,AA)
}

#reformat EFF column to read synonymous/nonsynonymous/utr
for(i in 1:length(vax.eff)){
  for(j in 1:length(vax.eff[[i]]$GENE)){
    if((vax.eff[[i]]$GENE[[j]] == "5'UTR") |
       (vax.eff[[i]]$GENE[[j]] == "3'UTR") |
       (vax.eff[[i]]$GENE[[j]] == "UTR")){
      vax.eff[[i]]$EFF[[j]] <- "UTR"
    }
    if(vax.eff[[i]]$EFF[[j]] == "synonymous_variant"){
      vax.eff[[i]]$EFF[[j]] <- "synonymous"
    }
  }
}
for(i in 1:length(vax.eff)){
  for(j in 1:length(vax.eff[[i]]$EFF)){
    if((vax.eff[[i]]$EFF[[j]] != "UTR") &&
       (vax.eff[[i]]$EFF[[j]] != "synonymous")){
      vax.eff[[i]]$EFF[[j]] <- "nonsynonymous"
    }
  }
}

#get rid of extraneous punctuation in AA column
for(i in 1:length(vax.eff)){
  for(j in 1:length(vax.eff[[i]]$AA)){
    vax.eff[[i]]$AA[[j]] <- strsplit(vax.eff[[i]]$AA[[j]],",")[[1]][1]
  }
}

#filter out seq artefacts
for(i in seq_along(vax.eff)){
  vax.eff[[i]] <- vax.eff[[i]] %>%
    filter(POS != 6696) %>%
    filter(POS != 11074) %>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)
}

#write individual annotation csv files for each user
for(i in 1:length(vax.eff)){
  filenames <- paste0("vaccinated_annotations/annotated_",names(vax.eff),".csv")
  write.csv(vax.eff[[i]],filenames[[i]])
}

#bind annotations for all users into one large table
vax.annotable <- bind_rows(vax.eff)
write.csv(vax.annotable,"vaccinated_annotations.csv")

#reformat variant plots to REF-POS-ALT SNP notation
for(i in 1:length(vax.vartables)){
  for(j in 1:length(vax.vartables[[i]]$POS)){
    vax.vartables[[i]]$POS[[j]] <- paste0(vax.vartables[[i]]$REF[[j]], 
                                            vax.vartables[[i]]$POS[[j]], vax.vartables[[i]]$ALT[[j]])
  }
}
for(i in 1:length(vax.vartables)){
  vax.vartables[[i]] <- vax.vartables[[i]] %>%
    rename("SNP"=POS)%>%
    select(!REF)%>%
    select(!ALT)
}

#transpose data (makes things easier to plot)
for(i in 1:length(vax.vartables)){
  vax.vartables[[i]] <- as.data.frame(t(vax.vartables[[i]]))
  colnames(vax.vartables[[i]]) <- vax.vartables[[i]][1,]
  vax.vartables[[i]] <- vax.vartables[[i]][-1,]
  rownames(vax.vartables[[i]]) <- paste0("freq_",1:dim(vax.vartables[[i]])[[1]])
}

#write csvs for variant frequency tracking tables
for(i in 1:length(vax.vartables)){
  filenames <- paste0("vaccinated_variant_tables/user_",
                      gsub("int_","",names(vax.vartables)),
                      ".csv")
  write.csv(vax.vartables[[i]],filenames[[i]])
}

#count the total number of intersecting (shared) mutations for each user
load("vax_intersecting.RData")
vax_intersecting_cut <- list()
vax_intersect_lengths <- list()
for(i in 1:length(vax_intersecting)){
  vax_intersecting_cut[[i]] <- vax_intersecting[[i]] %>%
    filter(POS!= 6696)%>%
    filter(POS!= 11074)%>%
    filter(POS != 15965)%>%
    filter(POS!= 29051)
  vax_intersect_lengths[[i]] <- length(vax_intersecting_cut[[i]]$POS)
}

#calculate percentage of intersecting (shared) mutants for each user
vax_intersect_percent <- list()
for(i in 1:length(uniquecounts_vax)){
  vax_intersect_percent[[i]] <- vax_intersect_lengths[[i]]/(vax_intersect_lengths[[i]]+uniquecounts_vax[[i]])
}

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
vaccinated_shared <- as.data.frame(matrix(nrow = 9,ncol = 4))
colnames(vaccinated_shared) <- c("Participant ID", "Shared", "Unique", "Total")
vaccinated_shared[,1] <- c("471876","475670","475951","481242","482828",
                                 "484249","486422","487250","487297")
vaccinated_shared[,2] <- as.numeric(unlist(vax_intersect_lengths))
vaccinated_shared[,3] <- as.numeric(unlist(uniquecounts_vax))
vaccinated_shared[,4] <- vaccinated_shared[,2] + vaccinated_shared[,3]
vaccinated_shared$`Participant ID` <- as.character(vaccinated_shared$`Participant ID`)
save(vaccinated_shared,file = "vaccinated_shared.RData")

#plot daily iSNV counts (FIG 2B)
load("snpcounts_dat_vax.RData")
vax_count_gather <- snpcounts_dat_vax%>%
  mutate("day"=1:8)%>%
  gather(key = "Participant ID",value="SNP_count",-day)

twelve.palette <- c("#9D6A90","#8E759F","#7981AA","#618CAF","#4696AE","#2E9EA7",
                    "#27A59B","#37AB8B","#52AF79","#6FB267","#8CB356","#AAB24A",
                    "#C8AF46")

ggplot(data=vax_count_gather,aes(x=`Participant ID`,y=SNP_count,fill = `Participant ID`))+
  geom_dotplot(binaxis = "y",binwidth = .12,stackdir = "center",color=NA)+
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Immune")+
  scale_y_log10(limits = c(1,1000))+
  geom_col(data = vaccinated_shared,aes(x=`Participant ID`,y=Shared), 
           fill=NA,color = "black")+
  geom_col(data=vaccinated_shared,aes(x=`Participant ID`,y=Total),
           fill=NA,color="grey40")+
  scale_fill_manual(values = twelve.palette)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size = 19),
        legend.title = element_text(size = 18),legend.text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle=50, margin = margin(t=3), 
                                   hjust = 1),
        axis.title.x = element_text(margin = margin(t=10)),
        plot.title = element_text(size = 22))
ggsave("figs/vaccinated_SNPs_shared_total_box.png")
save(vax_count_gather,file = "vax_count_gather.Rdata")


