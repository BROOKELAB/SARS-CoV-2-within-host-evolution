library(tidyverse)
library(rio)
library(rlist)
library(gdata)
library(gtools)
library(here)

dirlist.vax <- list.dirs("vaccinated_output_tables")[-1]

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
      select(POS,REF,ALT,GENE,EFF,AA) %>%
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
      filter(POS != 29700)
  }
  return(effuser)
} #gene annotation files

sleekuser.vax <- lapply(dirlist.vax,sleek)
everythinguser.vax <- lapply(dirlist.vax,everything)
effuser.vax <- lapply(dirlist.vax,eff)

names(sleekuser.vax) <- c("461913","471876","475670","475951","481242","481672",
                      "482828","484249","486422","487250","487297",
                      "487941")
names(everythinguser.vax) <- names(sleekuser.vax)
names(effuser.vax) <- names(sleekuser.vax)
names(dirlist.vax) <- names(sleekuser.vax)

#count number of SNPs for each user on each day
snpcount <- function(user){
  snplength <- vector(mode = "list", length(user))
  for (i in seq_along(user)){
    snplength[[i]] <- length(user[[i]]$POS)
  }
  return(snplength)
}
snpcounts.vax <- lapply(sleekuser.vax,snpcount)
snpcounts.vax <- unlist(snpcounts.vax)

user.info <- import("all_saliva_user_info.xlsx")
vax.info <- user.info[which(user.info$status == "immune"),]
vax.info <- vax.info %>%
  select(user_id,day_of_infection)
snpcounts.vax <- bind_cols(vax.info,snpcounts.vax)
colnames(snpcounts.vax)[3] <- "SNP_count"

save(snpcounts.vax, file = "vax_snpcounts.RData")
write.csv(snpcounts.vax,"vax_snpcounts.csv")

#pull out SNP positions
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
allpos.vax <- lapply(sleekuser.vax,positions)

#pull out variants that do NOT appear more than once within users
keepunique <- function(pos){
  splat <- list.rbind(pos)
  unique  <- distinct(splat)
  return(unique)
}
uniqueSNPs.vax <- list()
uniquecounts.vax <- list()
for(i in seq_along(allpos.vax)){
  uniqueSNPs.vax[[i]] <- keepunique(allpos.vax[[i]])
  uniquecounts.vax[[i]] <- length(uniqueSNPs.vax[[i]]$POS)
}

#pull out intersecting iSNVs
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

vax.intersecting <- lapply(allpos.vax,snp.intersect)
for(i in seq_along(vax.intersecting)){
  if(length(vax.intersecting[[i]]) > 0){
    vax.intersecting[[i]] <- vax.intersecting[[i]] %>%
      arrange(POS)
  }
}

save(vax.intersecting,file = "vax_intersecting.RData")

#count shared iSNVs present on each day
shared.lengths <- function(all,intersecting){
  subset <- list()
  subset.lengths <- list()
  for(i in seq_along(all)){
    subset[[i]] <- which(all[[i]]$POS %in% intersecting$POS)
    subset.lengths[[i]] <- length(subset[[i]])
  }
  subset.lengths <- unlist(subset.lengths)
  return(subset.lengths)
}
vax.daily.shared <- list()
for(i in seq_along(vax.intersecting)){
  vax.daily.shared[[i]] <- shared.lengths(allpos.vax[[i]],vax.intersecting[[i]])
}
names(vax.daily.shared) <- names(vax.intersecting)
save(vax.daily.shared, file = "vax_daily_shared.RData")

#create variant tables to track shared iSNV frequencies over time
varianttable <- function(intersectuser,everythinguser){
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

#remove users with no shared iSNVs
vax.intersecting <- vax.intersecting[-c(1,4,6,8,9,11,12)]
everythinguser.vax <- everythinguser.vax[-c(1,4,6,8,9,11,12)]
effuser.vax <- effuser.vax[-c(1,4,6,8,9,11,12)]

vax.vartables <- list()
for(i in seq_along(vax.intersecting)){
  vax.vartables[[i]] <- varianttable(vax.intersecting[[i]],everythinguser.vax[[i]])
}
names(vax.vartables) <- names(everythinguser.vax)

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
vax.vartables <- lapply(vax.vartables,NA_01)

#reformat SNP notation with no +/-'s (eg. REF = C, ALT= +T --> REF = C, ALT = CT)
for(i in seq_along(vax.vartables)){
  for(j in seq_along(vax.vartables[[i]]$ALT)){
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
vax.eff <- list()
for(i in seq_along(vax.vartables)){
  vax.eff[[i]] <- annotable(vax.vartables[[i]],effuser.vax[[i]])
}
names(vax.eff) <- names(effuser.vax)

#annotate UTR as 5' or 3' UTR
for(i in seq_along(vax.eff)){
  for(j in seq_along(vax.eff[[i]]$POS)){
    if(vax.eff[[i]]$POS[[j]] <= 265){
      vax.eff[[i]]$GENE[[j]] <- "5'UTR"
    }
    if(vax.eff[[i]]$POS[[j]] >= 29675){
      vax.eff[[i]]$GENE[[j]] <- "3'UTR"
    }
    if(vax.eff[[i]]$POS[[j]] >= 28260 && vax.eff[[i]]$POS[[j]] < 28274){
      vax.eff[[i]]$GENE[[j]] <- "UTR"
    }
  }
}

#add columns for ID and reformat to REF-POS-ALT SNP notation
for(i in seq_along(vax.eff)){
  vax.eff[[i]] <- vax.eff[[i]] %>%
    mutate("ID" = names(vax.eff)[[i]]) %>%
    relocate(ID,.before = POS) %>%
    mutate("SNP"=paste0(REF,POS,ALT))%>%
    relocate(SNP,.after = ID)%>%
    select(ID,SNP,POS,REF,ALT,GENE,EFF,AA)
}

#reformat EFF column to read synonymous/nonsynonymous/utr
for(i in seq_along(vax.eff)){
  for(j in seq_along(vax.eff[[i]]$GENE)){
    if((vax.eff[[i]]$GENE[[j]] == "5'UTR") |
       (vax.eff[[i]]$GENE[[j]] == "3'UTR") |
       (vax.eff[[i]]$GENE[[j]] == "UTR")) {
      vax.eff[[i]]$EFF[[j]] <- "UTR"
    }
    if(vax.eff[[i]]$EFF[[j]] == "synonymous_variant"){
      vax.eff[[i]]$EFF[[j]] <- "synonymous"
    }
  }
}
for(i in seq_along(vax.eff)){
  for(j in seq_along(vax.eff[[i]]$EFF)){
    if((vax.eff[[i]]$EFF[[j]] != "UTR") &&
       (vax.eff[[i]]$EFF[[j]] != "synonymous")){
      vax.eff[[i]]$EFF[[j]] <- "nonsynonymous"
    }
  }
}

#get rid of extraneous punctuation in AA column
for(i in seq_along(vax.eff)){
  for(j in seq_along(vax.eff[[i]]$AA)){
    vax.eff[[i]]$AA[[j]] <- strsplit(vax.eff[[i]]$AA[[j]],",")[[1]][1]
  }
}

#bind annotations for all users into one large table
vax.annotable <- bind_rows(vax.eff)
write.csv(vax.annotable,"vax_annotations.csv")

#reformat variant plots to REF-POS-ALT SNP notation
for(i in seq_along(vax.vartables)){
  for(j in seq_along(vax.vartables[[i]]$POS)){
    vax.vartables[[i]]$POS[[j]] <- paste0(vax.vartables[[i]]$REF[[j]], 
                                            vax.vartables[[i]]$POS[[j]], vax.vartables[[i]]$ALT[[j]])
  }
}

for(i in seq_along(vax.vartables)){
  vax.vartables[[i]] <- vax.vartables[[i]] %>%
    rename("SNP"=POS)%>%
    select(!REF)%>%
    select(!ALT)
}

#transpose data (makes things easier to plot)
for(i in seq_along(vax.vartables)){
  vax.vartables[[i]] <- as.data.frame(t(vax.vartables[[i]]))
  colnames(vax.vartables[[i]]) <- vax.vartables[[i]][1,]
}

#write csvs for variant frequency tracking tables
for(i in seq_along(vax.vartables)){
  filenames <- paste0("vaccinated_variant_tables/user_",names(vax.vartables),".csv")
  write.csv(vax.vartables[[i]],filenames[[i]])
}

#count the total number of intersecting mutations for each user
load("vax_intersecting.RData")
vax.intersect.lengths <- vector(mode = "list", length = length(vax.intersecting))
for(i in seq_along(vax.intersecting)){
  vax.intersect.lengths[[i]] <- length(vax.intersecting[[i]]$POS)
}

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
vax.shared <- as.data.frame(matrix(nrow = 12,ncol = 4))
colnames(vax.shared) <- c("Participant ID", "Shared", "Unique", "Total")
vax.shared[,1] <- names(dirlist.vax)
vax.shared[,2] <- as.numeric(unlist(vax.intersect.lengths))
vax.shared[,3] <- as.numeric(unlist(uniquecounts.vax))
vax.shared[,4] <- vax.shared[,2] + vax.shared[,3]
vax.shared$`Participant ID` <- as.character(vax.shared$`Participant ID`)
vax.shared <- as_tibble(vax.shared)
save(vax.shared, file = "vax_shared.RData")

#plot daily iSNV counts
load("vax_snpcounts.RData")
load("vax_shared.RData")
snpcounts.vax <- as_tibble(snpcounts.vax)
snpcounts.vax$user_id <- as.character(snpcounts.vax$user_id)

twelve.palette <- c("#9D6A90","#8E759F","#7981AA","#618CAF","#4696AE","#2E9EA7","#27A59B",
                    "#37AB8B","#52AF79","#6FB267","#8CB356","#AAB24A","#C8AF46")

vax.shared[c(1,6,12),] <- NA #only 1 datapoint

ggplot(data=snpcounts.vax,aes(x=user_id,y=SNP_count,fill = user_id))+
  geom_dotplot(binaxis = "y", binwidth = .08,stackdir = "center", color = NA)+
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Immune")+
  scale_y_log10()+
  geom_col(data = vax.shared,aes(x=`Participant ID`,y=Shared), 
           fill=NA,color = "black")+
  geom_col(data=vax.shared,aes(x=`Participant ID`,y=Total),
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

ggsave("figs/vax_SNPcounts.png")
