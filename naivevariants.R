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
  snplength <- vector(mode = "list", length(user))
  for (i in seq_along(user)){
    snplength[[i]] <- length(user[[i]]$POS)
  }
  return(snplength)
}
snpcounts <- lapply(sleekuser,snpcount)
snpcounts <- unlist(snpcounts)

user.info <- import("all_saliva_user_info.xlsx")
naive.info <- user.info[which(user.info$status == "unvaccinated"),]
naive.info <- naive.info %>%
  select(user_id,day_of_infection)
naive.snpcounts <- bind_cols(naive.info,snpcounts)
colnames(naive.snpcounts)[3] <- "SNP_count"

save(naive.snpcounts, file = "naive_snpcounts.RData")
write.csv(naive.snpcounts,"naive_snpcounts.csv")

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
  return(unique)
}
uniqueSNPs <- list()
uniquecounts <- list()
for(i in seq_along(allpos)){
  uniqueSNPs[[i]] <- keepunique(allpos[[i]])
  uniquecounts[[i]] <- length(uniqueSNPs[[i]]$POS)
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

naive.intersecting <- lapply(allpos,snp.intersect)
naive.intersecting <- lapply(naive.intersecting,arrange,POS)

save(naive.intersecting,file = "naive_intersecting.RData")

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
naive.daily.shared <- list()
for(i in seq_along(naive.intersecting)){
  naive.daily.shared[[i]] <- shared.lengths(allpos[[i]],naive.intersecting[[i]])
}
names(naive.daily.shared) <- names(naive.intersecting)
save(naive.daily.shared, file = "naive_daily_shared.RData")

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

naive.vartables <- list()
for(i in seq_along(naive.intersecting)){
  naive.vartables[[i]] <- varianttable(naive.intersecting[[i]],everythinguser[[i]])
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
}

#write csvs for variant frequency tracking tables
for(i in seq_along(naive.vartables)){
  filenames <- paste0("naive_variant_tables/user_",names(naive.vartables),".csv")
  write.csv(naive.vartables[[i]],filenames[[i]])
}

#count the total number of intersecting mutations for each user
load("naive_intersecting.RData")
naive.intersect.lengths <- vector(mode = "list", length = length(naive.intersecting))
for(i in seq_along(naive.intersecting)){
  naive.intersect.lengths[[i]] <- length(naive.intersecting[[i]]$POS)
}

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
naive.shared <- as.data.frame(matrix(nrow = 20,ncol = 4))
colnames(naive.shared) <- c("Participant ID", "Shared", "Unique", "Total")
naive.shared[,1] <- names(dirlist)
naive.shared[,2] <- as.numeric(unlist(naive.intersect.lengths))
naive.shared[,3] <- as.numeric(unlist(uniquecounts))
naive.shared[,4] <- naive.shared[,2] + naive.shared[,3]
naive.shared$`Participant ID` <- as.character(naive.shared$`Participant ID`)
naive.shared <- as_tibble(naive.shared)
save(naive.shared, file = "naive_shared.RData")

#plot daily iSNV counts
load("naive_snpcounts.RData")
load("naive_shared.RData")
naive.snpcounts <- as_tibble(naive.snpcounts)
naive.snpcounts$user_id <- as.character(naive.snpcounts$user_id)

twenty.palette <-c("#9D6A90","#94719A","#8978A2","#7C7FA9","#6D86AD","#5D8DAF",
                   "#4C93AF","#3C99AC","#2E9EA7","#26A39F","#2AA796","#36AB8C",
                   "#46AE81","#58B075","#6AB269","#7CB25E","#8FB355","#A2B24D",
                   "#B5B148","#C8AF46")

ggplot(data=naive.snpcounts,aes(x=user_id,y=SNP_count,fill = user_id))+
  geom_dotplot(binaxis = "y", binwidth = .08,stackdir = "center", color = NA)+
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  scale_y_log10()+
  geom_col(data = naive.shared,aes(x=`Participant ID`,y=Shared), 
           fill=NA,color = "black")+
  geom_col(data=naive.shared,aes(x=`Participant ID`,y=Total),
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

ggsave("figs/naive_SNPcounts.png")







