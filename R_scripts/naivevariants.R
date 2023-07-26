library(tidyverse)
library(rio)
library(rlist)
library(Nmisc)
library(gdata)
library(gtools)
library(here)

#load iVar and SnpEff files
dirlist <- list.dirs("naive_output_tables")[-c(1:3)]

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
      filter(ALT_FREQ>=.03)%>%
      filter(TOTAL_DP >=1000)%>%
      select(POS,REF,ALT,ALT_DP,TOTAL_DP,ALT_FREQ)%>%
      distinct()
    artifacts <- c(6696,11074,15965,29051,187,1059,2094,3037,
                   3130,6696,6990,8022,10323,10741,11074,13408,
                   14786,19684,20148,21137,24034,24378,25563,26144,
                   26461,26681,28077,28826,28854,29051,29700,29760)
    if(!identical(which(sleekuser[[i]]$POS %in% artifacts), integer(0))){
      sleekuser[[i]] <- sleekuser[[i]][-which(sleekuser[[i]]$POS %in% artifacts),]
    }
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

#remove samples with low overall coverage
load("naive_depth_files.RData")

rm.coverage <- function(user,depth){
  names(user) <- names(depth)
  if(length(which(flatten(depth) < 1000)) > 0){
    user <- user[-c(which(flatten(depth) < 1000))]
  }
  return(user)
}

sleekuser <- map2(sleekuser,naive.depth.files,rm.coverage)
everythinguser <- map2(everythinguser,naive.depth.files,rm.coverage)
effuser <- map2(effuser,naive.depth.files,rm.coverage)

#remove samples with low Ct
naive.ct <- import("naive_ct.xlsx")
naive.ct <- naive.ct %>%
  group_by(user_id, .add = TRUE) %>%
  group_split()
names(naive.ct) <- names(sleekuser)
sleekuser.uncut <- sleekuser

rm.ct <- function(user,ct){
  discard.list <- c(ct$vdl_barcode[(which(ct$ct > 26))])
  discard.list <- as.character(discard.list[which(discard.list %in% names(user))])
  user.cut <- discard_at(user, all_of(discard.list))
  return(user.cut)
}

sleekuser <- map2(sleekuser,naive.ct,rm.ct)

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
allpos.uncut <- lapply(sleekuser.uncut,positions)

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

#create variant tables to track shared iSNV frequencies over time
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

naive.vartables <- list()
for(i in seq_along(naive.intersecting)){
  naive.vartables[[i]] <- varianttable(naive.intersecting[[i]],everythinguser[[i]])
}
names(naive.vartables) <- names(sleekuser)

#replace NA with 0.01 (for plotting on log scale)
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
naive.vartables <- lapply(naive.vartables,NA_01)

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
naive.vartables <- lapply(naive.vartables, remove.consensus)

#generate list of shared variants
naive.intersecting <- list()
for(i in seq_along(naive.vartables)){
  naive.intersecting[[i]] <- naive.vartables[[i]][,c(1:3)]
}
names(naive.intersecting) <- paste0("user_",names(sleekuser))
save(naive.intersecting, file = "naive_intersecting.RData")

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
        if(dim(annotable)[1] != 0 &&
           (effuser[[i]]$POS[[j]] == annotable$POS[[x]]) &&
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
naive.eff <- naive.eff %>%
  discard(~ dim(.x)[1] == 0)
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
    dplyr::rename("SNP"=POS)%>%
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
naive.daily.shared <- map2(allpos.uncut, naive.intersecting, shared.lengths)
names(naive.daily.shared) <- names(naive.intersecting)
save(naive.daily.shared, file = "naive_daily_shared.RData")

#count the total number of shared iSNVs for each user
naive.intersect.lengths <- vector(mode = "list", length = length(naive.intersecting))
for(i in seq_along(naive.intersecting)){
  naive.intersect.lengths[[i]] <- length(naive.intersecting[[i]]$POS)
}

#count iSNVs that do NOT appear more than once within users
unique.lengths <- function(user, intersecting){
  unique <- list()
  lengths <- list()
  for(i in seq_along(user)){
    user[[i]] <- user[[i]] %>%
      mutate("SNP" = paste0(REF,POS,ALT)) %>%
      filter(ALT_FREQ <= 0.97)
    intersecting <-  mutate(intersecting, "SNP" = paste0(REF,POS,ALT))
    if(dim(intersecting)[[1]] > 0){
      unique[[i]] <- user[[i]][-which(user[[i]]$SNP %in% intersecting$SNP),]
    } else {
      unique[[i]] <- user[[i]]
    }
    lengths[[i]] <- length(unique[[i]]$POS)
  }
  lengths <- unlist(lengths)
  return(lengths)
}
naive.daily.unique <- map2(sleekuser.uncut, naive.intersecting, unique.lengths)

keep.unique <- function(user, intersecting){
  unique <- list()
  for(i in seq_along(user)){
    user[[i]] <- user[[i]] %>%
      mutate("SNP" = paste0(REF,POS,ALT)) %>%
      filter(ALT_FREQ <= 0.97)
    intersecting <-  mutate(intersecting, "SNP" = paste0(REF,POS,ALT))
    if(dim(intersecting)[[1]] > 0){
      unique[[i]] <- user[[i]][-which(user[[i]]$SNP %in% intersecting$SNP),]
    } else {
      unique[[i]] <- user[[i]]
    }
  }
  unique <- bind_rows(unique)
  count <- length(unique$POS)
  return(count)
}

naive.unique.counts <- map2(sleekuser.uncut, naive.intersecting, keep.unique)
naive.unique.counts <- unlist(naive.unique.counts)

#count number of SNPs for each user on each day
naive.daily.counts <- map2(naive.daily.shared, naive.daily.unique, ~.x+.y)
user.info <- import("all_saliva_user_info.xlsx")
naive.info <- user.info[which(user.info$status == "unvaccinated"),]
naive.info <- naive.info %>%
  select(user_id,day_of_infection)
naive.snpcounts <- bind_cols(naive.info,
                             unlist(naive.daily.shared),
                             unlist(naive.daily.counts))
colnames(naive.snpcounts)[3] <- "shared_SNPs"
colnames(naive.snpcounts)[4] <- "total_SNPs"

save(naive.snpcounts, file = "naive_snpcounts.RData")
write.csv(naive.snpcounts,"naive_snpcounts.csv")

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
naive.shared <- as.data.frame(matrix(nrow = 20,ncol = 4))
colnames(naive.shared) <- c("Participant ID", "Shared", "Unique", "Total")
naive.shared[,1] <- names(dirlist)
naive.shared[,2] <- as.numeric(unlist(naive.intersect.lengths))
naive.shared[,3] <- as.numeric(naive.unique.counts)
naive.shared[,4] <- naive.shared[,2] + naive.shared[,3]
naive.shared$`Participant ID` <- as.character(naive.shared$`Participant ID`)
naive.shared <- as_tibble(naive.shared)
save(naive.shared, file = "naive_shared.RData")


#plot iSNV information per participant
load("naive_snpcounts.RData")
load("naive_shared.RData")

naive.snpcounts <- as_tibble(naive.snpcounts)
naive.snpcounts$user_id <- as.character(naive.snpcounts$user_id)

twenty.palette <-c("#9D6A90","#94719A","#8978A2","#7C7FA9","#6D86AD","#5D8DAF",
                   "#4C93AF","#3C99AC","#2E9EA7","#26A39F","#2AA796","#36AB8C",
                   "#46AE81","#58B075","#6AB269","#7CB25E","#8FB355","#A2B24D",
                   "#B5B148","#C8AF46")

ggplot(data=naive.snpcounts,aes(x=user_id,y=total_SNPs))+ 
  geom_dotplot(binaxis = "y", binwidth = .08,stackdir = "center", color = NA, aes(fill = user_id))+
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  scale_y_log10()+
  geom_col(data=naive.shared,aes(x=`Participant ID`,y=Total),fill = NA,
           color="grey40")+
  geom_col(data=naive.shared,aes(x=`Participant ID`,y=Shared),fill = NA,
           color="black")+
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





