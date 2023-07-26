library(tidyverse)
library(rio)
library(rlist)
library(Nmisc)
library(gdata)
library(gtools)
library(here)

#load iVar and SnpEff files
dirlist.vax <- list.dirs("vaccinated_output_tables")[-1]

#old
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

sleekuser.vax <- lapply(dirlist.vax,sleek)
everythinguser.vax <- lapply(dirlist.vax,everything)
effuser.vax <- lapply(dirlist.vax,eff)

names(sleekuser.vax) <- c("461913","471876","475670","481242","481672",
                          "482828","484249","486422","487250","487297",
                          "487941")
names(everythinguser.vax) <- names(sleekuser.vax)
names(effuser.vax) <- names(sleekuser.vax)
names(dirlist.vax) <- names(sleekuser.vax)

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

#remove participants with only one time point
vax.intersecting <- vax.intersecting %>%
  discard(~dim(.x)[1] == 0) #users 1,5,11
everythinguser.vax <- everythinguser.vax[-c(1,5,11)]
effuser.vax <- effuser.vax[-c(1,5,11)]

vax.vartables <- list()
for(i in seq_along(vax.intersecting)){
  vax.vartables[[i]] <- varianttable(vax.intersecting[[i]],everythinguser.vax[[i]])
}
names(vax.vartables) <- names(sleekuser.vax)[-c(1,5,11)]

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
vax.vartables <- lapply(vax.vartables,NA_01)

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
vax.vartables <- lapply(vax.vartables, remove.consensus)

#generate list of shared variants
vax.intersecting <- list()
for(i in seq_along(vax.vartables)){
  vax.intersecting[[i]] <- vax.vartables[[i]][,c(1:3)]
}
names(vax.intersecting) <- paste0("user_",names(sleekuser.vax)[-c(1,5,11)])
save(vax.intersecting, file = "vax_intersecting.RData")

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
vax.eff <- list()
for(i in seq_along(vax.vartables)){
  vax.eff[[i]] <- annotable(vax.vartables[[i]],effuser.vax[[i]])
}
names(vax.eff) <- names(sleekuser.vax)[-c(1,5,11)]

#remove participants with no shared iSNVs
vax.eff <- vax.eff %>%
  discard(~dim(.x)[1] == 0)

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
    dplyr::rename("SNP"=POS)%>%
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
vax.daily.shared <- map2(allpos.vax[-c(1,5,11)], vax.intersecting, shared.lengths)
names(vax.daily.shared) <- names(vax.intersecting)
save(vax.daily.shared, file = "vax_daily_shared.RData")

#count the total number of shared iSNVs for each user
vax.intersect.lengths <- vector(mode = "list", length = length(vax.intersecting))
for(i in seq_along(vax.intersecting)){
  vax.intersect.lengths[[i]] <- length(vax.intersecting[[i]]$POS)
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
vax.daily.unique <- map2(sleekuser.vax[-c(1,5,11)], vax.intersecting, unique.lengths)

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
vax.unique.counts <- map2(sleekuser.vax[-c(1,5,11)], vax.intersecting, keep.unique)
vax.unique.counts <- unlist(vax.unique.counts)

#count number of SNPs for each user on each day
vax.daily.counts <- map2(vax.daily.shared, vax.daily.unique, ~.x+.y)
user.info <- import("all_saliva_user_info.xlsx")
vax.info <- user.info[which(user.info$status == "immune"),]
vax.info <- vax.info[-which(vax.info$user_id %in% c("461913","481672","487941")),]
vax.info <- vax.info %>%
  select(user_id,day_of_infection)
vax.snpcounts <- bind_cols(vax.info,
                             unlist(vax.daily.shared),
                             unlist(vax.daily.counts))
colnames(vax.snpcounts)[3] <- "shared_SNPs"
colnames(vax.snpcounts)[4] <- "total_SNPs"

save(vax.snpcounts, file = "vax_snpcounts.RData")
write.csv(vax.snpcounts,"vax_snpcounts.csv")

#compile all iSNV count data (shared iSNVs, unique iSNVs, and total iSNVs)
vax.shared <- as.data.frame(matrix(nrow = 8,ncol = 4))
colnames(vax.shared) <- c("Participant ID", "Shared", "Unique", "Total")
vax.shared[,1] <- names(dirlist.vax)[-c(1,5,11)]
vax.shared[,2] <- as.numeric(unlist(vax.intersect.lengths))
vax.shared[,3] <- as.numeric(vax.unique.counts)
vax.shared[,4] <- vax.shared[,2] + vax.shared[,3]
vax.shared$`Participant ID` <- as.character(vax.shared$`Participant ID`)
vax.shared <- as_tibble(vax.shared)
save(vax.shared, file = "vax_shared.RData")

#plot iSNV information per participant
load("vax_snpcounts.RData")
load("vax_shared.RData")

vax.snpcounts <- as_tibble(vax.snpcounts)
vax.snpcounts$user_id <- as.character(vax.snpcounts$user_id)

eight.palette <- c("#9D6A90","#807DA7","#5690AF","#2C9FA6","#34AA8D",
                   "#62B16E","#94B352","#C8AF46")

ggplot(data=vax.snpcounts,aes(x=user_id,y=total_SNPs))+ 
  geom_dotplot(binaxis = "y", binwidth = .08,stackdir = "center", color = NA, aes(fill = user_id))+
  xlab("Participant ID")+
  ylab("iSNV Count")+
  ggtitle("Unvaccinated")+
  scale_y_log10(limits = c(1,1000))+
  geom_col(data=vax.shared,aes(x=`Participant ID`,y=Total),fill = NA,
           color="grey40")+
  geom_col(data=vax.shared,aes(x=`Participant ID`,y=Shared),fill = NA,
           color="black")+
  scale_fill_manual(values = eight.palette)+
  theme_bw()+
  theme(axis.title = element_text(size=22),axis.text = element_text(size = 19),
        legend.title = element_text(size = 18),legend.text = element_text(size = 15),
        legend.position = "none",
        axis.text.x = element_text(angle=50, margin = margin(t=3), 
                                   hjust = 1),
        axis.title.x = element_text(margin = margin(t=10)),
        plot.title = element_text(size = 22))

ggsave("figs/vax_SNPcounts.png")




