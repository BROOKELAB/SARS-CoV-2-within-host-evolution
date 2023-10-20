library(tidyverse)
library(rio)
library(rlist)
library(Nmisc)
library(gdata)
library(gtools)
library(here)

#load iVar and SnpEff files
nasal.dirs <- list.dirs("nasal_output_tables")[-c(1,2)]
names(nasal.dirs) <- c("user_432686","user_433227","user_438577","user_442978",
                       "user_444332","user_444633","user_450241","user_451152",
                       "user_451709","user_453058","user_459597","user_471467",
                       "user_471588")

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
    artifacts <- c(6696,11074,15965,29051,78,187,635,1059,2094,3037,
                   3130,6696,6990,8022,10323,10741,11074,11567,13408,
                   14786,19684,20148,21137,24034,24378,25563,26144,
                   26461,26681,27964,28077,28253,28262,28281,28472,
                   28826,28854,29051,29700,29760)
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

sleekuser.nas <- lapply(nasal.dirs,sleek)
everythinguser.nas <- lapply(nasal.dirs,everything)

#remove samples with low overall coverage
load("nasal_depth_files.RData")

rm.coverage <- function(user,depth){
  names(user) <- names(depth)
  if(length(which(flatten(depth) < 1000)) > 0){
    user <- user[-c(which(flatten(depth) < 1000))]
  }
  return(user)
}

sleekuser.nas <- map2(sleekuser.nas,nasal.depth.files,rm.coverage)
everythinguser.nas <- map2(everythinguser.nas,nasal.depth.files,rm.coverage)

#remove samples with low Ct
nasal.ct <- import("all_nasal_user_info.xlsx")
nasal.ct <- nasal.ct %>%
  filter(user_ID %in% gsub("user_","",names(sleekuser.nas)))%>%
  arrange(user_ID) %>%
  mutate("sample" = NA) %>%
  select(user_ID,sample,sample_date,cn)

dates <- nasal.ct$sample_date
dates <- format(dates,"%m-%d-%Y")
dates <- as.character(dates)
m <- str_split(dates, "-") %>% map_chr(`[`, 1)
m <- gsub("0","",m)
d <- str_split(dates, "-") %>% map_chr(`[`, 2)
d1 <- gsub("0","",substr(d,1,1))
d2 <- substr(d,2,2)
y <- str_split(dates, "-") %>% map_chr(`[`, 3)
y2 <- substr(y,3,4)
re.dates <- paste0(m,"_",d1,d2,"_",y2)
nasal.ct$sample <- paste0(nasal.ct$user_ID,"_",re.dates)

nasal.ct <- nasal.ct %>%
  group_by(user_ID, .add = TRUE) %>%
  group_split()
names(nasal.ct) <- names(sleekuser.nas)
sleekuser.nas.uncut <- sleekuser.nas

rm.vec.ct <- function(user,ct){
  ct <- filter(ct, sample %in% names(user))
  vec <- which(ct$cn >= 26)
  return(vec)
}
nasal.vec.ct <- map2(sleekuser.nas,nasal.ct, rm.vec.ct)

for(i in seq_along(nasal.vec.ct)){
  if(length(nasal.vec.ct[[i]]) != 0){
    sleekuser.nas[[i]] <- sleekuser.nas[[i]][-c(nasal.vec.ct[[i]])]
  } else {
    sleekuser.nas[[i]] <- sleekuser.nas[[i]]
  }
}

#pull out SNP positions
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
allpos.nas <- lapply(sleekuser.nas,positions)
allpos.nas.uncut <- lapply(sleekuser.nas.uncut,positions)

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
nasal.intersecting <- lapply(allpos.nas,snp.intersect)
nasal.intersecting <- lapply(nasal.intersecting,arrange,POS)

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

nasal.vartables <- list()
for(i in seq_along(nasal.intersecting)){
  nasal.vartables[[i]] <- varianttable(nasal.intersecting[[i]],everythinguser.nas[[i]])
}
names(nasal.vartables) <- names(sleekuser.nas)

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
nasal.vartables <- lapply(nasal.vartables,NA_01)

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
nasal.vartables <- lapply(nasal.vartables, remove.consensus)

#generate list of shared variants
nasal.intersecting <- list()
for(i in seq_along(nasal.vartables)){
  nasal.intersecting[[i]] <- nasal.vartables[[i]][,c(1:3)]
}
names(nasal.intersecting) <- paste0("user_",names(sleekuser.nas))
save(nasal.intersecting, file = "nasal_intersecting.RData")

#reformat SNP notation with no +/-'s (eg. REF = C, ALT= +T --> REF = C, ALT = CT)
for(i in seq_along(nasal.vartables)){
  for(j in seq_along(nasal.vartables[[i]]$ALT)){
    if(substr(nasal.vartables[[i]]$ALT[[j]],0,1)=="+"){
      nasal.vartables[[i]]$ALT[[j]] <- paste0(nasal.vartables[[i]]$REF[[j]], 
                                              substr(nasal.vartables[[i]]$ALT[[j]],2,nchar(nasal.vartables[[i]]$ALT[[j]])))
    }
    if(substr(nasal.vartables[[i]]$ALT[[j]],0,1)=="-"){
      nasal.vartables[[i]]$REF[[j]] <- paste0(nasal.vartables[[i]]$REF[[j]],
                                              substr(nasal.vartables[[i]]$ALT[[j]],2,nchar(nasal.vartables[[i]]$ALT[[j]])))
      nasal.vartables[[i]]$ALT[[j]] <- substr(nasal.vartables[[i]]$REF[[j]],0,1)
    }
  }
}

#reformat variant plots to REF-POS-ALT SNP notation
for(i in seq_along(nasal.vartables)){
  for(j in seq_along(nasal.vartables[[i]]$POS)){
    nasal.vartables[[i]]$POS[[j]] <- paste0(nasal.vartables[[i]]$REF[[j]], 
                                            nasal.vartables[[i]]$POS[[j]], nasal.vartables[[i]]$ALT[[j]])
  }
}

for(i in seq_along(nasal.vartables)){
  nasal.vartables[[i]] <- nasal.vartables[[i]] %>%
    dplyr::rename("SNP"=POS)%>%
    select(!REF)%>%
    select(!ALT)
}

#transpose data (makes things easier to plot)
for(i in seq_along(nasal.vartables)){
  nasal.vartables[[i]] <- as.data.frame(t(nasal.vartables[[i]]))
  colnames(nasal.vartables[[i]]) <- nasal.vartables[[i]][1,]
}

#write csvs for variant frequency tracking tables
for(i in seq_along(nasal.vartables)){
  filenames <- paste0("nasal_variant_tables/",names(nasal.vartables),".csv")
  write.csv(nasal.vartables[[i]],filenames[[i]])
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
nasal.daily.shared <- map2(allpos.nas.uncut, nasal.intersecting, shared.lengths)
save(nasal.daily.shared, file = "nasal_daily_shared.RData")

#count the total number of shared iSNVs for each user
nasal.intersect.lengths <- vector(mode = "list", length = length(nasal.intersecting))
for(i in seq_along(nasal.intersecting)){
  nasal.intersect.lengths[[i]] <- length(nasal.intersecting[[i]]$POS)
}
names(nasal.intersect.lengths) <- names(nasal.daily.shared)

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
nasal.daily.unique <- map2(sleekuser.nas.uncut, nasal.intersecting, unique.lengths)

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

nasal.unique.counts <- map2(sleekuser.nas.uncut, nasal.intersecting, keep.unique)
nasal.unique.counts <- unlist(nasal.unique.counts)

#count number of SNPs for each user on each day
nasal.daily.counts <- map2(nasal.daily.shared, nasal.daily.unique, ~.x+.y)
nasal.info <- import("all_nasal_user_info.xlsx")
nasal.info <- nasal.info %>%
  select(user_ID,day_of_infection)
nasal.snpcounts <- bind_cols(nasal.info,
                             unlist(nasal.daily.shared),
                             unlist(nasal.daily.counts))

colnames(nasal.snpcounts)[3] <- "shared_SNPs"
colnames(nasal.snpcounts)[4] <- "total_SNPs"

save(nasal.snpcounts, file = "nasal_snpcounts.RData")
write.csv(nasal.snpcounts,"nasal_snpcounts.csv")

