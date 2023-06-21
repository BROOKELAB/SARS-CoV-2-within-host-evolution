library(tidyverse)
library(rio)
library(rlist)
library(gdata)
library(gtools)
library(here)

nasal.dirs <- list.dirs("nasal_output_tables")[-c(1,2)]
nasal.dirs26 <- list.dirs("nasal_output_tables_cut26")[-c(1:3)]
names(nasal.dirs) <- c("user_432686","user_433227","user_438577","user_442978",
                       "user_444332","user_444633","user_450241","user_451152",
                       "user_451709","user_453058","user_459597","user_471467",
                       "user_471588","user_475670")
nasal.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir, pattern = "ivar"))
  full <- mixedsort(sort(full)) #need this so it reads in date order
  files <- lapply(full, read.table,sep="\t",header =T)
  files <- lapply(files, select, POS,REF,ALT,ALT_FREQ,TOTAL_DP)
  return(files)
}
nasal.files <- lapply(nasal.dirs26, nasal.read)
nasal.files.uncut <- lapply(nasal.dirs, nasal.read)

sleekuser.nas <- list()
sleekuser.nas.uncut <- list()
nasal.sleek <- function(file){
  file <- file %>%
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
    distinct()
  return(file)
} #contains depth and frequency cutoffs
for(i in seq_along(nasal.files)){
  sleekuser.nas[[i]] <- lapply(nasal.files[[i]],nasal.sleek)
  sleekuser.nas.uncut[[i]] <- lapply(nasal.files.uncut[[i]],nasal.sleek)
}
names(sleekuser.nas) <- names(nasal.dirs)
names(sleekuser.nas.uncut) <- names(sleekuser.nas)
everythinguser.nas <- nasal.files.uncut

#count number of SNPs for each user
snpcount <- function(user){
  snplength <-list()
  for (i in seq_along(user)){
    snplength[[i]] <- length(user[[i]]$POS)
  }
  return(snplength)
}
snpcounts.nas <- lapply(sleekuser.nas.uncut,snpcount)

snpcounts.nas.table <- data.frame(matrix(ncol = 14,nrow = 9))
rownames(snpcounts.nas.table) <- paste("time point",1:9)
colnames(snpcounts.nas.table) <- names(nasal.dirs)
for(i in seq_along(snpcounts.nas)){
  for(j in seq_along(snpcounts.nas[[i]])){
    snpcounts.nas.table[j,i] <- snpcounts.nas[[i]][[j]]
  }
}

save(snpcounts.nas.table, file = "snpcounts_nas.RData")
write.csv(snpcounts.nas.table,"nasal_snpcounts.csv")

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
allpos.nas.uncut <- lapply(sleekuser.nas.uncut, positions)

keepunique <- function(pos){
  splat <- list.rbind(pos)
  unique  <- distinct(splat)
  return(unique)
}
uniqueSNPs.nas <- list()
uniquecounts.nas <- list()
for(i in seq_along(allpos.nas)){
  uniqueSNPs.nas[[i]] <- keepunique(allpos.nas.uncut[[i]])
  uniquecounts.nas[[i]] <- length(uniqueSNPs.nas[[i]]$POS)
}

#pull out intersecting variants
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

save(nasal.intersecting,file = "nasal_intersecting.RData")

#create variant tables
varianttable <-function(intersectuser,everythinguser){
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
nasal.vartables <- lapply(nasal.vartables,NA_01)

#reformat to REF-POS-ALT SNP notation
for(i in seq_along(nasal.vartables)){
  for(j in seq_along(nasal.vartables[[i]]$POS)){
    nasal.vartables[[i]]$POS[[j]] <- paste0(nasal.vartables[[i]]$REF[[j]], 
                                            nasal.vartables[[i]]$POS[[j]], 
                                            nasal.vartables[[i]]$ALT[[j]])
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


#write csvs
for(i in seq_along(nasal.vartables)){
  filenames <- paste0("nasal_variant_tables/",names(nasal.vartables),".csv")
  write.csv(nasal.vartables[[i]],filenames[[i]])
}





