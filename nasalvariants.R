library(tidyverse)
library(rio)
library(rlist)
library(gdata)
library(gtools)
library(here)

nasal.dirs <- list.dirs("nasal_output_tables")[-1]
names(nasal.dirs) <- c("user_432686","user_433227","user_438577","user_442978",
                       "user_444332","user_444633","user_450241","user_451152",
                       "user_451709","user_453058","user_459597","user_471467",
                       "user_471588","user_475670")

nasal.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir))
  full <- mixedsort(sort(full))
  files <- lapply(full, read.table,sep="\t")
  return(files)
}
nasal.paths <- function(dir){
  full <- list.files(dir)
  name <- gsub(".hard-filtered.vcf","",full)
  name <- mixedsort(sort(name))
  return(name)
}
nasal.read2 <- function(file){
  file <- file %>%
    select(V2,V4,V5,V8,V10)%>%
    mutate("ID"=NA)%>%
    rename("POS"=V2)%>%
    rename("REF"=V4)%>%
    rename("ALT"=V5)%>%
    rename("INFO1"=V8)%>%
    rename("INFO2"=V10)%>%
    mutate("ALT_FREQ"=NA)%>%
    mutate("DP"=NA)
  for(i in seq_along(file$INFO1)){
    file$DP[[i]] <- strsplit(file$INFO1[[i]],";")[[1]][1]
    file$DP[[i]] <- as.numeric(strsplit(file$DP[[i]],"=")[[1]][2])
    file$ALT_FREQ[[i]] <- strsplit(file$INFO2[[i]],":")[[1]][4]
  }
  file$DP <- as.numeric(file$DP)
  file <- file %>%
    select("POS","REF","ALT","ALT_FREQ","DP")%>%
    distinct()
  return(file)
}

nasal.names <- lapply(nasal.dirs,nasal.paths)
nasal.files <- lapply(nasal.dirs,nasal.read)
for(i in seq_along(nasal.files)){
  nasal.files[[i]] <- lapply(nasal.files[[i]],nasal.read2)
}

nasal.files$user_451709 <- nasal.files$user_451709[-c(2,5,7)] #removal due to potential contamination

sleekuser.nas <- list()
nasal.sleek <- function(file){
  file <- file %>%
    filter(ALT_FREQ >= .03) %>% 
    filter(ALT_FREQ <= .97)%>%
    filter(DP >= 200) %>%
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
    filter(POS != 25314) %>%
    filter(POS != 25317) %>%
    filter(POS != 25324) %>%
    filter(POS != 25563) %>%
    filter(POS != 26144) %>%
    filter(POS != 26461) %>%
    filter(POS != 26681) %>%
    filter(POS != 28077) %>%
    filter(POS != 28826) %>%
    filter(POS != 28854) %>%
    filter(POS != 29051) %>%
    filter(POS != 29700)
  return(file)
} #contains depth and frequency cutoffs
for(i in seq_along(nasal.files)){
  sleekuser.nas[[i]] <- lapply(nasal.files[[i]],nasal.sleek)
}
names(sleekuser.nas) <- names(nasal.dirs)
everythinguser.nas <- nasal.files

#count number of SNPs for each user
snpcount <- function(user){
  snplength <-list()
  for (i in seq_along(user)){
    user[[i]] <- user[[i]] %>%
      filter(POS != 6696)%>%
      filter(POS!= 11074)%>%
      filter(POS!= 29051)
    snplength[[i]] <- length(user[[i]]$POS)
  }
  return(snplength)
}
snpcounts.nas <- lapply(sleekuser.nas,snpcount)

snpcounts.nas <- data.frame(matrix(ncol = 14,nrow = 10))
rownames(snpcounts.nas) <- paste("time point",1:10)
colnames(snpcounts.nas) <- names(nasal.dirs)
for(i in seq_along(snpcounts.nas)){
  for(j in seq_along(snpcounts.nas[[i]])){
    snpcounts.nas[j,i] <- snpcounts.nas[[i]][[j]]
  }
}

write.csv(snpcounts.nas,"nasal_snpcounts.csv")

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

keepunique <- function(pos){
  splat <- list.rbind(pos)
  unique  <- distinct(splat)
  unique <- unique %>%
    filter(POS != 6696)%>%
    filter(POS != 11074)%>%
    filter(POS != 29051)
  return(unique)
}
uniqueSNPs.nas <- list()
uniquecounts.nas <- list()
for(i in seq_along(allpos.nas)){
  uniqueSNPs.nas[[i]] <- keepunique(allpos.nas[[i]])
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


#create variant plots
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
        if((!is.na(joiner[i,j])) && (joiner[i,j]<200)){
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
nasal.vartables <- list()
for(i in seq_along(nasal.intersecting)){
  nasal.vartables[[i]] <- varianttable(nasal.intersecting[[i]],everythinguser.nas[[i]])
}
names(nasal.vartables) <- names(sleekuser.nas)

#replace NA with 0.01
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
    rename("SNP"=POS)%>%
    select(!REF)%>%
    select(!ALT)
}

#transpose data (makes things easier to plot)
for(i in seq_along(nasal.vartables)){
  nasal.vartables[[i]] <- as.data.frame(t(nasal.vartables[[i]]))
  colnames(nasal.vartables[[i]]) <- nasal.vartables[[i]][1,]
  rownames(nasal.vartables[[i]]) <- paste0("freq_",1:dim(nasal.vartables[[i]])[[1]])
}

#write csvs
for(i in seq_along(nasal.vartables)){
  filenames <- paste0("nasal_variant_tables/",names(nasal.vartables),".csv")
  write.csv(nasal.vartables[[i]],filenames[[i]])
}





