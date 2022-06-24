library(tidyverse)
library(rio)
library(rlist)
library(gdata)
library(compiler)
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

sleekuser_nas <- list()
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
  sleekuser_nas[[i]] <- lapply(nasal.files[[i]],nasal.sleek)
}
names(sleekuser_nas) <- names(nasal.dirs)
everythinguser_nas <- nasal.files

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
snpcounts_nas <- lapply(sleekuser_nas,snpcount)

snpcounts_dat_nas <- data.frame(matrix(ncol = 14,nrow = 10))
rownames(snpcounts_dat_nas) <- paste("time point",1:10)
colnames(snpcounts_dat_nas) <- names(nasal.dirs)
for(i in seq_along(snpcounts_nas)){
  for(j in seq_along(snpcounts_nas[[i]])){
    snpcounts_dat_nas[j,i] <- snpcounts_nas[[i]][[j]]
  }
}

#save(snpcounts_dat_nas, file = "snpcounts_dat_nas.RData")
#write.csv(snpcounts_dat_nas,"nasal_snpcounts.csv")

#pull out SNP positions
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
allpos_nas <- lapply(sleekuser_nas,positions)

keepunique <- function(pos){
  splat <- list.rbind(pos)
  unique  <- distinct(splat)
  unique <- unique %>%
    filter(POS != 6696)%>%
    filter(POS != 11074)%>%
    filter(POS != 29051)
  return(unique)
}
uniqueSNPs_nas <- list()
uniquecounts_nas <- list()
for(i in seq_along(allpos_nas)){
  uniqueSNPs_nas[[i]] <- keepunique(allpos_nas[[i]])
  uniquecounts_nas[[i]] <- length(uniqueSNPs_nas[[i]]$POS)
}

#pull out intersecting variants
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

var_432686 <- list("time1"=var1(allpos_nas[[1]]),"time2"=var2(allpos_nas[[1]]),
                   "time3"=var3(allpos_nas[[1]]),"time4"=var4(allpos_nas[[1]]))
var_433227 <- list("time1"=var1(allpos_nas[[2]]),"time2"=var2(allpos_nas[[2]]),
                   "time3"=var3(allpos_nas[[2]]),"time4"=var4(allpos_nas[[2]]),
                   "time5"=var5(allpos_nas[[2]]),"time6"=var6(allpos_nas[[2]]))
var_438577 <- list("time1"=var1(allpos_nas[[3]]),"time2"=var2(allpos_nas[[3]]),
                   "time3"=var3(allpos_nas[[3]]),"time4"=var4(allpos_nas[[3]]),
                   "time5"=var5(allpos_nas[[3]]),"time6"=var6(allpos_nas[[3]]))
var_442978 <- list("time1"=var1(allpos_nas[[4]]),"time2"=var2(allpos_nas[[4]]),
                   "time3"=var3(allpos_nas[[4]]),"time4"=var4(allpos_nas[[4]]),
                   "time5"=var5(allpos_nas[[4]]),"time6"=var6(allpos_nas[[4]]))
var_444332 <- list("time1"=var1(allpos_nas[[5]]),"time2"=var2(allpos_nas[[5]]),
                   "time3"=var3(allpos_nas[[5]]),"time4"=var4(allpos_nas[[5]]),
                   "time5"=var5(allpos_nas[[5]]),"time6"=var6(allpos_nas[[5]]),
                   "time7"=var7(allpos_nas[[5]]),"time8"=var8(allpos_nas[[5]]),
                   "time9"=var9(allpos_nas[[5]]),"time10"=var10(allpos_nas[[5]]))
var_444633 <- list("time1"=var1(allpos_nas[[6]]),"time2"=var2(allpos_nas[[6]]),
                   "time3"=var3(allpos_nas[[6]]),"time4"=var4(allpos_nas[[6]]),
                   "time5"=var5(allpos_nas[[6]]))
var_450241 <- list("time1"=var1(allpos_nas[[7]]),"time2"=var2(allpos_nas[[7]]),
                   "time3"=var3(allpos_nas[[7]]))
var_451152 <- list("time1"=var1(allpos_nas[[8]]),"time2"=var2(allpos_nas[[8]]),
                   "time3"=var3(allpos_nas[[8]]),"time4"=var4(allpos_nas[[8]]),
                   "time5"=var5(allpos_nas[[8]]))
var_451709 <- list("time1"=var1(allpos_nas[[9]]),"time2"=var2(allpos_nas[[9]]),
                   "time3"=var3(allpos_nas[[9]]),"time4"=var4(allpos_nas[[9]]),
                   "time5"=var5(allpos_nas[[9]]),"time6"=var6(allpos_nas[[9]]),
                   "time7"=var7(allpos_nas[[9]]),"time8"=var8(allpos_nas[[9]]))
var_453058 <- list("time1"=var1(allpos_nas[[10]]),"time2"=var2(allpos_nas[[10]]),
                   "time3"=var3(allpos_nas[[10]]),"time4"=var4(allpos_nas[[10]]),
                   "time5"=var5(allpos_nas[[10]]))
var_459597 <- list("time1"=var1(allpos_nas[[11]]),"time2"=var2(allpos_nas[[11]]),
                   "time3"=var3(allpos_nas[[11]]),"time4"=var4(allpos_nas[[11]]),
                   "time5"=var5(allpos_nas[[11]]))
var_471467 <- list("time1"=var1(allpos_nas[[12]]),"time2"=var2(allpos_nas[[12]]),
                   "time3"=var3(allpos_nas[[12]]),"time4"=var4(allpos_nas[[12]]),
                   "time5"=var5(allpos_nas[[12]]))
var_471588 <- list("time1"=var1(allpos_nas[[13]]),"time2"=var2(allpos_nas[[13]]),
                   "time3"=var3(allpos_nas[[13]]),"time4"=var4(allpos_nas[[13]]))
var_475670 <- list("time1"=var1(allpos_nas[[14]]),"time2"=var2(allpos_nas[[14]]),
                   "time3"=var3(allpos_nas[[14]]),"time4"=var4(allpos_nas[[14]]),
                   "time5"=var5(allpos_nas[[14]]))

var_list_nas <- list("var_432686"=var_432686,"var_433227"=var_433227,
                 "var_438577"=var_438577,"var_442978"=var_442978,
                 "var_444332"=var_444332,"var_444633"=var_444633,
                 "var_450241"=var_450241,"var_451152"=var_451152,
                 "var_451709"=var_451709,"var_453058"=var_453058,
                 "var_459597"=var_459597,"var_471467"=var_471467,
                 "var_471588"=var_471588,"var_475670"=var_475670)

#remove self-comparison
selfevict <- function(var){
  for(i in 1:length(var)){
    var[[i]][i] <- NA
  }
  return(var)
}
var_list_nas <- lapply(var_list_nas,selfevict)

#write file for each user to a .csv
var_filenames <- paste0("nasal_intersect/",names(var_list_nas),".csv")
for(i in seq_along(var_list_nas)){
  capture.output(var_list_nas[[i]],file=var_filenames[[i]])
}

#pull out a list of mutations that appear > once for each user
keepcommon <- function(var){
  flat <- flatten(var)
  splat <- na.omit(rlist::list.rbind(flat))
  common <- distinct(splat)
  return(common)
}
nasal_intersecting <- lapply(var_list_nas,keepcommon)
names(nasal_intersecting) <- gsub("var","int",names(var_list_nas))
nasal_intersecting <- lapply(nasal_intersecting,arrange,POS)

#save(nasal_intersecting,file = "nasal_intersecting.RData")

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
varianttable <- cmpfun(varianttable)
nasal.vartables <- list()
for(i in seq_along(nasal_intersecting)){
  nasal.vartables[[i]] <- varianttable(nasal_intersecting[[i]],everythinguser_nas[[i]])
}
names(nasal.vartables) <- names(sleekuser_nas)

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





