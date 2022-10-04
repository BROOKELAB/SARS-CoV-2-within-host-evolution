library(tidyverse)
library(rio)
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
    dplyr::rename("POS"=V2)%>%
    dplyr::rename("REF"=V4)%>%
    dplyr::rename("ALT"=V5)%>%
    dplyr::rename("INFO1"=V8)%>%
    dplyr::rename("INFO2"=V10)%>%
    mutate("ALT_FREQ"=NA)%>%
    mutate("DP"=NA)
  for(i in seq_along(file$INFO1)){
    file$DP[[i]] <- strsplit(file$INFO1[[i]],";")[[1]][1]
    file$DP[[i]] <- as.numeric(strsplit(file$DP[[i]],"=")[[1]][2])
    file$ALT_FREQ[[i]] <- strsplit(file$INFO2[[i]],":")[[1]][4]
  }
  file$DP <- as.numeric(file$DP)
  file <- file %>%
    select("ID","POS","REF","ALT","ALT_FREQ","DP")%>%
    mutate("SNP"=paste0(REF,POS,ALT))%>%
    relocate(SNP,.before = ALT_FREQ)%>%
    distinct()
  return(file)
}

nasal.names <- lapply(nasal.dirs,nasal.paths)
nasal.files <- lapply(nasal.dirs,nasal.read)
for(i in seq_along(nasal.files)){
  nasal.files[[i]] <- lapply(nasal.files[[i]],nasal.read2)
}

for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    nasal.files[[i]][[j]]$ID <- nasal.names[[i]][[j]]
  }
}

load("nasal_coverage_depths.RData")
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    for(x in seq_along(nas.cov.files)){
      if(names(nas.cov.files)[[x]] == nasal.files[[i]][[j]]$ID[[1]]){
        nasal.files[[i]][[j]] <- nasal.files[[i]][[j]]%>%
          mutate("COVERAGE" = nas.cov.files[[x]])
      }
    }
  }
}
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    nasal.files[[i]][[j]] <- nasal.files[[i]][[j]]%>%
      filter(COVERAGE >= 200) #this can change
  }
}

nasal.clean <- function(file){
  file <- file %>%
    filter(ALT_FREQ >= .01) %>% 
    #filter(DP >= 200) %>%
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
    filter(POS != 15965) %>%
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
  nasal.files[[i]] <- lapply(nasal.files[[i]],nasal.clean)
}

saliva.dirs <- list.dirs("saliva_output_tables")[-1]
names(saliva.dirs) <- names(nasal.dirs)
saliva.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir,pattern="ivar"))
  files <- lapply(full, import)
  return(files)
}
saliva.files <- lapply(saliva.dirs,saliva.read)
saliva.paths <- function(dir){
  full <- list.files(dir,pattern="ivar")
  name <- gsub(".ivar.tsv","",full)
  return(name)
}
saliva.names <- lapply(saliva.dirs,saliva.paths)

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]]$REGION <- saliva.names[[i]][[j]]
  }
}

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]]%>%
      dplyr::rename("ID" = REGION)
  }
}

load("saliva_coverage_depths.RData")
for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    for(x in seq_along(cov.files)){
      if(names(cov.files)[[x]] == saliva.files[[i]][[j]]$ID[[1]]){
        saliva.files[[i]][[j]] <- saliva.files[[i]][[j]]%>%
          mutate("COVERAGE" = cov.files[[x]]$mean)
      }
    }
  }
}

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]] %>%
      filter(COVERAGE >= 1000)
  }
}
saliva.clean <- function(file){
  file <- file %>%
    dplyr::rename("DP"=TOTAL_DP)%>%
    select(POS,REF,ALT,ALT_FREQ,DP)
  file$DP <- as.numeric(file$DP)
  file <- file %>%
    #filter(DP >= 500)%>%
    filter(ALT_FREQ >= .01) %>% 
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
    filter(POS != 15965) %>%
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
    filter(POS != 29700) %>%
    distinct()
  return(file)
}
for(i in seq_along(saliva.files)){
  saliva.files[[i]] <- lapply(saliva.files[[i]],saliva.clean)
}

saliva.reformat <- function(saliva.user){
  for(i in seq_along(saliva.user)){
    for(j in seq_along(saliva.user[[i]]$ALT)){
      if(substr(saliva.user[[i]]$ALT[[j]],0,1)=="+"){
        saliva.user[[i]]$ALT[[j]] <- paste0(saliva.user[[i]]$REF[[j]], 
                                             substr(saliva.user[[i]]$ALT[[j]],2,nchar(saliva.user[[i]]$ALT[[j]])))
      }
      if(substr(saliva.user[[i]]$ALT[[j]],0,1)=="-"){
        saliva.user[[i]]$REF[[j]] <- paste0(saliva.user[[i]]$REF[[j]],
                                             substr(saliva.user[[i]]$ALT[[j]],2,nchar(saliva.user[[i]]$ALT[[j]])))
        saliva.user[[i]]$ALT[[j]] <- substr(saliva.user[[i]]$REF[[j]],0,1)
      }
    }
  }
  return(saliva.user)
}
saliva.files <- lapply(saliva.files,saliva.reformat)
for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]] %>% 
      mutate(SNP = paste0(REF,POS,ALT))%>%
      relocate(SNP,.before=ALT_FREQ)
  }
}

#cutting out files that exist in nasal but not saliva
#and vice versa
nasal.files$user_442978 <- nasal.files$user_442978[-4]
saliva.files$user_444332 <- saliva.files$user_444332[-c(2,8)]
nasal.files$user_444332 <- nasal.files$user_444332[-c(3,7,9:12)]
saliva.files$user_450241 <- saliva.files$user_450241[-c(4)]
saliva.files$user_451152 <- saliva.files$user_451152[-c(1,7,8)]
saliva.files$user_471467 <- saliva.files$user_471467[-c(2,7)]
saliva.files$user_471588 <- saliva.files$user_471588[-c(1,6)]

#cutting out the missing files from saliva set
#these numbers are all adjusted to account for deletion i did above
saliva.files$user_432686 <- saliva.files$user_432686[-c(5,6)]
saliva.files$user_433227 <- saliva.files$user_433227[-1]
saliva.files$user_450241 <- saliva.files$user_450241[-c(3,5)]
saliva.files$user_451709 <- saliva.files$user_451709[-7]
saliva.files$user_453058 <- saliva.files$user_453058[-c(2,6)]
saliva.files$user_459597 <- saliva.files$user_459597[-c(1,2)]
saliva.files$user_471467 <- saliva.files$user_471467[-c(1,7,8)]
saliva.files$user_471588 <- saliva.files$user_471588[-c(5:7)]

#remove empty files that have been excluded for low coverage
nasal.files$user_432686 <- nasal.files$user_432686[-c(2)] #2 nasal
saliva.files$user_432686 <- saliva.files$user_432686[-c(2)]
nasal.files$user_433227 <- nasal.files$user_433227[-c(1,6)] #all saliva
saliva.files$user_433227 <- saliva.files$user_433227[-c(1,6)] 
nasal.files$user_438577 <- nasal.files$user_438577[-c(2,6)] #all saliva
saliva.files$user_438577 <- saliva.files$user_438577[-c(2,6)] 
nasal.files$user_450241 <- nasal.files$user_450241[-c(2,3)] #nasal 2 nasal 3
saliva.files$user_450241 <- saliva.files$user_450241[-c(2,3)]
nasal.files$user_451709 <- nasal.files$user_451709[-c(6)] #saliva 6
saliva.files$user_451709 <- saliva.files$user_451709[-c(6)]
nasal.files$user_453058 <- nasal.files$user_453058[-c(1,4)] #nasal 1 saliva 4
saliva.files$user_453058 <- saliva.files$user_453058[-c(1,4)]
nasal.files$user_459597 <- nasal.files$user_459597[-c(4)] #saliva 4
saliva.files$user_459597 <- saliva.files$user_459597[-c(4)]
nasal.files$user_471588 <- nasal.files$user_471588[-c(4)] #nasal 4
saliva.files$user_471588 <- saliva.files$user_471588[-c(4)]
nasal.files$user_475670 <- nasal.files$user_475670[-c(3)] #saliva 3
saliva.files$user_475670 <- saliva.files$user_475670[-c(3)]

#get rid of day 2 and day 5 from 451709 for potential contam
#day 6 has been removed
saliva.files$user_451709 <- saliva.files$user_451709[-c(2,5)]
nasal.files$user_451709 <- nasal.files$user_451709[-c(2,5)]

#pull intersecting mutations for nasal samples
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
nasal.pos <- lapply(nasal.files,positions)

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

nas_432686 <- list("time1"=var1(nasal.pos[[1]]),"time2"=var2(nasal.pos[[1]]),
                   "time3"=var3(nasal.pos[[1]]))
nas_433227 <- list("time1"=var1(nasal.pos[[2]]),"time2"=var2(nasal.pos[[2]]),
                   "time3"=var3(nasal.pos[[2]]),"time4"=var4(nasal.pos[[2]]))
nas_438577 <- list("time1"=var1(nasal.pos[[3]]),"time2"=var2(nasal.pos[[3]]),
                   "time3"=var3(nasal.pos[[3]]),"time4"=var4(nasal.pos[[3]]))
nas_442978 <- list("time1"=var1(nasal.pos[[4]]),"time2"=var2(nasal.pos[[4]]),
                   "time3"=var3(nasal.pos[[4]]),"time4"=var4(nasal.pos[[4]]),
                   "time5"=var5(nasal.pos[[4]]))
nas_444332 <- list("time1"=var1(nasal.pos[[5]]),"time2"=var2(nasal.pos[[5]]),
                   "time3"=var3(nasal.pos[[5]]),"time4"=var4(nasal.pos[[5]]),
                   "time5"=var5(nasal.pos[[5]]),"time6"=var6(nasal.pos[[5]]))
nas_444633 <- list("time1"=var1(nasal.pos[[6]]),"time2"=var2(nasal.pos[[6]]),
                   "time3"=var3(nasal.pos[[6]]),"time4"=var4(nasal.pos[[6]]),
                   "time5"=var5(nasal.pos[[6]]))
nas_450241 <- list("time1"=var1(nasal.pos[[7]]))
nas_451152 <- list("time1"=var1(nasal.pos[[8]]),"time2"=var2(nasal.pos[[8]]),
                   "time3"=var3(nasal.pos[[8]]),"time4"=var4(nasal.pos[[8]]),
                   "time5"=var5(nasal.pos[[8]]))
nas_451709 <- list("time1"=var1(nasal.pos[[9]]),"time2"=var2(nasal.pos[[9]]),
                   "time3"=var3(nasal.pos[[9]]),"time4"=var4(nasal.pos[[9]]),
                   "time5"=var5(nasal.pos[[9]]))
nas_453058 <- list("time1"=var1(nasal.pos[[10]]),"time2"=var2(nasal.pos[[10]]),
                   "time3"=var3(nasal.pos[[10]]))
nas_459597 <- list("time1"=var1(nasal.pos[[11]]),"time2"=var2(nasal.pos[[11]]),
                   "time3"=var3(nasal.pos[[11]]),"time4"=var4(nasal.pos[[11]]))
nas_471467 <- list("time1"=var1(nasal.pos[[12]]),"time2"=var2(nasal.pos[[12]]),
                   "time3"=var3(nasal.pos[[12]]),"time4"=var4(nasal.pos[[12]]),
                   "time5"=var5(nasal.pos[[12]]))
nas_471588 <- list("time1"=var1(nasal.pos[[13]]),"time2"=var2(nasal.pos[[13]]),
                   "time3"=var3(nasal.pos[[13]]))
nas_475670 <- list("time1"=var1(nasal.pos[[14]]),"time2"=var2(nasal.pos[[14]]),
                   "time3"=var3(nasal.pos[[14]]),"time4"=var4(nasal.pos[[14]]))
nas_list <- list("nas_432686"=nas_432686,"nas_433227"=nas_433227,"nas_438577"=nas_438577,
                 "nas_442978"=nas_442978,"nas_444332"=nas_444332,"nas_444633"=nas_444633,
                 "nas_450241"=nas_450241,"nas_451152"=nas_451152,"nas_451709"=nas_451709,
                 "nas_453058"=nas_453058,"nas_459597"=nas_459597,"nas_471467"=nas_471467,
                 "nas_471588"=nas_471588,"nas_475670"=nas_475670)
selfevict <- function(var){
  for(i in 1:length(var)){
    var[[i]][i] <- NA
  }
  return(var)
}
nas_list <- lapply(nas_list,selfevict)

#need to remove nas_450241 because it only has one timepoint
nas_list <- nas_list[-7]

keepcommon <- function(var){
  flat <- flatten(var)
  splat <- na.omit(rlist::list.rbind(flat))
  common <- distinct(splat)
  return(common)
}
nasal_intersecting <- lapply(nas_list,keepcommon)
names(nasal_intersecting) <- names(nas_list)

nasal_intersecting <- lapply(nasal_intersecting,arrange,POS) 
for(i in seq_along(nasal_intersecting)){
  nasal_intersecting[[i]] <- nasal_intersecting[[i]] %>%
    mutate("SNP" = paste0(REF,POS,ALT))
}

#pull intersecting mutations for saliva samples
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      select(POS,REF,ALT)
  }
  return(allpos)
}
saliva.pos <- lapply(saliva.files,positions)

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

sal_432686 <- list("time1"=var1(saliva.pos[[1]]),"time2"=var2(saliva.pos[[1]]),
                   "time3"=var3(saliva.pos[[1]]))
sal_433227 <- list("time1"=var1(saliva.pos[[2]]),"time2"=var2(saliva.pos[[2]]),
                   "time3"=var3(saliva.pos[[2]]),"time4"=var4(saliva.pos[[2]]))
sal_438577 <- list("time1"=var1(saliva.pos[[3]]),"time2"=var2(saliva.pos[[3]]),
                   "time3"=var3(saliva.pos[[3]]),"time4"=var4(saliva.pos[[3]]))
sal_442978 <- list("time1"=var1(saliva.pos[[4]]),"time2"=var2(saliva.pos[[4]]),
                   "time3"=var3(saliva.pos[[4]]),"time4"=var4(saliva.pos[[4]]),
                   "time5"=var5(saliva.pos[[4]]))
sal_444332 <- list("time1"=var1(saliva.pos[[5]]),"time2"=var2(saliva.pos[[5]]),
                   "time3"=var3(saliva.pos[[5]]),"time4"=var4(saliva.pos[[5]]),
                   "time5"=var5(saliva.pos[[5]]),"time6"=var6(saliva.pos[[5]]))
sal_444633 <- list("time1"=var1(saliva.pos[[6]]),"time2"=var2(saliva.pos[[6]]),
                   "time3"=var3(saliva.pos[[6]]),"time4"=var4(saliva.pos[[6]]),
                   "time5"=var5(saliva.pos[[6]]))
sal_450241 <- list("time1"=var1(saliva.pos[[7]]))
sal_451152 <- list("time1"=var1(saliva.pos[[8]]),"time2"=var2(saliva.pos[[8]]),
                   "time3"=var3(saliva.pos[[8]]),"time4"=var4(saliva.pos[[8]]),
                   "time5"=var5(saliva.pos[[8]]))
sal_451709 <- list("time1"=var1(saliva.pos[[9]]),"time2"=var2(saliva.pos[[9]]),
                   "time3"=var3(saliva.pos[[9]]),"time4"=var4(saliva.pos[[9]]),
                   "time5"=var5(saliva.pos[[9]]))
sal_453058 <- list("time1"=var1(saliva.pos[[10]]),"time2"=var2(saliva.pos[[10]]),
                   "time3"=var3(saliva.pos[[10]]))
sal_459597 <- list("time1"=var1(saliva.pos[[11]]),"time2"=var2(saliva.pos[[11]]),
                   "time3"=var3(saliva.pos[[11]]),"time4"=var4(saliva.pos[[11]]))
sal_471467 <- list("time1"=var1(saliva.pos[[12]]),"time2"=var2(saliva.pos[[12]]),
                   "time3"=var3(saliva.pos[[12]]),"time4"=var4(saliva.pos[[12]]),
                   "time5"=var5(saliva.pos[[12]]))
sal_471588 <- list("time1"=var1(saliva.pos[[13]]),"time2"=var2(saliva.pos[[13]]),
                   "time3"=var3(saliva.pos[[13]]))
sal_475670 <- list("time1"=var1(saliva.pos[[14]]),"time2"=var2(saliva.pos[[14]]),
                   "time3"=var3(saliva.pos[[14]]),"time4"=var4(saliva.pos[[14]]))

sal_list <- list("sal_432686"=sal_432686,"sal_433227"=sal_433227,"sal_438577"=sal_438577,
                 "sal_442978"=sal_442978,"sal_444332"=sal_444332,"sal_444633"=sal_444633,
                 "sal_450241"=sal_450241,"sal_451152"=sal_451152,"sal_451709"=sal_451709,
                 "sal_453058"=sal_453058,"sal_459597"=sal_459597,"sal_471467"=sal_471467,
                 "sal_471588"=sal_471588,"sal_475670"=sal_475670)
selfevict <- function(var){
  for(i in 1:length(var)){
    var[[i]][i] <- NA
  }
  return(var)
}
sal_list <- lapply(sal_list,selfevict)

#need to remove sal_450241 because it only has one timepoint
sal_list <- sal_list[-7]

keepcommon <- function(var){
  flat <- flatten(var)
  splat <- na.omit(rlist::list.rbind(flat))
  common <- distinct(splat)
  return(common)
}
saliva_intersecting <- lapply(sal_list,keepcommon)
names(saliva_intersecting) <- names(sal_list)
saliva_intersecting <- lapply(saliva_intersecting,arrange,POS) #450241 has none
for(i in seq_along(saliva_intersecting)){
  saliva_intersecting[[i]] <- saliva_intersecting[[i]] %>%
    mutate("SNP" = paste0(REF,POS,ALT))
}

#pare down both sets of files to include just intersecting SNPs
#first remove 450241 from nasal files
nasal.files <- nasal.files[-7]
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    nasal.files[[i]][[j]] <- nasal.files[[i]][[j]] %>%
      mutate(KEEP=NA)
    for(k in seq_along(nasal.files[[i]][[j]]$SNP)){
      if(nasal.files[[i]][[j]]$SNP[[k]] %in% nasal_intersecting[[i]]$SNP){
        nasal.files[[i]][[j]]$KEEP[[k]] <- "YES"
      }else{
        nasal.files[[i]][[j]]$KEEP[[k]] <- "NO"
      }
    }
  }
}

for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    nasal.files[[i]][[j]] <- nasal.files[[i]][[j]]%>%
      filter(KEEP=="YES")
  }
}

#remove 450241 from saliva files
saliva.files <- saliva.files[-7]
for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]] %>%
      mutate(KEEP=NA)
    for(k in seq_along(saliva.files[[i]][[j]]$SNP)){
      if(saliva.files[[i]][[j]]$SNP[[k]] %in% saliva_intersecting[[i]]$SNP){
        saliva.files[[i]][[j]]$KEEP[[k]] <- "YES"
      }else{
        saliva.files[[i]][[j]]$KEEP[[k]] <- "NO"
      }
    }
  }
}

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]]%>%
      filter(KEEP=="YES")
  }
}


#label high freq SNVs as 1
for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    for(k in seq_along(saliva.files[[i]][[j]]$POS)){
      if(saliva.files[[i]][[j]]$ALT_FREQ[[k]] >= .99){
        saliva.files[[i]][[j]]$ALT_FREQ[[k]] <- 1
      }
    }
  }
}
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    for(k in seq_along(nasal.files[[i]][[j]]$POS)){
      if(nasal.files[[i]][[j]]$ALT_FREQ[[k]] >= .99){
        nasal.files[[i]][[j]]$ALT_FREQ[[k]] <- 1
      }
    }
  }
}


#label low dp as X
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    for(k in seq_along(nasal.files[[i]][[j]]$POS))
    if(nasal.files[[i]][[j]]$DP[[k]] < 200){
      nasal.files[[i]][[j]]$ALT_FREQ[[k]] <- "X"
    }
  }
}
for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    for(k in seq_along(saliva.files[[i]][[j]]$POS))
      if(saliva.files[[i]][[j]]$DP[[k]] < 500){
        saliva.files[[i]][[j]]$ALT_FREQ[[k]] <- "X"
      }
  }
}

save(nasal.files, file = "nasal_files.RData")
save(saliva.files, file = "saliva_files.RData")

#comparison tables
env.comp <- function(nasal,saliva){
  comp <- list()
  for(i in seq_along(saliva)){
    saliva[[i]]$ALT_FREQ <- as.character(saliva[[i]]$ALT_FREQ)
  }
  for(i in seq_along(nasal)){
    comp[[i]] <- full_join(nasal[[i]],saliva[[i]],by=c("POS","REF","ALT","SNP"))
  }
  return(comp)
}
comp.tables <- list()
for(i in seq_along(nasal.files)){
  comp.tables[[i]] <- env.comp(nasal.files[[i]],saliva.files[[i]])
}

for(i in seq_along(comp.tables)){
  for(j in seq_along(comp.tables[[i]])){
    comp.tables[[i]][[j]] <- comp.tables[[i]][[j]] %>%
      dplyr::rename("NASAL_FREQ"=ALT_FREQ.x)%>%
      dplyr::rename("SALIVA_FREQ"=ALT_FREQ.y)%>%
      select(POS,REF,ALT,SNP,NASAL_FREQ,SALIVA_FREQ)
  }
}

comp.list <- list()
for(i in seq_along(comp.tables)){
  comp.list[[i]] <- comp.tables[[i]] %>%
    purrr::reduce(left_join, by=c("POS","REF","ALT","SNP"))
}

vec <- list()
vec2 <- list()
recols <- list()
for(i in seq_along(comp.list)){
  vec[[i]] <- 1:((length(comp.list[[i]])-4)/2)
  vec2[[i]] <- rep(vec[[i]],each=2)
  recols[[i]] <- colnames(comp.list[[i]])[-c(1:4)]
  recols[[i]] <- gsub(".x","",recols[[i]])
  recols[[i]] <- gsub(".y","",recols[[i]])
  recols[[i]] <- paste0(recols[[i]],"_",vec2[[i]])
  colnames(comp.list[[i]])[-c(1:4)] <- recols[[i]]
}
names(comp.list) <- names(nasal.dirs)[-7]

#replace NA values (due to freq) with 0
#replace X values (due to dp) with NA
for(i in seq_along(comp.list)){
  comp.list[[i]] <- comp.list[[i]] %>%
    replace(is.na(.),0)
  comp.list[[i]][comp.list[[i]]=="X"] <- NA
}


#match NA values between environments

freq1 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_1[[i]])){
      comptable$SALIVA_FREQ_1[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_1[[i]])){
      comptable$NASAL_FREQ_1[[i]] <- NA
    }
  }
  return(comptable)
}
freq2 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_2[[i]])){
      comptable$SALIVA_FREQ_2[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_2[[i]])){
      comptable$NASAL_FREQ_2[[i]] <- NA
    }
  }
  return(comptable)
}
freq3 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_3[[i]])){
      comptable$SALIVA_FREQ_3[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_3[[i]])){
      comptable$NASAL_FREQ_3[[i]] <- NA
    }
  }
  return(comptable)
}
freq4 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_4[[i]])){
      comptable$SALIVA_FREQ_4[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_4[[i]])){
      comptable$NASAL_FREQ_4[[i]] <- NA
    }
  }
  return(comptable)
} 
freq5 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_5[[i]])){
      comptable$SALIVA_FREQ_5[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_5[[i]])){
      comptable$NASAL_FREQ_5[[i]] <- NA
    }
  }
  return(comptable)
}
freq6 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_6[[i]])){
      comptable$SALIVA_FREQ_6[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_6[[i]])){
      comptable$NASAL_FREQ_6[[i]] <- NA
    }
  }
  return(comptable)
} 
freq7 <- function(comptable){
  for(i in seq_along(comptable$POS)){
    if(is.na(comptable$NASAL_FREQ_7[[i]])){
      comptable$SALIVA_FREQ_7[[i]] <- NA
    }
    if(is.na(comptable$SALIVA_FREQ_7[[i]])){
      comptable$NASAL_FREQ_7[[i]] <- NA
    }
  }
  return(comptable)
}

#432686
comp.list$user_432686 <- freq1(comp.list$user_432686)
comp.list$user_432686 <- freq2(comp.list$user_432686)
comp.list$user_432686 <- freq3(comp.list$user_432686)
#433227
comp.list$user_433227 <- freq1(comp.list$user_433227)
comp.list$user_433227 <- freq2(comp.list$user_433227)
comp.list$user_433227 <- freq3(comp.list$user_433227)
comp.list$user_433227 <- freq4(comp.list$user_433227)
#438577
comp.list$user_438577 <- freq1(comp.list$user_438577)
comp.list$user_438577 <- freq2(comp.list$user_438577)
comp.list$user_438577 <- freq3(comp.list$user_438577)
comp.list$user_438577 <- freq4(comp.list$user_438577)
#442978
comp.list$user_442978 <- freq1(comp.list$user_442978)
comp.list$user_442978 <- freq2(comp.list$user_442978)
comp.list$user_442978 <- freq3(comp.list$user_442978)
comp.list$user_442978 <- freq4(comp.list$user_442978)
comp.list$user_442978 <- freq5(comp.list$user_442978)
#444332
comp.list$user_444332 <- freq1(comp.list$user_444332)
comp.list$user_444332 <- freq2(comp.list$user_444332)
comp.list$user_444332 <- freq3(comp.list$user_444332)
comp.list$user_444332 <- freq4(comp.list$user_444332)
comp.list$user_444332 <- freq5(comp.list$user_444332)
comp.list$user_444332 <- freq6(comp.list$user_444332)
#444633
comp.list$user_444633 <- freq1(comp.list$user_444633)
comp.list$user_444633 <- freq2(comp.list$user_444633)
comp.list$user_444633 <- freq3(comp.list$user_444633)
comp.list$user_444633 <- freq4(comp.list$user_444633)
comp.list$user_444633 <- freq5(comp.list$user_444633)
#451152
comp.list$user_451152 <- freq1(comp.list$user_451152)
comp.list$user_451152 <- freq2(comp.list$user_451152)
comp.list$user_451152 <- freq3(comp.list$user_451152)
comp.list$user_451152 <- freq4(comp.list$user_451152)
comp.list$user_451152 <- freq5(comp.list$user_451152)
#451709
comp.list$user_451709 <- freq1(comp.list$user_451709)
comp.list$user_451709 <- freq2(comp.list$user_451709)
comp.list$user_451709 <- freq3(comp.list$user_451709)
comp.list$user_451709 <- freq4(comp.list$user_451709)
comp.list$user_451709 <- freq5(comp.list$user_451709)
#453058
comp.list$user_453058 <- freq1(comp.list$user_453058)
comp.list$user_453058 <- freq2(comp.list$user_453058)
comp.list$user_453058 <- freq3(comp.list$user_453058)
#459597
comp.list$user_459597 <- freq1(comp.list$user_459597)
comp.list$user_459597 <- freq2(comp.list$user_459597)
comp.list$user_459597 <- freq3(comp.list$user_459597)
comp.list$user_459597 <- freq4(comp.list$user_459597)
#471467
comp.list$user_471467 <- freq1(comp.list$user_471467)
comp.list$user_471467 <- freq2(comp.list$user_471467)
comp.list$user_471467 <- freq3(comp.list$user_471467)
comp.list$user_471467 <- freq4(comp.list$user_471467)
comp.list$user_471467 <- freq5(comp.list$user_471467)
#471588
comp.list$user_471588 <- freq1(comp.list$user_471588)
comp.list$user_471588 <- freq2(comp.list$user_471588)
comp.list$user_471588 <- freq3(comp.list$user_471588)
#475670
comp.list$user_475670 <- freq1(comp.list$user_475670)
comp.list$user_475670 <- freq2(comp.list$user_475670)
comp.list$user_475670 <- freq3(comp.list$user_475670)
comp.list$user_475670 <- freq4(comp.list$user_475670)

save(comp.list, file = "comparison_list.RData")

#replace NA with NaN (to read into matlab)
for(i in seq_along(comp.list)){
  comp.list[[i]][is.na(comp.list[[i]])] <- NaN
}

filenames <- paste0("comparison_tables/",names(comp.list),".csv")
for(i in seq_along(comp.list)){
  write.csv(comp.list[[i]],filenames[[i]])
}

for(i in seq_along(comp.tables)){
  for(j in seq_along(comp.tables[[i]])){
    for(x in seq_along(comp.tables[[i]][[j]]$POS)){
      if(is.na(comp.tables[[i]][[j]]$NASAL_FREQ[[x]])){
        comp.tables[[i]][[j]]$NASAL_FREQ[[x]] <- 0
      }
      if(is.na(comp.tables[[i]][[j]]$SALIVA_FREQ[[x]])){
        comp.tables[[i]][[j]]$SALIVA_FREQ[[x]] <- 0
      }
    }
  }
}

for(i in seq_along(comp.tables)){
  for(j in seq_along(comp.tables[[i]])){
    for(x in seq_along(comp.tables[[i]][[j]]$POS)){
      if(comp.tables[[i]][[j]]$NASAL_FREQ[[x]] == "X"){
        comp.tables[[i]][[j]]$NASAL_FREQ[[x]] <- NA
      }
      if(comp.tables[[i]][[j]]$SALIVA_FREQ[[x]] == "X"){
        comp.tables[[i]][[j]]$SALIVA_FREQ[[x]] <- NA
      }
    }
  }
}


unlist.comp <- list()
for(i in seq_along(comp.tables)){
  unlist.comp[[i]] <- rlist::list.rbind(comp.tables[[i]])
}
unlist.comp <- rlist::list.rbind(unlist.comp)
unlist.comp <- select(unlist.comp,NASAL_FREQ,SALIVA_FREQ)
unlist.comp <- na.omit(unlist.comp)
write.csv(unlist.comp,"nasal_saliva_freqs.csv")

#distribution of saliva iSNVs missing from nasal
load("comparison_list.RData")


missing.distribution <- function(user){
  saliva.only <- user[,grep("SALIVA",colnames(user))]
  nasal.only <- user[,grep("NASAL",colnames(user))]
  for(i in seq_along(colnames(saliva.only))){
    saliva.only[,i] <- coalesce(saliva.only[,i],"X")
    nasal.only[,i] <- coalesce(nasal.only[,i],"X")
  }
  
  dis.matrix <- as.data.frame(matrix(nrow = nrow(saliva.only), 
                                     ncol = ncol(saliva.only),
                                     data = NA))
  colnames(dis.matrix) <- colnames(saliva.only)
  
  
  for(i in seq_along(rownames(saliva.only))){
    for(j in seq_along(colnames(saliva.only))){
      if(nasal.only[i,j] == 0){
        dis.matrix[i,j] <- saliva.only[i,j]
      }
    }
  }
  dis.vec <- na.omit(as.numeric(unlist(flatten(dis.matrix))))
  return(dis.vec)
}

distributions <- lapply(comp.list,missing.distribution)
distributions <- unlist(distributions)
#hist(unlist(distributions))
test <- as.data.frame(distributions)
  
ggplot(test, aes(x = distributions))+
  geom_histogram(aes(y=stat(density)), bins =100)+
  geom_density()+
  xlab("Frequency")+
  ylab("Density")+
  ggtitle("Frequency distribution of saliva iSNVs \nmissing from nasal compartment")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 19),
          plot.title = element_text(size = 22))
ggsave("figs/saliva_distribution.png")




