library(tidyverse)
library(rio)
library(rlist)
library(gtools)
library(here)

#saliva
cov.names <- list.files("coverage_summaries/")
cov.files <- paste0("coverage_summaries/",cov.names)
cov.files <- lapply(cov.files,import)
for(i in seq_along(cov.files)){
  cov.files[[i]] <- cov.files[[i]] %>%
    filter(chrom == "total") %>%
    select(mean)
}
cov.names <- gsub(".mosdepth.summary.txt","",cov.names)
cov.names <- gsub(".trim.genome","",cov.names)
names(cov.files) <- cov.names
save(cov.files,file = "saliva_coverage_depths.RData")

cov.means <- unlist(cov.files)
total.mean <- mean(cov.means)
total.med <- median(cov.means)
total.sd <- sd(cov.means)
cov.tibble <- as_tibble(cov.means)
colnames(cov.tibble) <- "means"
ggplot(cov.tibble,aes(x=means))+
  geom_histogram()+
  scale_x_discrete(name = "means", limits = c(0,5000,10000,15000,20000,30000,
  40000,50000,60000))

#nasal
nas.cov.names <- list.files("nasal_coverage_summaries/")
nas.cov.names <- mixedsort(sort(nas.cov.names))
nas.cov.files <- paste0("nasal_coverage_summaries/",nas.cov.names)
nas.cov.files <- lapply(nas.cov.files,import)
names(nas.cov.files) <- nas.cov.names
for(i in seq_along(nas.cov.files)){
  nas.cov.files[[i]] <- nas.cov.files[[i]]$median_coverage$NC_045512.2
}
nas.cov.names <- gsub("_report_metrics.json","",nas.cov.names)
names(nas.cov.files) <- nas.cov.names
save(nas.cov.files,file = "nasal_coverage_depths.RData")

nas.cov.meds <- unlist(nas.cov.files)
nas.cov.tibble <- as_tibble(nas.cov.meds)
colnames(nas.cov.tibble) <- "meds"
ggplot(nas.cov.tibble,aes(x=meds))+
  geom_histogram()



#testing testing
saliva.dirs <- list.dirs("saliva_output_tables")[-1]
names(saliva.dirs) <- names(nasal.dirs)
saliva.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir,pattern="ivar"))
  files <- lapply(full, import)
  return(files)
}
saliva.paths <- function(dir){
  full <- list.files(dir,pattern="ivar")
  name <- gsub(".ivar.tsv","",full)
  return(name)
}
saliva.names <- lapply(saliva.dirs,saliva.paths)
saliva.files <- lapply(saliva.dirs,saliva.read)

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]]$REGION <- saliva.names[[i]][[j]]
  }
}

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]]%>%
      rename("ID" = REGION)
  }
}

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

test <- saliva.files[[2]][[2]]

#testing 2
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
      filter(COVERAGE >= 500)
  }
}

test <- nasal.files[[1]][[2]]
nasal.clean <- function(file){
  file <- file %>%
    filter(ALT_FREQ >= .03) %>% 
    filter(DP >= 500) %>%
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
}
for(i in seq_along(nasal.files)){
  nasal.files[[i]] <- lapply(nasal.files[[i]],nasal.clean)
}


