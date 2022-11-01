library(tidyverse)
library(rio)
library(gtools)
library(here)

nasal.dirs <- list.dirs("nasal_output_tables_cut/")[-1]
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
    filter(ALT_FREQ >= .03) %>% 
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
} #contains frequency cutoffs

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
    select(POS,REF,ALT,ALT_FREQ,DP,COVERAGE)
  file$DP <- as.numeric(file$DP)
  file <- file %>%
    #filter(DP >= 500)%>%
    filter(ALT_FREQ >= .03) %>% 
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

#remove empty files that have been excluded for low coverage
saliva.files$user_432686 <- saliva.files$user_432686[-2]
nasal.files$user_432686 <- nasal.files$user_432686[-2]
saliva.files$user_433227 <- saliva.files$user_433227[-c(1,6)]
nasal.files$user_433227 <- nasal.files$user_433227[-c(1,6)]
saliva.files$user_438577 <- saliva.files$user_438577[-c(2,6)]
nasal.files$user_438577 <- nasal.files$user_438577[-c(2,6)]
saliva.files$user_450241 <- saliva.files$user_450241[-c(2,3)]
nasal.files$user_450241 <- nasal.files$user_450241[-c(2,3)]
saliva.files$user_451709 <- saliva.files$user_451709[-6]
nasal.files$user_451709 <- nasal.files$user_451709[-6]
saliva.files$user_453058 <- saliva.files$user_453058[-c(1,4)]
nasal.files$user_453058 <- nasal.files$user_453058[-c(1,4)]
saliva.files$user_459597 <- saliva.files$user_459597[-4]
nasal.files$user_459597 <- nasal.files$user_459597[-4]
saliva.files$user_471588 <- saliva.files$user_471588[-4]
nasal.files$user_471588 <- nasal.files$user_471588[-4]

#get rid of times 2,5,7 from 451709 for potential contam
#point 6 has already been removed (7 is now 6)
saliva.files$user_451709 <- saliva.files$user_451709[-c(2,5,6)]
nasal.files$user_451709 <- nasal.files$user_451709[-c(2,5,6)]

#get rid of user 450241 because only 1 time pt
saliva.files <- saliva.files[-7]
nasal.files <- nasal.files[-7]

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

nasal.intersecting <- lapply(nasal.pos,snp.intersect)
nasal.intersecting <- lapply(nasal.intersecting,arrange,POS)

names(nasal.intersecting) <- names(nasal.files)

for(i in seq_along(nasal.intersecting)){
  nasal.intersecting[[i]] <- nasal.intersecting[[i]] %>%
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

saliva.intersecting <- lapply(saliva.pos,snp.intersect)
saliva.intersecting <- lapply(saliva.intersecting,arrange,POS)

names(saliva.intersecting) <- names(saliva.files)

for(i in seq_along(saliva.intersecting)){
  saliva.intersecting[[i]] <- saliva.intersecting[[i]] %>%
    mutate("SNP" = paste0(REF,POS,ALT))
}

#pare down both sets of files to include just intersecting SNPs
for(i in seq_along(nasal.files)){
  for(j in seq_along(nasal.files[[i]])){
    nasal.files[[i]][[j]] <- nasal.files[[i]][[j]] %>%
      mutate(KEEP=NA)
    for(k in seq_along(nasal.files[[i]][[j]]$SNP)){
      if(nasal.files[[i]][[j]]$SNP[[k]] %in% nasal.intersecting[[i]]$SNP){
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

for(i in seq_along(saliva.files)){
  for(j in seq_along(saliva.files[[i]])){
    saliva.files[[i]][[j]] <- saliva.files[[i]][[j]] %>%
      mutate(KEEP=NA)
    for(k in seq_along(saliva.files[[i]][[j]]$SNP)){
      if(saliva.files[[i]][[j]]$SNP[[k]] %in% saliva.intersecting[[i]]$SNP){
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

#remove consensus variants
remove.consensus <- function(user){
  info <- user[,1:4]
  freqs <- user[,grep("FREQ",colnames(user))]
  for(i in seq_along(rownames(freqs))){
    if(all(freqs[i,] > .97 | is.na(freqs[i,]))){
      freqs[i,] <- "X"
    }
  }
  final <- cbind(info,freqs)
  final <- final %>%
    filter(NASAL_FREQ_1 != "X" | is.na(NASAL_FREQ_1))
  return(final)
}
comp.list <- lapply(comp.list,remove.consensus)

save(comp.list, file = "comparison_list.RData")

#replace NA with NaN (to read into matlab)
for(i in seq_along(comp.list)){
  comp.list[[i]][is.na(comp.list[[i]])] <- NaN
}

filenames <- paste0("comparison_tables/",names(comp.list),".csv")
for(i in seq_along(comp.list)){
  write.csv(comp.list[[i]],filenames[[i]])
}

#frequency comparisons
freq.comp <- function(user){
  nasal.user <- user[,grep("NASAL",colnames(user))]
  nasal.user <- (pivot_longer(nasal.user, cols = everything()))[,2]
  colnames(nasal.user) <- "NASAL"
  saliva.user <- user[,grep("SALIVA",colnames(user))]
  saliva.user <- (pivot_longer(saliva.user, cols = everything()))[,2]
  colnames(saliva.user) <- "SALIVA"
  saliva.nasal <- bind_cols(saliva.user,nasal.user)
  return(saliva.nasal)
}

unlist.comp <- lapply(comp.list,freq.comp)
unlist.comp <- bind_rows(unlist.comp)
write.csv(unlist.comp,"nasal_saliva_freqs.csv")

#distribution of saliva iSNVs missing from nasal
load("comparison_list_unfiltered.RData") #version without frequency thresholds
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
distributions <- as.data.frame(distributions)
quantile(distributions$distributions)

ggplot(distributions, aes(x = distributions))+
  geom_histogram(aes(y=stat(density)), bins =100)+
  geom_density()+
  geom_vline(xintercept = 0.03, color = "red")+
  xlab("Frequency")+
  ylab("Density")+
  ggtitle("Frequency distribution of saliva iSNVs \nmissing from nasal compartment")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
          axis.text = element_text(size = 19),
          plot.title = element_text(size = 22))
ggsave("figs/saliva_distribution.png")




