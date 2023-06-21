library(tidyverse)
library(rio)
library(gtools)
library(here)

#read nasal files
nasal.dirs <- list.dirs("nasal_output_tables_cutsaliva")[-c(1,2)]
nasal.dirs26 <- list.dirs("nasal_output_tables_cutsaliva_cut26")[-c(1,2)]
names(nasal.dirs) <- c("user_432686","user_433227","user_438577","user_442978",
                       "user_444332","user_444633","user_451152","user_451709",
                       "user_453058","user_459597","user_471467","user_471588",
                       "user_475670")
names(nasal.dirs26) <- names(nasal.dirs)
nasal.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir, pattern = "ivar"))
  full <- mixedsort(sort(full)) #need this so it reads in date order
  files <- lapply(full, read.table,sep="\t",header =T)
  return(files)
}
nasal.files <- lapply(nasal.dirs26, nasal.read)
nasal.files.uncut <- lapply(nasal.dirs, nasal.read)

seq.filter <- function(file){
  file <- file %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP) %>%
    filter(ALT_FREQ >= 0.03)%>%
    #filter(TOTAL_DP >= 1000)%>% 
    #not filtering out depth yet because we want to distinguish between
    #iSNVs not present and iSNVs present but at a low depth
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
}
for(i in seq_along(nasal.files)){
  nasal.files[[i]] <- lapply(nasal.files[[i]], seq.filter)
  nasal.files[[i]] <- lapply(nasal.files[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  nasal.files[[i]] <- lapply(nasal.files[[i]], relocate, SNP, .after = ALT)
}
for(i in seq_along(nasal.files.uncut)){
  nasal.files.uncut[[i]] <- lapply(nasal.files.uncut[[i]], seq.filter)
  nasal.files.uncut[[i]] <- lapply(nasal.files.uncut[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  nasal.files.uncut[[i]] <- lapply(nasal.files.uncut[[i]], relocate, SNP, .after = ALT)
}

#read saliva files
saliva.dirs <- list.dirs("saliva_output_tables")[-1]
saliva.dirs26 <- list.dirs("saliva_output_tables_cut26")[-c(1,2)]
names(saliva.dirs) <- names(nasal.dirs)
names(saliva.dirs26) <- names(nasal.dirs)
saliva.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir,pattern="ivar"))
  files <- lapply(full, import)
  return(files)
}
saliva.files <- lapply(saliva.dirs26,saliva.read)
saliva.files.uncut <- lapply(saliva.dirs, saliva.read)

seq.filter <- function(file){
  file <- file %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP) %>%
    filter(ALT_FREQ >= 0.03)%>%
    #filter(TOTAL_DP >= 1000)%>% 
    #not filtering out depth yet because we want to distinguish between
    #iSNVs not present and iSNVs present but at a low depth
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
}
for(i in seq_along(saliva.files)){
  saliva.files[[i]] <- lapply(saliva.files[[i]], seq.filter)
  saliva.files[[i]] <- lapply(saliva.files[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  saliva.files[[i]] <- lapply(saliva.files[[i]], relocate, SNP, .after = ALT)
}
for(i in seq_along(saliva.files.uncut)){
  saliva.files.uncut[[i]] <- lapply(saliva.files.uncut[[i]], seq.filter)
  saliva.files.uncut[[i]] <- lapply(saliva.files.uncut[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  saliva.files.uncut[[i]] <- lapply(saliva.files.uncut[[i]], relocate, SNP, .after = ALT)
}

#pull intersecting iSNVs for nasal samples
#(iSNVs over 3% freq, present on at least 2 days in samples with Ct<26)
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      filter(ALT_FREQ >= 0.03) %>%
      filter(TOTAL_DP >= 1000) %>%
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
names(nasal.intersecting) <- gsub("user","nasal",names(nasal.files))

for(i in seq_along(nasal.intersecting)){
  nasal.intersecting[[i]] <- nasal.intersecting[[i]] %>%
    mutate("SNP" = paste0(REF,POS,ALT))
}

#pull intersecting mutations for saliva samples
positions <- function(sleekuser){
  allpos <- list()
  for (i in seq_along(sleekuser)){
    allpos[[i]] <- sleekuser[[i]] %>%
      filter(ALT_FREQ >= 0.03) %>%
      filter(TOTAL_DP >= 1000) %>%
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
names(saliva.intersecting) <- gsub("user","saliva",names(saliva.files))
saliva.intersecting <- lapply(saliva.intersecting,arrange,POS) 

for(i in seq_along(saliva.intersecting)){
  saliva.intersecting[[i]] <- saliva.intersecting[[i]] %>%
    mutate("SNP" = paste0(REF,POS,ALT))
}

combo.intersecting <- map2(nasal.intersecting,saliva.intersecting,bind_rows)
combo.intersecting <- lapply(combo.intersecting,distinct)

#pare down both datasets to include just intersecting iSNVs
#will include all instances of intersecting iSNVs 
#(even if they are outside of freq range on certain days)
intersect.only <- function(file, intersecting){
  for(i in seq_along(file)){
    file[[i]] <- file[[i]][which(file[[i]]$SNP %in% intersecting$SNP),]
  }
  return(file)
}

nasal.cut <- list()
for(i in seq_along(nasal.files)){
  nasal.cut[[i]] <- intersect.only(nasal.files.uncut[[i]], combo.intersecting[[i]])
}
names(nasal.cut) <- names(nasal.files)

saliva.cut <- list()
for(i in seq_along(saliva.files)){
  saliva.cut[[i]] <- intersect.only(saliva.files.uncut[[i]], combo.intersecting[[i]])
}
names(saliva.cut) <- names(saliva.files)

#label low dp as X
dp.relabel <- function(cut.file){
  for(i in seq_along(cut.file)){
    cut.file[[i]]$ALT_FREQ[which(cut.file[[i]]$TOTAL_DP < 1000)] <- "X"
  }
  return(cut.file)
}

nasal.cut <- lapply(nasal.cut, dp.relabel)
saliva.cut <- lapply(saliva.cut,dp.relabel)

#label alt_freq < 0.03 as 0 and alt_freq > 0.97 as 1
freq.relabel <- function(cut.file){
  for(i in seq_along(cut.file)){
    for(j in seq_along(cut.file[[i]]$POS)){
      if(cut.file[[i]]$ALT_FREQ[[j]] != "X"){
        if(cut.file[[i]]$ALT_FREQ[[j]] < 0.03){
          cut.file[[i]]$ALT_FREQ[[j]] <- 0
        }
        if(cut.file[[i]]$ALT_FREQ[[j]] > 0.97){
          cut.file[[i]]$ALT_FREQ[[j]] <- 1
        }
      }
    }
  }
  return(cut.file)
}

nasal.cut <- lapply(nasal.cut, freq.relabel)
saliva.cut <- lapply(saliva.cut, freq.relabel)

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

for(i in seq_along(nasal.cut)){
  comp.tables[[i]] <- env.comp(nasal.cut[[i]],saliva.cut[[i]])
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
    purrr::reduce(full_join, by=c("POS","REF","ALT","SNP")) 
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
names(comp.list) <- names(nasal.dirs)

#replace NA values (due to freq) with 0
#replace X values (due to dp) with NA
for(i in seq_along(comp.list)){
  comp.list[[i]] <- comp.list[[i]] %>%
    replace(is.na(.),0)
  comp.list[[i]][comp.list[[i]]=="X"] <- NA
}

#remove consensus variants 
remove.consensus <- function(user){
  info <- user[,1:4]
  freqs <- user[,grep("FREQ",colnames(user))]
  for(i in seq_along(rownames(freqs))){
    if(length(which(freqs[i,] < 1)) < 2){
      freqs[i,] <- "X"
    }
  }
  final <- cbind(info,freqs)
  final <- final %>%
    filter(NASAL_FREQ_1 != "X" | is.na(NASAL_FREQ_1))
  return(final)
}
comp.list <- lapply(comp.list,remove.consensus)

#write comparison table files
filenames <- paste0("comparison_tables/",names(comp.list),".csv")
for(i in seq_along(comp.list)){
  write.csv(comp.list[[i]],filenames[[i]])
}

save(comp.list, file = "comparison_list.RData")

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
unlist.comp <- na.omit(unlist.comp)
write.csv(unlist.comp,"nasal_saliva_freqs.csv")




