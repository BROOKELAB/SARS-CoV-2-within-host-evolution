library(tidyverse)
library(rio)
library(gtools)
library(Nmisc)
library(here)

#read saliva files
saliva.dirs <- list.dirs("saliva_output_tables")[-c(1)]
names(saliva.dirs) <- c("user_432686","user_433227","user_438577","user_442978",
                       "user_444332","user_444633","user_451152","user_451709",
                       "user_453058","user_459597","user_471467","user_471588")

saliva.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir,pattern="ivar"))
  files <- lapply(full, import)
  return(files)
}

saliva.files <- lapply(saliva.dirs, saliva.read)

seq.filter <- function(file){
  file <- file %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP) %>%
    filter(ALT_FREQ >= 0.03)%>%
    distinct()
  artifacts <- c(6696,11074,15965,29051,187,1059,2094,3037,
                 3130,6696,6990,8022,10323,10741,11074,13408,
                 14786,19684,20148,21137,24034,24378,25563,26144,
                 26461,26681,28077,28826,28854,29051,29700,29760)
  if(!identical(which(file$POS %in% artifacts), integer(0))){
    file <- file[-which(file$POS %in% artifacts),]
  }
  #not filtering out depth yet because we want to distinguish between
  #iSNVs not present and iSNVs present but at a low depth
}
for(i in seq_along(saliva.files)){
  saliva.files[[i]] <- lapply(saliva.files[[i]], seq.filter)
  saliva.files[[i]] <- lapply(saliva.files[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  saliva.files[[i]] <- lapply(saliva.files[[i]], relocate, SNP, .after = ALT)
}

#read nasal files
nasal.dirs <- list.dirs("nasal_output_tables_cutsaliva")[-c(1,2)]

names(nasal.dirs) <- names(saliva.dirs)

nasal.read <- function(dir){
  full <- paste0(dir,"/",list.files(dir, pattern = "ivar"))
  full <- mixedsort(sort(full)) #need this so it reads in date order
  files <- lapply(full, read.table,sep="\t",header =T)
  return(files)
}
nasal.files <- lapply(nasal.dirs, nasal.read)

seq.filter <- function(file){
  file <- file %>%
    select(POS,REF,ALT,ALT_FREQ,TOTAL_DP) %>%
    filter(ALT_FREQ >= 0.03)%>%
    distinct()
  artifacts <- c(6696,11074,15965,29051,187,1059,2094,3037,
                 3130,6696,6990,8022,10323,10741,11074,13408,
                 14786,19684,20148,21137,24034,24378,25563,26144,
                 26461,26681,28077,28826,28854,29051,29700,29760)
  if(!identical(which(file$POS %in% artifacts), integer(0))){
    file <- file[-which(file$POS %in% artifacts),]
  }
  #not filtering out depth yet because we want to distinguish between
  #iSNVs not present and iSNVs present but at a low depth
}
for(i in seq_along(nasal.files)){
  nasal.files[[i]] <- lapply(nasal.files[[i]], seq.filter)
  nasal.files[[i]] <- lapply(nasal.files[[i]], mutate,"SNP" = paste0(REF,POS,ALT))
  nasal.files[[i]] <- lapply(nasal.files[[i]], relocate, SNP, .after = ALT)
}

#name saliva samples
saliva.ids <- lapply(paste0("sampleIDs_saliva/",
                          list.files("sampleIDs_saliva/")), import,
                   header = F)
names(saliva.ids) <- names(saliva.files)
for(i in seq_along(saliva.files)){
  names(saliva.files[[i]]) <- saliva.ids[[i]]$V1
}

#name nasal samples
nasal.ids <- lapply(paste0("sampleIDs_nasalcutsaliva/",
                          list.files("sampleIDs_nasalcutsaliva/")), import,
                   header = F)
names(nasal.ids) <- names(nasal.files)
for(i in seq_along(nasal.files)){
  names(nasal.files[[i]]) <- nasal.ids[[i]]$V1
}

#create version of nasal.files that has same sample names as saliva.files
nasal.files.snames <- nasal.files
for(i in seq_along(nasal.files.snames)){
  names(nasal.files.snames[[i]]) <- saliva.ids[[i]]$V1
}

#remove saliva files with low mean coverage
load("naive_depth_files.RData")
naive.depth.files <- naive.depth.files[c(which(names(naive.depth.files) %in% 
                                                 names(saliva.files)))]

rm.vec <- function(user,depth){
  depth <- keep_at(depth, names(user))
  vec <- which(flatten(depth) < 1000)
  return(vec)
}

saliva.vec <- map2(saliva.files,naive.depth.files, rm.vec)

for(i in seq_along(saliva.vec)){
  if(length(saliva.vec[[i]]) != 0){
    saliva.files[[i]] <- saliva.files[[i]][-c(saliva.vec[[i]])]
  } else {
    saliva.files[[i]] <- saliva.files[[i]]
  }
}

#remove nasal files with low mean coverage
load("nasal_depth_files.RData")
nasal.depth.files <- nasal.depth.files[-c(7,14)]

rm.vec <- function(user,depth){
  depth <- keep_at(depth, names(user))
  vec <- which(flatten(depth) < 1000)
  return(vec)
}

nasal.vec <- map2(nasal.files,nasal.depth.files, rm.vec)

for(i in seq_along(nasal.vec)){
  if(length(nasal.vec[[i]]) != 0){
    nasal.files[[i]] <- nasal.files[[i]][-c(nasal.vec[[i]])]
    nasal.files.snames[[i]] <- nasal.files.snames[[i]][-c(nasal.vec[[i]])]
  } else {
    nasal.files[[i]] <- nasal.files[[i]]
    nasal.files.snames[[i]] <- nasal.files.snames[[i]]
  }
}

#remove saliva files with low Ct
saliva.ct <- import("naive_ct.xlsx")
saliva.ct <- saliva.ct %>%
  filter(user_id %in% gsub("user_","",names(saliva.files)))%>%
  arrange(user_id) %>%
  group_by(user_id, .add = TRUE) %>%
  group_split()

names(saliva.ct) <- names(saliva.files)
saliva.uncut <- saliva.files

rm.vec.ct <- function(user,ct){
  ct <- filter(ct, vdl_barcode %in% names(user))
  vec <- which(ct$ct >= 26)
  return(vec)
}

saliva.vec.ct <- map2(saliva.files,saliva.ct, rm.vec.ct)

for(i in seq_along(saliva.vec.ct)){
  if(length(saliva.vec.ct[[i]]) != 0){
    saliva.files[[i]] <- saliva.files[[i]][-c(saliva.vec.ct[[i]])]
  } else {
    saliva.files[[i]] <- saliva.files[[i]]
  }
}

#remove nasal files with low Ct
nasal.ct <- import("all_nasal_user_info.xlsx")
nasal.ct <- nasal.ct %>%
  filter(user_ID %in% gsub("user_","",names(nasal.files)))%>%
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
names(nasal.ct) <- names(nasal.files)
nasal.uncut <- nasal.files
nasal.uncut.snames <- nasal.files.snames

rm.vec.ct <- function(user,ct){
  ct <- filter(ct, sample %in% names(user))
  vec <- which(ct$cn >= 26)
  return(vec)
}

nasal.vec.ct <- map2(nasal.files,nasal.ct, rm.vec.ct)

for(i in seq_along(nasal.vec.ct)){
  if(length(nasal.vec.ct[[i]]) != 0){
    nasal.files[[i]] <- nasal.files[[i]][-c(nasal.vec.ct[[i]])]
    nasal.files.snames[[i]] <- nasal.files.snames[[i]][-c(nasal.vec.ct[[i]])]
  } else {
    nasal.files[[i]] <- nasal.files[[i]]
    nasal.files.snames[[i]] <- nasal.files.snames[[i]]
  }
}

#remove saliva files that do not exist in nasal dataset and vice-versa

rm.saliva <- function(saliva, nasal){
  s.keep <- which(names(saliva) %in% names(nasal))
  saliva <- saliva[s.keep]
  return(saliva)
}

saliva.files <- map2(saliva.files, nasal.files.snames, rm.saliva)
saliva.uncut <- map2(saliva.uncut, nasal.uncut.snames, rm.saliva)

rm.nasal <- function(nasal, saliva){
  n.keep <- which(names(nasal) %in% names(saliva))
  nasal <- nasal[n.keep]
  return(nasal)
}

nasal.files <- map2(nasal.files.snames, saliva.files, rm.nasal)
nasal.uncut <- map2(nasal.uncut.snames, saliva.uncut, rm.nasal)

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
  nasal.cut[[i]] <- intersect.only(nasal.uncut[[i]], combo.intersecting[[i]])
}
names(nasal.cut) <- names(nasal.files)

saliva.cut <- list()
for(i in seq_along(saliva.files)){
  saliva.cut[[i]] <- intersect.only(saliva.uncut[[i]], combo.intersecting[[i]])
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
save(unlist.comp, file = "nasal_saliva_freqs.RData")
write.csv(unlist.comp,"nasal_saliva_freqs.csv")

cor.test(as.numeric(unlist.comp$SALIVA), 
         as.numeric(unlist.comp$NASAL), method = "p")
#cor = 0.7025951 #p-value = 1.672e-14

