library(tidyverse)
library(rio)
library(here)

##unvaccinated##
#read in genomic coverage info for each sample
naive.depth.dirs <- list.dirs("naive_depths")[-1]
naive.depth.files <- lapply(naive.depth.dirs, list.files)
naive.depth.names <- list()
for(i in seq_along(naive.depth.files)){
  naive.depth.names[[i]] <- gsub(".mosdepth.summary.txt","",naive.depth.files[[i]])
  naive.depth.names[[i]] <- gsub(".trim.genome","",naive.depth.names[[i]])
}

for(i in seq_along(naive.depth.dirs)){
  naive.depth.files[[i]] <- lapply(paste0(naive.depth.dirs[[i]],"/",naive.depth.files[[i]]), 
                                   import)
  names(naive.depth.files[[i]]) <- naive.depth.names[[i]]
}
names(naive.depth.files) <- gsub("naive_depths/","",naive.depth.dirs)

#pull out mean coverage values
for(i in seq_along(naive.depth.files)){
  naive.depth.files[[i]] <- lapply(naive.depth.files[[i]], select, mean)
  naive.depth.files[[i]] <- lapply(naive.depth.files[[i]], distinct)
}
save(naive.depth.files, file = "naive_depth_files.RData")


#per-participant genomic coverage tables
naive.depth.tables <- lapply(naive.depth.files, bind_rows)
naive.depth.tables <- lapply(naive.depth.tables, as_tibble)
for(i in seq_along(naive.depth.tables)){
  naive.depth.tables[[i]] <- naive.depth.tables[[i]] %>%
    mutate("sample" = names(naive.depth.files[[i]])) %>%
    mutate("participant" = names(naive.depth.files)[[i]])%>%
    relocate(sample, .before = mean) %>%
    relocate(participant, .before = sample)%>%
    dplyr::rename("mean_coverage" = mean)
}

naive.depth.table <- bind_rows(naive.depth.tables)
write.csv(naive.depth.table, "naive_depth_table.csv")

##vaccinated##
#read in genomic coverage info for each sample
vax.depth.dirs <- list.dirs("vax_depths")[-1]
vax.depth.files <- lapply(vax.depth.dirs, list.files)
vax.depth.names <- list()
for(i in seq_along(vax.depth.files)){
  vax.depth.names[[i]] <- gsub(".mosdepth.summary.txt","",vax.depth.files[[i]])
  vax.depth.names[[i]] <- gsub(".trim.genome","",vax.depth.names[[i]])
}

for(i in seq_along(vax.depth.dirs)){
  vax.depth.files[[i]] <- lapply(paste0(vax.depth.dirs[[i]],"/",vax.depth.files[[i]]), 
                                   import)
  names(vax.depth.files[[i]]) <- vax.depth.names[[i]]
}
names(vax.depth.files) <- gsub("vax_depths/","",vax.depth.dirs)

#pull out mean coverage values
for(i in seq_along(vax.depth.files)){
  vax.depth.files[[i]] <- lapply(vax.depth.files[[i]], select, mean)
  vax.depth.files[[i]] <- lapply(vax.depth.files[[i]], distinct)
}
save(vax.depth.files, file = "vax_depth_files.RData")

#per-participant genomic coverage tables
vax.depth.tables <- lapply(vax.depth.files, bind_rows)
vax.depth.tables <- lapply(vax.depth.tables, as_tibble)
for(i in seq_along(vax.depth.tables)){
  vax.depth.tables[[i]] <- vax.depth.tables[[i]] %>%
    mutate("sample" = names(vax.depth.files[[i]])) %>%
    mutate("participant" = names(vax.depth.files)[[i]])%>%
    relocate(sample, .before = mean) %>%
    relocate(participant, .before = sample)%>%
    dplyr::rename("mean_coverage" = mean)
}

vax.depth.table <- bind_rows(vax.depth.tables)
write.csv(vax.depth.table, "vax_depth_table.csv")

##nasal##
#read in genomic coverage info for each sample
nasal.depth.dirs <- list.dirs("nasal_depths")[-1]
nasal.depth.files <- lapply(nasal.depth.dirs, list.files)
nasal.depth.names <- list()
for(i in seq_along(nasal.depth.files)){
  nasal.depth.names[[i]] <- gsub(".mosdepth.summary.txt","",nasal.depth.files[[i]])
  nasal.depth.names[[i]] <- gsub(".trim.genome","",nasal.depth.names[[i]])
}

for(i in seq_along(nasal.depth.dirs)){
  nasal.depth.files[[i]] <- lapply(paste0(nasal.depth.dirs[[i]],"/",nasal.depth.files[[i]]), 
                                   import)
  names(nasal.depth.files[[i]]) <- nasal.depth.names[[i]]
}
names(nasal.depth.files) <- gsub("nasal_depths/","",nasal.depth.dirs)

#pull out mean coverage values
for(i in seq_along(nasal.depth.files)){
  nasal.depth.files[[i]] <- lapply(nasal.depth.files[[i]], select, mean)
  nasal.depth.files[[i]] <- lapply(nasal.depth.files[[i]], distinct)
}
save(nasal.depth.files, file = "nasal_depth_files.RData")


#per-participant genomic coverage tables
nasal.depth.tables <- lapply(nasal.depth.files, bind_rows)
nasal.depth.tables <- lapply(nasal.depth.tables, as_tibble)
for(i in seq_along(nasal.depth.tables)){
  nasal.depth.tables[[i]] <- nasal.depth.tables[[i]] %>%
    mutate("sample" = names(nasal.depth.files[[i]])) %>%
    mutate("participant" = names(nasal.depth.files)[[i]])%>%
    relocate(sample, .before = mean) %>%
    relocate(participant, .before = sample)%>%
    dplyr::rename("mean_coverage" = mean)
}

nasal.depth.table <- bind_rows(nasal.depth.tables)
write.csv(nasal.depth.table, "nasal_depth_table.csv")



