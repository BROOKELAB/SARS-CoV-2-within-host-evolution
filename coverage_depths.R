library(tidyverse)
library(rio)
library(rlist)
library(gtools)
library(here)

#saliva
cov.names <- list.files("saliva_coverage_summaries/")
cov.files <- paste0("saliva_coverage_summaries/",cov.names)
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
  geom_histogram(binwidth = 3000)+
  xlab("Mean coverage")+
  ylab("Count")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22))
ggsave("figs/saliva_coverage.png")

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
  geom_histogram(binwidth = 150)+
  xlab("Median coverage")+
  ylab("Count")+
  theme_bw()+
  theme(axis.text = element_text(size = 19),
        axis.title = element_text(size = 22))
ggsave("figs/nasal_coverage.png")
