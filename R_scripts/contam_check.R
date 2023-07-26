library(tidyverse)
library(rio)
library(rlist)
library(Nmisc)
library(gdata)
library(gtools)
library(here)

#combine all sample info
naive.dirlist <- list.dirs("naive_output_tables")[-c(1:3)]
vax.dirlist <- list.dirs("vaccinated_output_tables")[-1]
dirlist <- c(naive.dirlist, vax.dirlist)
snp.filter <- function(dir){
  files <- list.files(dir,pattern = "ivar")
  files <- gtools::mixedsort(sort(files))
  full <- paste0(dir,"/",files)
  user <- list()
  for(i in 1:length(full)){
    user[[i]] <- import(full[i])
  }
  sleekuser <- list()
  for(i in 1:length(user)){
    sleekuser[[i]] <- user[[i]] %>%
      filter(ALT_FREQ>=.03)%>%
      filter(TOTAL_DP >=1000)%>%
      select(POS,REF,ALT,ALT_DP,TOTAL_DP,ALT_FREQ)%>%
      distinct()
  }
  return(sleekuser)
} 
users <- lapply(dirlist,snp.filter)
names(users) <- c("432686","432870","433227","435786","435805","438577",
                      "442978","444332","444446","444633","445602","449650",
                      "450241","450348","451152","451709","453058","459597",
                      "471467","471588","461913","471876","475670","481242",
                      "481672","482828","484249","486422","487250","487297",
                      "487941")

#name samples
load("naive_depth_files.RData")
load("vax_depth_files.RData")
depth.files <- c(naive.depth.files, vax.depth.files)
name.files <- function(user, depth){
  names(user) <- names(depth)
  return(user)
}
users <- map2(users,depth.files, name.files)

#pull all mutations from samples of interest
#2301664 (user 432870)
#2309019 (user 433227)
#2322121 (user 433227)
#2453549 (user 444332)

sample.list <- list(users$`432870`$`2301664`, users$`433227`$`2309019`,
                    users$`433227`$`2322121`, users$`444332`$`2453549`)
names(sample.list) <- c("2301664","2309019","2322121","2453549")
for(i in seq_along(sample.list)){
  sample.list[[i]] <- mutate(sample.list[[i]], "SAMPLE_ID" = as.numeric(names(sample.list)[[i]]))
}

#pull all mutations from all samples 
all <- flatten(users)
for(i in seq_along(all)){
  all[[i]] <- mutate(all[[i]], "SAMPLE_ID" = as.numeric(names(all)[[i]]))
}
all <- bind_rows(all)
all <- all %>%
  arrange(POS) %>%
  mutate("SNP" = paste0(REF,POS,ALT))%>%
  select(SNP, ALT_FREQ, SAMPLE_ID)


#check for unique SNPs
contam.check <- function(sample) {
  all <- all %>%
    select(SNP, ALT_FREQ)
  sample <- sample %>%
    mutate("SNP" = paste0(REF,POS,ALT)) %>%
    select(SNP, ALT_FREQ)
  compare <- left_join(sample, all, by = "SNP")
  colnames(compare)[c(2:3)] <- c("sample_freq", "everything_freq")
  counts <- table(compare$SNP)
  onlys <- counts[counts == 1]
  compare <- compare[which(compare$SNP %in% names(onlys)),]
  compare <- select(compare, SNP, sample_freq)
  return(compare)
}

checked <- lapply(sample.list,contam.check)

snp.counts <- list()
snp.min <- list()
snp.max <- list()
for(i in seq_along(checked)){
  snp.counts[[i]] <- length(checked[[i]]$sample)
  snp.min[[i]] <- min(checked[[i]]$sample)
  snp.max[[i]] <- max(checked[[i]]$sample)
}

check.table <- as.data.frame(matrix(data = NA, nrow = 4, ncol = 5))
colnames(check.table) <- c("Participant","Sample","Unique SNPs","Min. SNP frequency","Max. SNP frequency")
check.table[,1] <- c("432870","433227","433227","444332")
check.table[,2] <- names(checked)
check.table[,3] <- unlist(snp.counts)
check.table[,4] <- unlist(snp.min)
check.table[,5] <- unlist(snp.max)

write.csv(check.table, file = "contam_check_table.csv")

#check the correlation between sample of interest and every other sample in the dataset
all.list <- flatten(users)
all.list <- lapply(all.list, mutate, "SNP" = paste0(REF,POS,ALT))
all.list <- lapply(all.list, select, SNP, ALT_FREQ)

check.cor <- function(sample){
  sample <- mutate(sample, "SNP" = paste0(REF,POS,ALT))
  sample <- select(sample, SNP, ALT_FREQ)
  
  tables <- list()
  for(i in seq_along(all.list)){
    tables[[i]] <- full_join(sample, all.list[[i]], by = "SNP")
    colnames(tables[[i]])[c(2:3)] <- c("SOI", "test_sample")
    tables[[i]][is.na(tables[[i]])] <- 0
  }
  names(tables) <- names(all.list)
  
  cor.list <- list()
  p.list <- list()
  for(i in seq_along(tables)){
    cor.list[[i]] <- suppressWarnings(cor.test(tables[[i]]$SOI, tables[[i]]$test_sample, method = "p")$estimate)
    p.list[[i]] <- suppressWarnings(cor.test(tables[[i]]$SOI, tables[[i]]$test_sample, method = "p")$p.value)
  }
  cor.list <- unlist(cor.list)
  cor.list[which(is.na(cor.list))] <- 0
  p.list <- unlist(p.list)
  p.list[which(is.na(p.list))] <- 0
  cor.table <- as.data.frame(matrix(data = NA, nrow = length(cor.list),ncol = 3))
  colnames(cor.table) <- c("sample","cor","p_val")
  cor.table$sample <- names(tables)
  cor.table$cor <- cor.list
  cor.table$p_val <- p.list
  return(cor.table)
}

cor.tables <- lapply(sample.list, check.cor)

a <- 0.05/length(cor.tables[[1]]$sample)
  
for(i in seq_along(cor.tables)){
  cor.tables[[i]] <- mutate(cor.tables[[i]], "SIG" = NA)
  cor.tables[[i]]$SIG[which(cor.tables[[i]]$p_val < a)] <- "P < \u03b1"
  cor.tables[[i]]$SIG[which(cor.tables[[i]]$p_val >= a)] <- "P >= \u03b1"
}

#which comparisons are significant
load("run_cors.RData")
info <- import("all_saliva_user_info.xlsx")

match.sig <- function(cortable){
  high <- cortable[which(cortable$cor > min(run.cors)),]
  sig <- high[which(high$p_val < a),]
  sig <- sig %>%
    mutate("participant" = NA)%>%
    relocate(participant, .before = sample)%>%
    arrange(desc(cor))
  for(i in seq_along(sig$sample)){
    sig$participant[[i]] <- info$user_id[which(info$sample_barcode == sig$sample[[i]])]
  }
  return(sig)
}

sig.tables <- cor.tables
sig.tables <- lapply(sig.tables, match.sig)


#plot
load("run_cors.RData")

names(cor.tables) <- c("432870 - Day 7", "433227 - Day 5", "433227 - Day 7",
                       "444332 - Day 7")
names(sig.tables) <- names(cor.tables)
cor.plots <- list()
for(i in seq_along(cor.tables)){
  cor.plots[[i]] <- ggplot(data = cor.tables[[i]], 
                           aes(x = 1:length(sample), 
                               y = cor,
                               col = SIG))+
    geom_point()+
    geom_hline(yintercept = min(run.cors), linetype = 2)+
    ggtitle(names(cor.tables)[[i]])+
    ylab("Correlation")+
    scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
    scale_color_manual(values = c("red","black"))+
    theme_bw()+
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = 22, margin = margin(r=15)),
          axis.text.y = element_text(size = 19),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 19),
          plot.title = element_text(size = 22))+
    guides(color = guide_legend(title = "Significance"))
}




