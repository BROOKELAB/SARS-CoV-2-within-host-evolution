library(tidyverse)
library(rio)
library(trackViewer)
library(here)

#naive lolliplot
thirteen.palette <- c("#9D6A90","#8E759F","#7981AA","#618CAF","#4696AE","#2E9EA7",
                      "#27A59B","#37AB8B","#52AF79","#6FB267","#8CB356","#AAB24A",
                      "#C8AF46")

genome <- GRanges(seqnames = "SARS-CoV-2",
                  ranges = IRanges(c(1,266,21563,25393,26245,26523,27202,
                                     27394,27756, 27894,28274,29558,29675),
                                   width = c(265,21290,3822,828,228,669,186,
                                             366,132,366,1260,117,229),
                                   names = c("5'UTR","ORF1ab","S","ORF3a","E","M","ORF6",
                                             "ORF7a","ORF7b","ORF8","N","ORF10","3'UTR")),
                  fill = thirteen.palette,
                  height = .05)
naive_annotations <- import("naive_annotations.csv")
naive_annotations <- naive_annotations %>%
  select(-c(V1,ID,AA))%>%
  arrange(POS)
naive_counts <- table(factor(naive_annotations$SNP, levels=unique(naive_annotations$SNP)))
naive_annotations <- distinct(naive_annotations)
naive_annotations <- cbind(naive_annotations,naive_counts)
naive_snp_table <- naive_annotations %>%
  select(-Var1)%>%
  dplyr::rename("COUNT" = Freq)%>%
  mutate("color" = NA)

naiveSNPs <- as.vector(naive_snp_table$POS)

for(i in 1:length(naive_snp_table$EFF)){
  if(naive_snp_table$EFF[[i]]=="synonymous"){
    naive_snp_table$color[[i]] <- "grey80"
  }
  if(naive_snp_table$EFF[[i]]=="nonsynonymous"){
    naive_snp_table$color[[i]] <- "grey30"
  }
}

SNPranges <- GRanges("SARS-CoV-2",
                     IRanges(naiveSNPs,width = 1,names = NULL),
                     score = naive_snp_table$COUNT,
                     color = naive_snp_table$color)
lolliplot(SNPranges,genome,
          xaxis = c(0,5000,10000,15000,20000,25000,30000),
          yaxis = F,
          xaxis.gp = gpar(fontsize = 19),
          cex=.5)

#vaccinated lolliplot
thirteen.palette <- c("#9D6A90","#8E759F","#7981AA","#618CAF","#4696AE","#2E9EA7",
                      "#27A59B","#37AB8B","#52AF79","#6FB267","#8CB356","#AAB24A",
                      "#C8AF46")
genome <- GRanges(seqnames = "SARS-CoV-2",
                  ranges = IRanges(c(1,266,21563,25393,26245,26523,27202,
                                     27394,27756, 27894,28274,29558,29675),
                                   width = c(265,21290,3822,828,228,669,186,
                                             366,132,366,1260,117,229),
                                   names = c("5'UTR","ORF1ab","S","ORF3a","E","M","ORF6",
                                             "ORF7a","ORF7b","ORF8","N","ORF10","3'UTR")),
                  fill = thirteen.palette,
                  height = .05)

vax_annotations <- import("vax_annotations.csv")
vax_annotations <- vax_annotations %>%
  select(-c(V1,ID,AA))%>%
  arrange(POS)
vax_counts <- table(factor(vax_annotations$SNP, levels=unique(vax_annotations$SNP)))
vax_annotations <- distinct(vax_annotations)
vax_annotations <- cbind(vax_annotations,vax_counts)
vax_snp_table <- vax_annotations %>%
  select(-Var1)%>%
  dplyr::rename("COUNT" = Freq)%>%
  mutate("color" = NA)

vaxSNPs <- as.vector(vax_snp_table$POS)
for(i in 1:length(vax_snp_table$EFF)){
  if(vax_snp_table$EFF[[i]]=="synonymous"){
    vax_snp_table$color[[i]] <- "grey80"
  }
  if(vax_snp_table$EFF[[i]]=="nonsynonymous"){
    vax_snp_table$color[[i]] <- "grey30"
  }
}

SNPranges <- GRanges("SARS-CoV-2",
                     IRanges(vaxSNPs,width = 1,names = NULL),
                     score = vax_snp_table$COUNT,
                     color = vax_snp_table$color)
lolliplot(SNPranges,genome,
          xaxis = c(0,5000,10000,15000,20000,25000,30000),
          yaxis = F,
          xaxis.gp = gpar(fontsize = 19),
          cex=.5)

#dispersion calculations
snp.window <- function(pos,ann){
  range <- pos:(pos+100)
  snps <- ann[which(ann$POS %in% range),]
  snp.count <- sum(snps$Freq)
  snp.table <- as.data.frame(cbind(range[1],range[100],snp.count))
  colnames(snp.table) <- c("Start","End","Count")
  return(snp.table)
}

snp.density.naive <- list()
for(i in seq_along(1:29803)){
  snp.density.naive[[i]] <- snp.window(i,naive_annotations)
}
snp.density.naive <- bind_rows(snp.density.naive)
snp.mean.naive <- mean(snp.density.naive$Count)

poisson.snp <- function(region,mean){
  e <- 2.718
  x <- region
  p <- ((e^-mean) * (mean^x))/ (factorial(x))
  return(p)
}

snp.density.naive <- snp.density.naive %>%
  mutate("P" = poisson.snp(Count, snp.mean.naive))

#bonferroni correction
a <- 0.05
obs <- 29803
a.a <- a/obs

which(snp.density.naive$P < a.a) #none

ggplot(data = snp.density.naive, aes(x = Start, y = log(P)))+
  geom_segment(data = snp.density.naive, aes(x = Start, xend = End,
                                      y = log(P), yend = log(P),
                                      color = log(P)),
               size = 4)+
  xlab("Nucleotide position")+
  ggtitle("Unvaccinated")+
  ylim(c(-30,0))+
  scale_color_viridis_c(limits = c(-30,0))+
  geom_hline(yintercept = log(a.a), linetype = 2)+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t=10)),
        axis.title.y = element_text(margin = margin(r=12)),
        axis.text = element_text(size = 19),
        plot.title = element_text(size=22),
        legend.position = "none")
ggsave("figs/naive_poisson.png")

#vaccinated
snp.density.vax <- list()
for(i in seq_along(1:29803)){
  snp.density.vax[[i]] <- snp.window(i,vax_annotations)
}
snp.density.vax <- bind_rows(snp.density.vax)
snp.mean.vax <- mean(snp.density.vax$Count)

snp.density.vax <- snp.density.vax %>%
  mutate("P" = poisson.snp(Count, snp.mean.vax))

which(snp.density.vax$P < a.a)

ggplot(data = snp.density.vax, aes(x = Start, y = log(P)))+
  geom_segment(data = snp.density.vax, aes(x = Start, xend = End,
                                             y = log(P), yend = log(P),
                                             color = log(P)),
               size = 4)+
  xlab("Nucleotide position")+
  ggtitle("Vaccinated")+
  scale_y_continuous(limits = c(-30, 0), breaks = c(-30,-20,-10,0))+
  scale_color_viridis_c(limits = c(-30,0))+
  geom_hline(yintercept = log(a.a), linetype = 2)+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t=10)),
        axis.title.y = element_text(margin = margin(r=12)),
        axis.text = element_text(size = 19),
        plot.title = element_text(size=22),
        legend.position = "none")
ggsave("figs/vax_poisson.png")









