library(tidyverse)
library(rio)
library(trackViewer)
library(here)

#naive lolliplot
ten.palette <- c("#9D6A90","#8879A3","#6A88AE","#4795AE","#29A1A4","#30A991",
                 "#52AF79","#78B261","#A0B24E","#C8AF46")

genome <- GRanges(seqnames = "SARS-CoV-2",
                  ranges = IRanges(c(1,266,21563,25393,26245,27394,27894,28274,29558,29675),
                                   width = c(265,21290,3822,828,228,366,366,1260,117,229),
                                   names = c("5'UTR","ORF1ab","S","ORF3a","E","ORF7a","ORF8",
                                             "N","ORF10","3'UTR")),
                  fill = ten.palette,
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
         xaxis.gp = gpar(fontsize = 19),
         yaxis = F,
         cex=.8)

#vaccinated lolliplot
ten.palette <- c("#9D6A90","#8879A3","#6A88AE","#4795AE","#29A1A4","#30A991",
                 "#52AF79","#78B261","#A0B24E","#C8AF46")

genome <- GRanges(seqnames = "SARS-CoV-2",
                  ranges = IRanges(c(1,266,21563,25393,26245,27394,27894,28274,29558,29675),
                                   width = c(265,21290,3822,828,228,366,366,1260,117,229),
                                   names = c("5'UTR","ORF1ab","S","ORF3a","E","ORF7a","ORF8",
                                             "N","ORF10","3'UTR")),
                  fill = ten.palette,
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
          cex=.8)

#histogram of snp counts 
naive.counts <- naive_snp_table %>%
  select(COUNT)
naive.counts <- as.data.frame(table(naive.counts))
colnames(naive.counts) <- c("Count","Freq")
ggplot(naive.counts,aes(x=Count,y=Freq))+
  geom_col()+
  xlab("Number of Participants Sharing an iSNV")+
  ylab("Frequency")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19))
ggsave("figs/naive_variantcount.png")

vax.counts <- vax_snp_table %>%
  select(COUNT)
vax.counts <- as.data.frame(table(vax.counts))
colnames(vax.counts) <- c("Count","Freq")
ggplot(vax.counts,aes(x=Count,y=Freq))+
  geom_col()+
  xlab("Variant Count Across Participants")+
  ylab("Frequency")+
  theme_bw()+
  theme(axis.title = element_text(size = 22),axis.text = element_text(size = 19))
ggsave("figs/vax_variantcount.png")





