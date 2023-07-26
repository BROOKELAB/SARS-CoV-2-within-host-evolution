library(tidyverse)
library(here)

#for saliva
s.stderr <- function(ct,freq){
  copies <- (10 ^ (14.24 - (0.28*ct)))*0.01
  err <- sqrt((freq*(1-freq))/copies)
  return(err)
}

freq.seq <- seq(0,1,0.01)

s.ct22 <- map2_dbl(22,freq.seq,s.stderr)
s.ct24 <- map2_dbl(24,freq.seq,s.stderr)
s.ct26 <- map2_dbl(26,freq.seq,s.stderr)
s.ct28 <- map2_dbl(28,freq.seq,s.stderr)
s.ct30 <- map2_dbl(30,freq.seq,s.stderr)

s.ct.err <- as_tibble(cbind(freq.seq,s.ct22,s.ct24,s.ct26,s.ct28,s.ct30))
colnames(s.ct.err) <- c("Frequency", "22","24","26","28","30")
s.ct.err <- gather(s.ct.err,key = "Ct", value = "error", - Frequency)

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data = s.ct.err, aes(x = Frequency, y = error, color = Ct)) +
  geom_line(size = 0.9)+
  xlab("iSNV frequency")+
  ylab("Error")+
  ggtitle("Saliva samples")+
  scale_color_manual(values = five.palette)+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size =19),
        plot.title = element_text(size =22))
ggsave("figs/saliva_error.png")

#for nasal
n.stderr <- function(ct,freq){
  copies <- (10 ^ (11.35 - (0.25*ct)))*0.01
  err <- sqrt((freq*(1-freq))/copies)
  return(err)
}

freq.seq <- seq(0,1,0.01)

n.ct22 <- map2_dbl(22,freq.seq,n.stderr)
n.ct24 <- map2_dbl(24,freq.seq,n.stderr)
n.ct26 <- map2_dbl(26,freq.seq,n.stderr)
n.ct28 <- map2_dbl(28,freq.seq,n.stderr)
n.ct30 <- map2_dbl(30,freq.seq,n.stderr)

n.ct.err <- as_tibble(cbind(freq.seq,n.ct22,n.ct24,n.ct26,n.ct28,n.ct30))
colnames(n.ct.err) <- c("Frequency", "22","24","26","28","30")
n.ct.err <- gather(n.ct.err,key = "CN", value = "error", - Frequency)

five.palette <- c("#9D6A90","#618BAF","#27A59B","#6EB267","#C8AF46")
ggplot(data = n.ct.err, aes(x = Frequency, y = error, color = CN)) +
  geom_line(size = 0.9)+
  xlab("iSNV frequency")+
  ylab("Error")+
  ggtitle("Nasal samples")+
  scale_color_manual(values = five.palette)+
  theme_bw()+
  theme(axis.title = element_text(size = 22),
        axis.title.x = element_text(margin = margin(t = 12)),
        axis.title.y = element_text(margin = margin(r = 15)),
        axis.text = element_text(size = 19),
        legend.title = element_text(size = 22),
        legend.text = element_text(size =19),
        plot.title = element_text(size=22))
ggsave("figs/nasal_error.png")


