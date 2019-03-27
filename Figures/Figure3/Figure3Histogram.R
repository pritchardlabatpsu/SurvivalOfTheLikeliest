
setwd("C:/Users/Scott/Box Sync/MutBias/Figures/Figure3/")
library(ggplot2)
rm(list=ls())

imat.df = read.csv("Combined_data_frame_IC_Mutprob_abundance.csv")

imat.df$Compound = factor(imat.df$Compound,levels=imat.df$Compound[order(-imat.df$Abundance)])

ggplot(imat.df,aes(x=Compound,y=Abundance))+theme_bw()+
  geom_bar(aes(fill=Abundance),stat="identity")+
  xlab("Mutant")+ylab("Clinical Abundance")+
  scale_fill_gradient(low="#F8C391",high="#F39237")+
  theme(
    axis.title = element_text(size=20,face="bold"),
    axis.text = element_text(size=16,angle=90,color="black"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  guides(fill=F)

# ggsave("MutantHistogram.pdf",height=6,width=8)
