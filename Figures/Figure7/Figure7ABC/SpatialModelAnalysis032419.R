# Analyze spatial modeling results

rm(list=ls())
setwd("C:/Users/Scott/Box Sync/MutBias/SpatialModel/032419/")
library(ggplot2)
library(reshape)

# Import simulation results

res.df = read.csv("SpatialModelLoopResults032419.csv",header=F,na.strings="NaN")
colnames(res.df) = c("nviabS","nviabA","nviabB","ndead")
res.df = res.df[!is.na(res.df$nviabS),]
rownames(res.df) = 1:nrow(res.df)

res.df$ncells = rowSums(res.df)
res.df = res.df[res.df$ncells<500,] # Filter out patients that did not respond

res.df$propA = res.df$nviabA/rowSums(res.df[,c("nviabA","nviabB")])
res.df$propB = res.df$nviabB/rowSums(res.df[,c("nviabA","nviabB")])

fracAonly = sum(res.df$nviabB==0)
fracBonly = sum(res.df$nviabA==0)

# Organize data

res.df = res.df[order(res.df$propB),]

res.melt.df = data.frame(patient=1:nrow(res.df))
res.melt.df = cbind(res.melt.df,res.df[,c("propA","propB")])
res.melt.df = melt(res.melt.df,id="patient")

# Plot results
ggplot(res.melt.df,aes(x=patient,y=value))+theme_bw()+
  geom_bar(stat="identity",aes(fill=variable),width=1)+
  xlab("Patients")+
  ylab("Proportion of Tumor upon Relapse")+
  ggtitle("Population Restriction: Spatial Heterogeneity")+
  scale_fill_manual("Allele",label=c("A","B"),values=c("#2e368f","#ec2027"))+
  scale_x_continuous(breaks=c())+
  theme(
    plot.title=element_text(hjust=0.5,size=20,face="bold"),
    axis.text=element_text(size=16,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    legend.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16,face="bold")
  )

# ggsave("SpatialModelResults.pdf",width=8,height=6)

