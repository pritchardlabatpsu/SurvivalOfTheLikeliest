# Analyze infectious disease simulation results

library(ggplot2)

rm(list=ls())
setwd("C:/Users/Scott/Box Sync/MutBias/InfDisease/InfDisease030919/")

# Parse through files

sim.files = grep("InfDiseaseDelt",list.files(),value=T)

log.delt.vec = seq(-1,-7,length.out=100)
summary.df = data.frame(log.delt=log.delt.vec)

for (i in 1:length(sim.files)) {
  
  sim.file = sim.files[i]
  
  sim.df = read.csv(sim.file,header=F)
  pres.df = sim.df/rowSums(sim.df)
  
  pS = mean(pres.df$V1)
  pA = mean(pres.df$V2)
  pB = mean(pres.df$V3)
  
  sdS = sd(pres.df$V1)
  sdA = sd(pres.df$V2)
  sdB = sd(pres.df$V3)
  
  summary.df$pS[i] = pS
  summary.df$pA[i] = pA
  summary.df$pB[i] = pB
  
  summary.df$sdS[i] = sdS
  summary.df$sdA[i] = sdA
  summary.df$sdB[i] = sdB
  
}

ggplot(summary.df[summary.df$log.delt<=-2,],aes(x=log.delt))+theme_bw()+
  # geom_line(aes(y=pS,color="S"))+
  geom_line(aes(y=pA,color="A"),size=3)+
  geom_line(aes(y=pB,color="B"),size=3)+
  scale_color_manual("Allele",values=c("#2e368f","#ec2027"))+
  scale_x_continuous(breaks=c(-2:-7),labels=parse(text=c("10^-2","10^-3","10^-4","10^-5","10^-6","10^-7")))+
  ggtitle("Population Restriction: Transmission Bottleneck")+
  xlab("Fractional Bottleneck")+
  ylab("Allele Prevalence")+
  theme(
    plot.title=element_text(hjust=0.5,size=20,face="bold"),
    axis.text=element_text(size=16,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    legend.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16,face="bold")
  )

# ggsave("BottleneckResults.pdf",width=8,height=6)
