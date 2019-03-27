# Analyze adjuvant therapy simulation results

library(MASS)
library(ggplot2)
library(DescTools)

rm(list=ls())
setwd("C:/Users/Scott/Box Sync/MutBias/AdjuvantTherapy/AT022719/")

## Generated parameters

pars.df = read.csv("ATParameters022719.csv",header=F)
colnames(pars.df) = c("mut.bias","rel.alpha")

Rsites = nrow(pars.df)
vars = as.character(0:Rsites)

## Parse through simulation results

all.files = list.files()
sim.files = grep("ATSimsT[0-9]{2}R[0-9]{2}_022719.csv",all.files,value=T)

T.vec = as.numeric(unique(substr(sim.files,8,9)))
R.vec = as.numeric(unique(substr(sim.files,11,12)))

df1 = expand.grid(log.remain=R.vec,log.treat = T.vec)
df1 = df1[,c(2,1)]
df2 = as.data.frame(matrix(nrow=nrow(df1),ncol=3))
colnames(df2) = c("mut.bias.pval","rel.alpha.pval","pval.ratio")
summary.df = cbind(df1,df2)

for (i in 1:length(sim.files)) {
  
  sim.file = sim.files[i]
  
  log.treat = as.numeric(unique(substr(sim.file,8,9)))
  log.remain = as.numeric(unique(substr(sim.file,11,12)))

  sims.df = read.csv(sim.file,header=F)
  colnames(sims.df) = "dom.allele"
  
  prev = table(sims.df$dom.allele)
  prev = prev[vars]
  names(prev) = vars
  prev[is.na(prev)] = 0
  
  simsum.df = cbind(
    pars.df,
    data.frame(unname(prev[-1]))
  )
  
  
  # fit.nbin = glm.nb(Freq~mut.bias+rel.alpha,data=simsum.df)
  # fit.sum = summary(fit.nbin)
  # mut.bias.pval = fit.sum$coefficients["mut.bias","Pr(>|z|)"]
  # rel.alpha.pval = fit.sum$coefficients["rel.alpha","Pr(>|z|)"]
  # 
  # GKG.bias = GoodmanKruskalGamma(x=simsum.df$mut.bias,y=simsum.df$Freq)
  # GKG.negalpha = GoodmanKruskalGamma(x=-simsum.df$rel.alpha,y=simsum.df$Freq)
  
  cor.bias = cor(simsum.df$mut.bias,simsum.df$Freq)
  cor.negalpha = cor(-simsum.df$rel.alpha,simsum.df$Freq)
  
  spear.bias = cor(simsum.df$mut.bias,simsum.df$Freq,method="spearman")
  spear.negalpha = cor(-simsum.df$rel.alpha,simsum.df$Freq,method="spearman")
  
  summary.row = summary.df$log.treat==log.treat & summary.df$log.remain==log.remain
  
  # summary.df$mut.bias.pval[summary.row] = mut.bias.pval
  # summary.df$rel.alpha.pval[summary.row] = rel.alpha.pval
  # summary.df$pval.ratio[summary.row] = mut.bias.pval/rel.alpha.pval
  # 
  # summary.df$GKG.bias[summary.row] = GKG.bias
  # summary.df$GKG.negalpha[summary.row] = GKG.negalpha
  
  summary.df$cor.bias[summary.row] = cor.bias
  summary.df$cor.negalpha[summary.row] = cor.negalpha
  
  summary.df$spear.bias[summary.row] = spear.bias
  summary.df$spear.negalpha[summary.row] = spear.negalpha
  
  if ((log.treat==9&log.remain==5)) {

    max.freq = round(max(simsum.df$Freq/sum(simsum.df$Freq)),2)+0.01

    ggplot(simsum.df,aes(x=-rel.alpha+max(rel.alpha),y=Freq/sum(Freq)))+theme_bw()+
      geom_point(size=7.5,color="#2e368f")+
      geom_rect(xmin=0.0,xmax=0.06,ymin=0,ymax=0.04,alpha=0.045,fill="#2e368f")+
      annotate("text",x=0.03,y=0.02,size=12,parse=T,
               label=paste("rho:",as.character(round(spear.negalpha,2)),sep=""))+
      # ggtitle(paste("Spearman Correlation:",as.character(round(spear.negalpha,2)),sep=" "))+
      scale_x_continuous("Resistance",breaks=c(0,0.25),limits=c(0,0.25))+
      scale_y_continuous("Frequency",breaks=c(0,max.freq),limits=c(0,max.freq))+
      theme(
        plot.title=element_text(hjust=0.5,size=20,face="bold",color="#2e368f"),
        axis.text=element_text(size=20,face="bold",color="black"),
        axis.title=element_text(size=30,face="bold")
      )
    
    ggsave(paste("CorResT",as.character(log.treat),"R",as.character(log.remain),".pdf",sep=""),
           width=8,height=4)

    ggplot(simsum.df,aes(x=mut.bias,y=Freq/sum(Freq)))+theme_bw()+
      geom_point(size=7.5,color="#ec2027")+
      geom_rect(xmin=0,xmax=0.05,ymin=0.14,ymax=0.18,alpha=0.045,fill="#ec2027")+
      annotate("text",x=0.025,y=0.16,size=12,parse=T,
               label=paste("rho:",as.character(round(spear.bias,2)),sep=""))+
      # ggtitle(paste("Spearman Correlation:",as.character(round(spear.bias,2)),sep=" "))+
      scale_x_continuous("Mutational Probability",breaks=c(0,0.2),limits=c(0,0.2))+
      scale_y_continuous("Frequency",breaks=c(0,max.freq),limits=c(0,max.freq))+
      theme(
        plot.title=element_text(hjust=0.5,size=20,face="bold",color="#ec2027"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=30,face="bold")
      )
    
    ggsave(paste("CorBiasT",as.character(log.treat),"R",as.character(log.remain),".pdf",sep=""),
           width=8,height=4)

  }
  
  if ((log.treat==12&log.remain==8)) {
    
    max.freq = round(max(simsum.df$Freq/sum(simsum.df$Freq)),2)+0.01
    
    ggplot(simsum.df,aes(x=-rel.alpha+max(rel.alpha),y=Freq/sum(Freq)))+theme_bw()+
      geom_point(size=7.5,color="#2e368f")+
      geom_rect(xmin=0.0,xmax=0.06,ymin=0.38,ymax=0.48,alpha=0.045,fill="#2e368f")+
      annotate("text",x=0.03,y=0.43,size=12,parse=T,
               label=paste("rho:",as.character(round(spear.negalpha,2)),sep=""))+
      # ggtitle(paste("Spearman Correlation:",as.character(round(spear.negalpha,2)),sep=" "))+
      scale_x_continuous("Resistance",breaks=c(0,0.25),limits=c(0,0.25))+
      scale_y_continuous("Frequency",breaks=c(0,max.freq),limits=c(0,max.freq))+
      theme(
        plot.title=element_text(hjust=0.5,size=20,face="bold",color="#2e368f"),
        axis.text=element_text(size=20,face="bold",color="black"),
        axis.title=element_text(size=30,face="bold")
      )
    
    ggsave(paste("CorResT",as.character(log.treat),"R",as.character(log.remain),".pdf",sep=""),
           width=8,height=4)
    
    ggplot(simsum.df,aes(x=mut.bias,y=Freq/sum(Freq)))+theme_bw()+
      geom_point(size=7.5,color="#ec2027")+
      geom_rect(xmin=0.0,xmax=0.05,ymin=0.38,ymax=0.48,alpha=0.045,fill="#ec2027")+
      annotate("text",x=0.025,y=0.43,size=12,parse=T,
               label=paste("rho:",as.character(round(spear.bias,2)),sep=""))+
      # ggtitle(paste("Spearman Correlation:",as.character(round(spear.bias,2)),sep=" "))+
      scale_x_continuous("Mutational Probability",breaks=c(0,0.2),limits=c(0,0.2))+
      scale_y_continuous("Frequency",breaks=c(0,max.freq),limits=c(0,max.freq))+
      theme(
        plot.title=element_text(hjust=0.5,size=20,face="bold",color="#ec2027"),
        axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=30,face="bold")
      )
    
    ggsave(paste("CorBiasT",as.character(log.treat),"R",as.character(log.remain),".pdf",sep=""),
           width=8,height=4)
    
  }

}

# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=log10(mut.bias.pval)))
# 
# ggplot(summary.df[summary.df$log.treat!=9|summary.df$log.remain!=5,],aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=log10(rel.alpha.pval)))
# 
# ggplot(summary.df[summary.df$log.treat!=9|summary.df$log.remain!=5,],aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=log10(pval.ratio)))
# 
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=GKG.bias))+
#   scale_fill_continuous(limits=c(-.1,1),low="pink",high="dark blue")
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=GKG.negalpha))+
#   scale_fill_continuous(limits=c(-.1,1),low="pink",high="dark blue")
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=GKG.bias-GKG.negalpha))
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=cor.bias))
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=cor.negalpha))

ggplot(summary.df[summary.df$log.remain!=4,],aes(x=log.treat,y=log.remain))+theme_bw()+
  geom_raster(aes(fill=spear.bias-spear.negalpha))+
  xlab("Population Size Before Resection")+
  ylab("Population Size After Resection")+
  ggtitle("Population Restriction: Adjuvant Therapy")+
  scale_x_continuous(breaks=9:12,labels=parse(text=c("10^9","10^10","10^11","10^12")))+
  scale_y_continuous(breaks=4:8,labels=parse(text=c("10^4","10^5","10^6","10^7","10^8")))+
  scale_fill_continuous("Difference in\nCorrelation",
                        low="#2e368f",high="#ec2027")+
  theme(
    plot.title=element_text(hjust=0.5,size=20,face="bold"),
    axis.text=element_text(size=16,face="bold"),
    axis.title=element_text(size=16,face="bold"),
    legend.title=element_text(size=18,face="bold"),
    legend.text=element_text(size=16,face="bold")
  )

# ggsave("ATSimResults.pdf",width=8,height=6)
  

# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=spear.bias))
# 
# ggplot(summary.df,aes(x=log.treat,y=log.remain))+theme_bw()+
#   geom_raster(aes(fill=spear.negalpha))