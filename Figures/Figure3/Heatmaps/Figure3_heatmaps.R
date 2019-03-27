
setwd("C:/Users/Scott/Box Sync/MutBias/R/")


##### IC50 Heat Map

rm(list=ls())

IC50df = read.csv('IC50HeatMap.csv')

doses = 10000*2.^-(0:9) # nM

pIC50 = data.frame()

for (i in 1:nrow(IC50df)) {
  
  pIC50 = rbind(
    pIC50,
    data.frame(genotype=IC50df[i,1],n=i,dose=doses,val=as.numeric(IC50df[i,13:22]))
  )
  
}

# ggplot(pIC50,aes(x=log10(dose),y=-n))+theme_bw()+
#   geom_raster(aes(fill=val))+
#   xlab('Dose [nM]')+ylab('Resistance Variant')+
#   ggtitle('Dose Response for BCR-ABL Variants')+
#   scale_x_continuous(breaks=c(1,2,3,4),labels=parse(text=c('10^1','10^2','10^3','10^4')))+
#   scale_y_continuous(breaks=c())+
#   scale_fill_continuous(name='Viability',breaks=c(0,1))+
#   theme(
#     plot.title=element_text(hjust=0.5,size=30,face='bold'),
#     axis.title=element_text(size=26,face='bold'),
#     axis.text=element_text(size=20,color='black'),
#     legend.title=element_text(size=20),
#     legend.text=element_text(size=15)
#   )


# alt: sort by IC50

x=doses

for (i in 1:nrow(IC50df)) {
  
  IC50df[i,"noriginal"] = i
  
  y=as.numeric(IC50df[i,13:22])
  
  # nlm has trouble with some, so eyeball IC50
  
  # if (i == 29) {
  #   IC50df[i,"IC50"] = 3.4
  # }
  if (i == 108) {
    IC50df[i,"IC50"] = 2.8
  }
  else if (i==172) {
    IC50df[i,"IC50"] = 3.2
  }
  else if (i==180) {
    IC50df[i,"IC50"] = 3.1
  }
  else if (i==185) {
    IC50df[i,"IC50"] = 3.2
  }
  else if (i==188) {
    IC50df[i,"IC50"] = 3.4
  }
  else if (i==192) {
    IC50df[i,"IC50"] = 3.7
  }
  else if (i==193) {
    IC50df[i,"IC50"] = 2.7
    # plot(log10(x),y)
    # lines(log10(x),predict(m),lty=2,col="red",lwd=3)
  }
  else if (i==194) {
    IC50df[i,"IC50"] = 2.6
    # plot(log10(x),y)
    # lines(log10(x),predict(m),lty=2,col="red",lwd=3)
  }
  else if (i==198) {
    IC50df[i,"IC50"] = 3.8
  }
  else if (i==206) {
    IC50df[i,"IC50"] = 3.75
  }
  else {
    m=nls(y~1/(1+(log10(x)/a)^b),start=c(a=2.5,b=2),control=list(maxiter=500))
    # plot(log10(x),y)
    # lines(log10(x),predict(m),lty=2,col="red",lwd=3)
    IC50 = coef(m)[1]
    IC50df[i,"IC50"] = IC50
  }
  
  
}

sIC50 = IC50df[order(IC50df$IC50),]

spIC50 = data.frame()

for (i in 1:nrow(sIC50)) {
  
  spIC50 = rbind(
    spIC50,
    data.frame(genotype=sIC50[i,1],n=i,noriginal=sIC50[i,'noriginal'],dose=doses,val=as.numeric(sIC50[i,13:22]))
  )
  
}

ggplot(spIC50,aes(x=log10(dose),y=n))+theme_bw()+
  geom_raster(aes(fill=val))+
  xlab('Dose [nM]')+ylab('Resistance Variant')+
  ggtitle('Imatinib Dose Response\nfor BCR-ABL Variants')+
  scale_x_continuous(breaks=c(1,2,3,4),labels=parse(text=c('10^1','10^2','10^3','10^4')))+
  scale_y_continuous(breaks=c(0,max(spIC50$n)),labels=c('Low\nResistance\n','High\nResistance\n'))+
  scale_fill_continuous(name='Viability',breaks=c(0,1))+
  theme(
    plot.title=element_text(hjust=0.5,size=26,face='bold'),
    axis.title=element_text(size=26,face='bold'),
    axis.text=element_text(size=20,color='black'),
    axis.text.y=element_text(hjust=0.5,angle=90),
    legend.title=element_text(size=20),
    legend.text=element_text(size=15)
  )

ggsave("IC50HeatMap2.png")

#### Growth Heat Map

rm(list=ls())

growth = read.csv('GrowthCurveHeatMap.csv')
subgrowth = growth[,1:12]

# Sort by growth rate

days = c(0,3,5,7:13)

for (i in 1:nrow(subgrowth)) {
  
  cells = as.double(subgrowth[i,3:12])
  mdl = lm(log(cells[5:10]) ~ days[5:10])
  m = mdl$coefficients[2]
  b = mdl$coefficients[1]
  subgrowth[i,'r'] = mdl$coefficients[2]
  
  # plot(x=days,y=cells)
  # 
  # tval = seq(0,13,0.1)
  # count = exp(b)*exp(m*tval)
  # lines(tval,count)
  
}

subgrowth$rrank = rank(subgrowth$r,ties.method='first')
subgrowth$finalrank = rank(subgrowth$Day13,ties.method='first')

pgrowth = data.frame()

for (i in 1:nrow(subgrowth)) {
  
  pgrowth = rbind(
    pgrowth,
    data.frame(genotype=subgrowth[i,1],n=i,day=days,val=as.numeric(subgrowth[i,3:12]),rrank=subgrowth[i,'rrank'],finalrank = subgrowth[i,'finalrank'])
  )
  
}

ggplot(pgrowth,aes(x=as.factor(day),y=finalrank))+theme_bw()+
  geom_raster(aes(fill=log10(val)))+
  scale_fill_gradient(name='Cell\nCount',low='white',high='red',breaks=c(6,8),labels=parse(text=c('10^6','10^8')))+
  scale_x_discrete(breaks=c(0,13))+
  scale_y_continuous(breaks=c(1,max(pgrowth$finalrank)-1),
                   labels=c("Low\nGrowth Rate","High\nGrowth Rate"))+
  ggtitle('BaF3 Cell Counts\nfor BCR-ABL Variants')+
  xlab('Time [days]')+ylab('Resistance Variant')+
  theme(
    plot.title=element_text(hjust=0.5,size=26,face='bold'),
    axis.title=element_text(size=26,face='bold'),
    axis.text=element_text(size=20,color='black'),
    axis.text.y=element_text(angle=90,hjust=0.5),
    legend.title=element_text(size=20),
    legend.text=element_text(size=15),
    legend.title.align=0.5
  )

ggsave("GrowthHeatMap2.png")
