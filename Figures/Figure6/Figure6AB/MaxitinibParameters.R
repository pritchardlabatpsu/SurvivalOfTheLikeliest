
rm(list=ls())

setwd("C:/Users/Scott/Box Sync/MutBias/Figures/Figure5/Figure5DE/")

imat.df = read.csv("Combined_data_frame_IC_Mutprob_abundance.csv")

mut.prob = imat.df$Mutation.Probability/sum(imat.df$Mutation.Probability)
imat.var = imat.df$Compound
sort.mut.prob = mut.prob[order(substr(imat.var,2,5))]
sort.imat.var = imat.var[order(substr(imat.var,2,5))]
names(sort.mut.prob) = as.character(sort.imat.var)


alpha.df = read.csv("ImatinibAlphaGeneratorResults011619.csv",header=F)
alphas = alpha.df[,4:ncol(alpha.df)]
mean.alphas = colMeans(alphas)

alpha.sen = mean.alphas[1]
alpha.res = mean.alphas[2:length(mean.alphas)]
names(alpha.res) = as.character(sort.imat.var)

res.df = data.frame(variant=sort.imat.var,prob=sort.mut.prob,alpha=alpha.res)
sort.res.df = res.df[order(-res.df$prob),]

out.df = rbind(
  data.frame(variant="WT",prob=NaN,alpha=alpha.sen),
  sort.res.df
)


# write.csv(out.df,"MaxitinibParameters012619.csv")
