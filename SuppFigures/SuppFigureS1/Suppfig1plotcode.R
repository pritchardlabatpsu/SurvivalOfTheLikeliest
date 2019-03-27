setwd("~/MutationFreq_ABL_IC50/Figures/Figure 2")
#install.packages("GGally")
library("ggplot2")
library("GGally")
source("ggcorplot.R")
icdat=read.csv("Supp_fig1_data.csv",header=T,stringsAsFactors = F)
ggpairs(icdat[,2:6],lower=list(continuous="smooth"),
        diag=list(continuous="bar" ), 
        upper=list(), axisLabels='show')
