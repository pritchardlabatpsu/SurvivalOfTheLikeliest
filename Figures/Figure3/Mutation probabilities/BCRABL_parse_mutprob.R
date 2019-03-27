rm(list=ls())

rm(list=ls())

setwd("")

##########################Parse BCR ABL seq to data frame################################################################
#Do for any patient
ba=readLines('BCRABLseqcodon.txt')
ba=ba[!ba==""]#Removes lines of whitespace
bacod=ba[grep("^[ATCG]{3}",ba)]#separates DNA
batran=ba[-grep("^[ATCG]{3}",ba)]#separates AA 3 letter code
#Coding seq to 1 string
coding=paste(bacod,collapse='')#makes 1 string
coding=gsub(" ","",coding,fixed=T)
#Translation to 1 string
tran=paste(batran,collapse='')
tran=gsub(" ","",tran,fixed=T)#removes whitespace
codframe=data.frame(unlist(lapply(1:nchar(coding),function(i){substring(coding,i,i)})))
x=nrow(codframe)/3
codframe[,2]=rep(c(1:x),each=3)
tranframe=data.frame(unlist(lapply(1:x,function(i){rep(substring(tran,((i*3)-2),(i*3)),times=3)})))
codframe[,3]=tranframe
colnames(codframe)[1:3]=c("NT","Codon","AA")
###########################################################################################################################
####codframe is a datframe that contains the gene BCR-ABL##################################################################
###########################################################################################################################
#read in data
mutprop=read.csv("Mutsubmat_broadExac.csv",stringsAsFactors=F, header=T)
codtable=read.csv('codtable.csv',stringsAsFactors=F,header=T)
##############################################################################################################################
seqcapture=codframe[3271:4470,]
seqcapture[1:1200,"codpos"]=rep(c(1,2,3),times=400)
seqcapture$NT=as.character(seqcapture$NT)
seqcapture[1:nrow(seqcapture),"ablcod"]=seqcapture$Codon-901
#make a separate row for each substitution
seqmutframe=seqcapture[rep(seq_len(nrow(seqcapture)), each=3),]
seqmutframe$NT=as.character(seqmutframe$NT)
As=nrow(seqcapture[seqcapture$NT=="A",])
Cs=nrow(seqcapture[seqcapture$NT=="C",])
Gs=nrow(seqcapture[seqcapture$NT=="G",])
Ts=nrow(seqcapture[seqcapture$NT=="T",])

seqmutframe[1:nrow(seqmutframe),"probtot"]=0
for (ii in 1:nrow(seqcapture)){
  if(seqcapture$NT[ii]=="A"){seqmutframe[((3*ii)-2):((3*ii)),"probtot"]=mutprop[mutprop$ref==seqcapture$NT[ii],"prop"]*(1/As)
  seqmutframe[((3*ii)-2):((3*ii)),"var"]=mutprop[mutprop$ref==seqcapture$NT[ii],"var"]}
  if(seqcapture$NT[ii]=="C"){seqmutframe[((3*ii)-2):((3*ii)),"probtot"]=mutprop[mutprop$ref==seqcapture$NT[ii],"prop"]*(1/Cs)
  seqmutframe[((3*ii)-2):((3*ii)),"var"]=mutprop[mutprop$ref==seqcapture$NT[ii],"var"]}
  if(seqcapture$NT[ii]=="T"){seqmutframe[((3*ii)-2):((3*ii)),"probtot"]=mutprop[mutprop$ref==seqcapture$NT[ii],"prop"]*(1/Ts)
  seqmutframe[((3*ii)-2):((3*ii)),"var"]=mutprop[mutprop$ref==seqcapture$NT[ii],"var"]}
  if(seqcapture$NT[ii]=="G"){seqmutframe[((3*ii)-2):((3*ii)),"probtot"]=mutprop[mutprop$ref==seqcapture$NT[ii],"prop"]*(1/Gs)
  seqmutframe[((3*ii)-2):((3*ii)),"var"]=mutprop[mutprop$ref==seqcapture$NT[ii],"var"]} 
}
for(jj in 1:nrow(seqmutframe)){
  getcodonnum=seqmutframe[jj,'Codon']
  codWtNt=seqcapture[seqcapture$Codon==getcodonnum,"NT"]
  codposvar=as.list(c(seqmutframe[jj,'codpos'],seqmutframe[jj,'var']))
  if(codposvar[1]=="1"){codvarNt=c(codposvar[2],codWtNt[2],codWtNt[3])
  seqmutframe[jj,'varAA']=codtable[codtable$cod1==codvarNt[1]&codtable$cod2==codvarNt[2]&codtable$cod3==codvarNt[3],'AA']
  }
  if(codposvar[1]=="2"){codvarNt=c(codWtNt[1],codposvar[2],codWtNt[3])
  seqmutframe[jj,'varAA']=codtable[codtable$cod1==codvarNt[1]&codtable$cod2==codvarNt[2]&codtable$cod3==codvarNt[3],'AA']
  }
  if(codposvar[1]=="3"){codvarNt=c(codWtNt[1],codWtNt[2],codposvar[2])
  seqmutframe[jj,'varAA']=codtable[codtable$cod1==codvarNt[1]&codtable$cod2==codvarNt[2]&codtable$cod3==codvarNt[3],'AA']
  }
}


write.csv(seqmutframe,"MutProbABL1.csv")

