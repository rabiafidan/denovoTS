#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

a=read.delim("20130606_g1k_3202_samples_ped_population.txt",sep=" ", header=T)
a2=data.frame(apply(a, 2,as.character),stringsAsFactors=FALSE)
a=a2
rm(a2)

children=a[ a$MotherID!="0" & a$FatherID !=0,  ]

cnames=c() #children names
mnames=c() #mother names
fnames=c() #father names
weird=c()
pops=c() #population
sup_pop=c() #super population

for (i in 1:dim(children)[1]){
  c=children$SampleID[i]
  p=children$Population[i]
  sp=children$Superpopulation[i]
  m=a[a$SampleID ==children[i,"MotherID"],"SampleID"]
  f=a[a$SampleID == children[i,"FatherID"],"SampleID"]
  if(length(m)==1 & length(f)==1){
    mnames=c(mnames,m)
    fnames=c(fnames,f)
    cnames=c(cnames,children$SampleID[i])
    pops=c(pops,p)
    sup_pop=c(sup_pop,sp)
  } else{weird=c(weird,children$SampleID[i]) }
}

write.table(data.frame("a"=cnames),file="children.txt",col.names = F,row.names = F,quote = F)
write.table(data.frame("a"=fnames),file="father.txt",col.names = F,row.names = F,quote = F)
write.table(data.frame("a"=mnames),file="mother.txt",col.names = F,row.names = F,quote = F)
write.table(data.frame("a"=pops),file="pops.txt",col.names = F,row.names = F,quote = F)
write.table(data.frame("a"=sup_pop),file="sup_pops.txt",col.names = F,row.names = F,quote = F)
write.table(data.frame("a"=unique(c(cnames,fnames,mnames)),file="inds.txt",col.names = F,row.names = F,quote = F)
