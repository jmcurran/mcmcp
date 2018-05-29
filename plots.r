results.10K<-read.csv("E:/curran/Code/MCMCP/results10K.csv",sep=" ",head=F)
results.5K<-read.csv("E:/curran/Code/MCMCP//results5K.csv",sep=" ",head=F)
#loci<-c("D8S1179","D21S11","D7S820","CSF1PO","D3S1358","TH01","D13S317","D16S539","D2S1338","D19S433","vWA","TPOX","D18S51","AMELO","D5S818","FGA")
loci<-c("D3","vWA","FGA","D8","D21","D18","D5","D13","D7","D16","THO","TPOX","CSF")
n.loci<-length(loci)

profile<-read.csv("E:/curran/Code/MCMCP/wangTW.csv",head=F)
y.max<-max(profile[,2])

par(mfrow=c(4,4))
for(i in 1:n.loci){
  lwr<-(i-1)*4+1
  upr<-i*4
  peaks<-profile[lwr:upr,]
  heights<-peaks$V2[peaks$V2!=0]
  alleles<-peaks$V1[peaks$V2!=0]
  barplot(heights,names.arg=alleles,main=loci[i],ylim=c(0,1.1*y.max))
  box()
}



par(mfrow=c(4,4))
for(i in 1:n.loci){
  tbl.5k<-table(results.5K[,i])
  tbl.5k<-tbl.5k/sum(tbl.5k)
  tbl.10k<-table(results.10K[,i])
  tbl.10k<-tbl.10k/sum(tbl.10k)
  barplot(rbind(tbl.5k,tbl.10k),beside=T,ylim=c(0,1),las=2,col=c("red","blue"),main=loci[i],ylab="Probability")
}


par(mfrow=c(4,4))
for(i in 1:n.loci){
  tbl.5k<-table(results.5K[,i+n.loci])
  tbl.5k<-tbl.5k/sum(tbl.5k)
  tbl.10k<-table(results.10K[,i+n.loci])
  tbl.10k<-tbl.10k/sum(tbl.10k)
  barplot(rbind(tbl.5k,tbl.10k),beside=T,ylim=c(0,1),las=2,col=c("red","blue")
          main=loci[i],ylab="Probability")
}

i<-8
tbl.5k<-table(results.5K[,i+n.loci])
  tbl.5k<-tbl.5k/sum(tbl.5k)*100
  tbl.10k<-table(results.10K[,i+n.loci])
  tbl.10k<-tbl.10k/sum(tbl.10k)*100
  barplot(rbind(tbl.5k,tbl.10k),beside=T,ylim=c(0,100),las=2,col=c("red","blue"),main=loci[i])
box()
legend(1,95,fill=c("red","blue"),legend=c("50 million","100 million"))


