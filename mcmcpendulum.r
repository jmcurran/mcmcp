alleles<-list(L1=c(17,20,22,23),L2=c(16,17,18),L3=c(9,11,12,13),L4=c(13,14,15), L5=c(6,7,8,9.3))
areas<-list(L1=c(2863, 1630, 1024, 1005),L2=c(1910,953,398),L3=c(1633,1906,475,874),L4=c(1753,2212,857),L5=c(1514,651,1410,512))


G<-matrix(c(
1,1,2,3,
1,2,1,3,
1,2,2,3,
1,2,3,3,
1,3,1,2,
1,3,2,2,
1,3,2,3,
2,2,1,3,
2,3,1,1,
2,3,1,2,
2,3,1,3,
3,3,1,2), byrow=T, nrow=12)

for(i in 1:12)
	G[i,]<-G[i,c(3,4,1,2)]


Total<-lapply(areas,sum)

rg<-function(loc)
{
	if(loc==1 | loc==3 | loc==5){
		i1<-floor(3*runif(1))+1
		i2<-floor((3-i1)*runif(1))+i1+1
		g1<-c(i1,i2)
		g1<-c(g1,(1:4)[-g1])
		return(g1)
	}else{              
		i1<-floor(12*runif(1))+1
		g2<-G[i1,]
		return(g2)
	}
}

g0<-list(g1=rg(1),g2=rg(2),g3=rg(3),g4=rg(4),g5=rg(5))
mx0<-runif(1,0,0.5)
l0<-dchisq(loglik(1,mx0,g0$g1)+loglik(2,mx0,g0$g2)+loglik(3,mx0,g0$g3)+loglik(4,mx0,g0$g4)+loglik(5,mx0,g0$g5),19,log=T)

mx.sample<-NULL
g.sample<-NULL
naccept<-0
N<-0

while(naccept<=11000)
{
	g1<-g0
	l<-1+floor(5*runif(1))
	if(l==1)
		g1$g1=rg(1)
	else if(l==2)
		g1$g2=rg(2)
	else if(l==3)
		g1$g3=rg(3)
	else if(l==4)
		g1$g4=rg(4)
	else
		g1$g5=rg(5)
   
	mx1<-runif(1,0,0.5)
   	l1<-dchisq(loglik(1,mx1,g1$g1)+loglik(2,mx1,g1$g2)+loglik(3,mx1,g1$g3)+loglik(4,mx1,g1$g4)+loglik(5,mx1,g1$g5),19,log=T)

   	p<-min(0,l1-l0)
   	u<-log(runif(1))
    
   	if(u<p){
		g0<-g1
       	mx0<-mx1
       	l0<-l1
      	
      if(naccept>1000)
      {
           mx.sample<-c(mx.sample,mx1)
           g.sample<-rbind(g.sample,c(alleles$L1[g1$g1[1:2]],alleles$L2[g1$g2[1:2]]
												,alleles$L3[g1$g3[1:2]],alleles$L4[g1$g4[1:2]],alleles$L5[g1$g5[1:2]]))
      }
       naccept<-naccept+1
   }
   N<-N+1
}

naccept/N
