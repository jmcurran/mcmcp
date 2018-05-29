m<-0.25
q95<-0.4
qhat<-1
alpha.lower<-1
alpha.upper<-100
alpha.mid<-alpha.lower+0.5*(alpha.upper-alpha.lower)

beta.lower<-(alpha.lower*(1-m)+2*m-1)/m
beta.mid<-(alpha.mid*(1-m)+2*m-1)/m
beta.upper<-(alpha.upper*(1-m)+2*m-1)/m

q.lower<-qbeta(0.95,alpha.lower,beta.lower)
q.mid<-qbeta(0.95,alpha.mid,beta.mid)
q.upper<-qbeta(0.95,alpha.upper,beta.upper)

while(abs(q.lower-q.upper)>1e-4){
	if(q.mid<q95){
		alpha.upper<-alpha.mid
	}else{
		alpha.lower<-alpha.mid
	}

	alpha.mid<-alpha.lower+0.5*(alpha.upper-alpha.lower)
		
	beta.lower<-(alpha.lower*(1-m)+2*m-1)/m
	beta.mid<-(alpha.mid*(1-m)+2*m-1)/m
	beta.upper<-(alpha.upper*(1-m)+2*m-1)/m
	
	q.lower<-qbeta(0.95,alpha.lower,beta.lower)
	q.mid<-qbeta(0.95,alpha.mid,beta.mid)
	q.upper<-qbeta(0.95,alpha.upper,beta.upper)
}
	
