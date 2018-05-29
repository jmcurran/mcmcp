#TWO_CASE1
x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1/1.4
x<-c(x1,x2)
y<-c(y1,y2)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[400]),c(0,y,0),col="light blue")
axis(1,at=c(0,6),labels=c("1","2"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(2.5,0.3,expression(frac(P[1],P[2])>=ClearMajor2Allele),adj=0)
text(2.5,0.25,"Major not known",adj=0)
text(2.5,0.22,"Minor not known",adj=0)

#THREE_CASE1
x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1/2.4
x3<-x2+6
x<-c(x1,x2,x3)
y<-c(y1,y2,y2*0.95)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[600]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12),labels=c("1","2","3"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(6,0.3,expression(frac(P[1],P[2]+P[3])>=ClearMajor3Allele),adj=0)
text(6,0.25,expression(frac(P[3],P[2])<=PrefAmp),adj=0)
text(6,0.22,"Major not known",adj=0)
text(6,0.20,"Minor not known",adj=0)

#THREE_CASE2
x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1/(2.4*8/3)
x3<-x2+6
x<-c(x1,x2,x3)
y<-c(y1,y2,y2*5/3)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[600]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12),labels=c("1","2","3"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(5,0.3,expression(frac(P[1],P[2]+P[3])>=ClearMajor3Allele),adj=0)
text(5,0.25,expression(frac(P[3],P[2])>PrefAmp),adj=0)
text(5,0.22,"Major not known",adj=0)
text(5,0.20,"Minor not known",adj=0)

#THREE_CASE3
x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1/2.4
x3<-x2+6
x<-c(x1,x2,x3)
y<-c(y1,y2,y2*2)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[600]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12),labels=c("1","2","3"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(5,0.3,expression(frac(P[1],P[2]+P[3])<ClearMajor3Allele),adj=0)
text(5,0.25,expression(frac(P[3],P[2])>PrefAmp),adj=0)
text(5,0.22,"Major not known",adj=0)
text(5,0.20,"Minor not known",adj=0)

#THREE_CASE4
x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1*0.6
x3<-x2+6
x<-c(x1,x2,x3)
y<-c(y1,y2,y2*1.1)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[600]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12),labels=c("1","2","3"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(5,0.4,expression(frac(P[1],P[2]+P[3])<ClearMajor3Allele),adj=0)
text(5,0.35,expression(paste("or"~~frac(P[3],P[2])>PrefAmp)),adj=0)
text(5,0.32,"or major known",adj=0)
text(5,0.30,"or minor known",adj=0)


#FOUR_CASE1

x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1
y3<-y1*.5
y4<-y1*.5
x3<-x2+6
x4<-x3+6
x<-c(x1,x2,x3,x4)
y<-c(y1,y2,y3,y4)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[800]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12,18),labels=c("1","2","3","4"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(18,0.1,expression(P[4]))

text(9,0.32,expression(frac(P[1]+P[2],P[3]+P[4])>=MajMinDiff),adj=0)
text(9,0.27,expression(frac(P[2],P[1])>PrefAmp),adj=0)
text(9,0.24,"Major not known",adj=0)
text(9,0.22,"Minor not known",adj=0)

#FOUR_CASE2

x1<-seq(-3,3,length=200)
y1<-dnorm(x1)
x2<-x1+6
y2<-y1
y1<-0.5*y1
y3<-y1*.8
y4<-y2*.8
x3<-x2+6
x4<-x3+6
x<-c(x1,x2,x3,x4)
y<-c(y1,y2,y3,y4)
par(xaxs="i",yaxs="i")
plot(x,y,type="l",axes=F,xlab="",ylab="",ylim=c(0,0.42))
polygon(c(x[1],x,x[800]),c(0,y,0),col="light blue")
axis(1,at=c(0,6,12,18),labels=c("1","2","3","4"))
abline(h=0)

text(0,0.1,expression(P[1]))
text(6,0.1,expression(P[2]))
text(12,0.1,expression(P[3]))
text(18,0.1,expression(P[4]))

text(9,0.32,expression(frac(P[1]+P[2],P[3]+P[4])<MajMinDiff),adj=0)
text(9,0.27,expression(paste("or",~~frac(P[2],P[1])<=PrefAmp)),adj=0)
text(9,0.24,"or major not known",adj=0)
text(9,0.22,"or minor not known",adj=0)
