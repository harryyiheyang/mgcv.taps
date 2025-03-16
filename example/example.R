require(mgcv)
devtools::load_all()
#set.seed(100)
dat <- gamSim(1,n=1000,scale=1,dist = "binary")

linear_changepoint=function(x,para){
softthres=function(x,b){
  y=x-b
  y[y<0]=0
  return(y)
}
knot=para
p=length(knot)
G=matrix(0,length(x),p)
for(i in 1:p){
G[,i]=softthres(x,knot[i])
}
return(cbind(1,x,G))
}

b<-bam(y~s(x0,bs="Amatern",k=10,m=300,xt=list(getA=linear_changepoint,para=0.5))+s(x1,bs="Amatern",k=10,m=300,xt=list(getA=linear_changepoint,para=0.4))
     +s(x2,bs="Amatern",k=10,m=300,xt=list(getA=linear_changepoint,para=c(0.24,0.5,0.7)))+s(x3,bs="Amatern",m=300,k=10),data=dat,method="REML",family="binomial",select=F)
fitb=summary(b)
plot(b,pages=1)
summary(b)
mgcv_smooth_wald(b)

b1<-bam(y~s(x0,bs="cosine")+s(x1,bs="cosine")
     +s(x2,bs="cosine")+s(x3,bs="cosine"),data=dat,method="REML",family="binomial",select=F)
plot(b1,pages=1)
summary(b1)
mgcv_smooth_wald(b1)
