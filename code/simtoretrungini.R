# Simulate a field with dwindling SD, does our model fit it okay?






rm(list=ls())

root="C:/Users/royburst/Documents/tmb_inla_sim"
setwd(root)
source('code/utils.R')

library(seegMBG)
library(boot)
library(raster)
library(ggplot2)

###########
# data simulation
set.seed(123445)
simobj= mortsim(nu         = 2  ,       
                betas      = c(-5,0)      ,  
                scale      = .2            ,  
                Sigma2     = c(3,1,.5,.25),  
                rho        = 0.9          ,  
                l          = 51           ,  
                n_clusters = 75          ,  
                n_periods  = 4            ,  
                mean.exposure.months = 1000,  
                extent = c(0,1,0,1)       ,  
                ncovariates = 1           ,  
                seed   = NULL             ,
                returnall=TRUE            )
                
colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))
plot(simobj$r.true.mr,zlim=c(0,0.65),col=colz,nc=4,legend=F,xaxt='n',yaxt='n')         
plot(simobj$r.true.mr[[1]],zlim=c(0,0.65),col=colz,nc=4,legend=T,xaxt='n',yaxt='n')         

par(mfrow=c(1,4))
for(i in 1:4)
hist(as.vector(simobj$r.true.mr[[i]]),xlim=c(0,1),main="",xlab="",breaks=50,freq=F)

# plot simulated data
simobj$d$observed_mr=as.numeric(simobj$d$mr)
ggplot(simobj$d,aes(x=x,y=y))+ #,colour=factor(period)))+ # scale_colour_manual(values=rev(c('#f6ea8c', '#f26d5b', '#c03546', '#492540')))+
  geom_jitter(aes(size=observed_mr),shape=21,stroke=1.3,alpha=.5)+ #,width=.01,height=.01)+
  theme_bw()+ 
  scale_size(range = c(0, 15)) +
  facet_wrap(~period,nc=4)+ylab('')+xlab('')+
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())

# fit the model in INLA
mod=fitinla(simobj=simobj)
names(mod$mean_ras)<-paste("Period",1:4)
colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))
plot(mod$mean_ras,zlim=c(0,0.65),col=colz,nc=4,legend=F,xaxt='n',yaxt='n')         
plot(mod$mean_ras[[1]],zlim=c(0,0.6),col=colz,nc=4,legend=T,xaxt='n',yaxt='n')         


## What inequality underlie these true surfaces?

set.seed(12345)
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

for(i in 1:4){
    x=as.vector(simobj$r.true.mr[[i]])
    #x=plogis(x)
    
    res$gini[res$t==i]           =  
      IID(x,alpha=1,beta=1)
    res$madmed[res$t==i]         =     
      madmed(x,weights=rep(1,length(x)))
    res$range1090[res$t==i]       = 
      unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
    res$rangeratio1090[res$t==i]  =  
      unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
    res$variance[res$t==i]       =      
      var(x)
    res$MAD[res$t==i]            =  
      mad(x)
    res$cv[res$t==i]            =  
      mean(x)/max(x)
    
  
}
restrue=res
t(round(res,5))


## How do the inequalities in the data compare to the true surface?
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

for(i in 1:4){
  x=as.vector(simobj$d$observed_mr[simobj$d$period==i])
  #x=plogis(x)
  
  res$gini[res$t==i]           =  
    IID(x,alpha=1,beta=1)
  res$madmed[res$t==i]         =     
    madmed(x,weights=rep(1,length(x)))
  res$range1090[res$t==i]       = 
    unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
  res$rangeratio1090[res$t==i]  =  
    unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
  res$variance[res$t==i]       =      
    var(x)
  res$MAD[res$t==i]            =  
    mad(x)
  res$cv[res$t==i]            =  
    mean(x)/max(x)
  
  
}
resdata=res
t(round(res,5))


## How do the inequalities in the data compare to the estimated surface?
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

for(i in 1:4){
  x=as.vector(mod$mean_ras[[i]])
  #x=plogis(x)
  
  res$gini[res$t==i]           =  
    IID(x,alpha=1,beta=1)
  res$madmed[res$t==i]         =     
    madmed(x,weights=rep(1,length(x)))
  res$range1090[res$t==i]       = 
    unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
  res$rangeratio1090[res$t==i]  =  
    unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
  res$variance[res$t==i]       =      
    var(x)
  res$MAD[res$t==i]            =  
    mad(x)
  res$cv[res$t==i]            =  
    mean(x)/max(x)
  
  
}
resest=res
t(round(res,5))

resest$source   = 'Estimate'
resdata$source  = 'Observed Data'
restrue$source  = 'Truth'

res=rbind(resest,resdata,restrue)

ggplot(res,aes(x=t,y=gini,colour=source))+
  geom_line(size=2)+theme_bw()

require(gridExtra)
grid.arrange(
ggplot(res,aes(x=t,y=madmed,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ggplot(res,aes(x=t,y=range1090,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ggplot(res,aes(x=t,y=rangeratio1090,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ggplot(res,aes(x=t,y=variance,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ggplot(res,aes(x=t,y=MAD,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ggplot(res,aes(x=t,y=cv,colour=source))+
  geom_line(size=2)+theme_bw() + theme(legend.position="none"),
ncol=3
)



######################################
### TIME WORKS, HOW ABOUT SPACE with different SDs?
## Assume instead of 4 times we have 1 area with each of these, does our model fit okay?
xmin=0;xmax=1;ymin=0;ymax=1

rr=simobj$r.true.mr 
r1=raster(extent(xmin,xmax/2,ymax/2,ymax),nrows=51,ncols=51)
r2=raster(extent(xmax/2,xmax,ymax/2,ymax),nrows=51,ncols=51)
r3=raster(extent(xmin,xmax/2,ymin,ymax/2),nrows=51,ncols=51)
r4=raster(extent(xmax/2,xmax,ymin,ymax/2),nrows=51,ncols=51)
r1[is.na(r1)]=1
r2[is.na(r2)]=1
r3[is.na(r3)]=1
r4[is.na(r4)]=1
r1=insertRaster(r1,cbind(as.vector(rr[[1]])))
r2=insertRaster(r2,cbind(as.vector(rr[[2]])))
r3=insertRaster(r3,cbind(as.vector(rr[[3]])))
r4=insertRaster(r4,cbind(as.vector(rr[[4]])))

newr=merge(r1,r2,r3,r4)



# add country IDS
cntrysOrig = raster(matrix(c(1,3,2,4),nrow=sqrt(4),ncol=sqrt(4)))
cntrys     = resample(cntrysOrig,newr,method='ngb')  
cntryPoly  = rasterToPolygons(cntrysOrig)

#plot
plot(newr,zlim=c(0,0.65),col=colz,legend=T,xaxt='n',yaxt='n')         
lines(cntryPoly)
text(coordinates(cntryPoly), label=cntryPoly@data[,1])

# sample 300 clusters
set.seed(12344)
d=data.table(x=runif(300,0,1),y=runif(300,0,1))

d$p=raster::extract(newr,cbind(d$x,d$y)) 
d$exposures=round(abs(rnorm(n=300,mean=1000,sd=1000/5)))
d$deaths <- rbinom(300,size=d$exposures, prob=d$p)
d$mr = d$death/d$exposure
d$observed_mr=d$mr

ggplot(d,aes(x=x,y=y))+ 
  geom_jitter(aes(size=observed_mr),shape=21,stroke=1.3,alpha=.5)+ #,width=.01,height=.01)+
  theme_bw()+ 
  scale_size(range = c(0, 20)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_hline(yintercept=.5)+  geom_vline(xintercept=.5)


# fit inla model.. 
fit=fitinla2()

plot(fit,zlim=c(0,.85),col=colz,legend=T,xaxt='n',yaxt='n')         
lines(cntryPoly)
text(coordinates(cntryPoly), label=cntryPoly@data[,1])


# res stack by country
fits = stack(insertRaster(rr,cbind(as.vector(crop(fit,extent(xmin,xmax/2,ymax/2,ymax))))),
             insertRaster(rr,cbind(as.vector(crop(fit,extent(xmax/2,xmax,ymax/2,ymax))))),
             insertRaster(rr,cbind(as.vector(crop(fit,extent(xmin,xmax/2,ymin,ymax/2))))),
             insertRaster(rr,cbind(as.vector(crop(fit,extent(xmax/2,xmax,ymin,ymax/2))))))

names(fits)<-names(rr)<-paste('Country',1:4)


##################################
# see if we get good ginis
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

for(i in 1:4){
  x=as.vector(rr[[i]])
  #x=plogis(x)
  
  res$gini[res$t==i]           =  
    IID(x,alpha=1,beta=1)
  res$madmed[res$t==i]         =     
    madmed(x,weights=rep(1,length(x)))
  res$range1090[res$t==i]       = 
    unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
  res$rangeratio1090[res$t==i]  =  
    unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
  res$variance[res$t==i]       =      
    var(x)
  res$MAD[res$t==i]            =  
    mad(x)
  res$cv[res$t==i]            =  
    mean(x)/max(x)
  
  
}
restrue=res
t(round(res,5))


## How do the inequalities in the data compare to the true surface?
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

d$country=raster::extract(cntrys,cbind(d$x,d$y))
for(i in 1:4){
  x=as.vector(d$observed_mr[d$country==i])
  #x=plogis(x)
  
  res$gini[res$t==i]           =  
    IID(x,alpha=1,beta=1)
  res$madmed[res$t==i]         =     
    madmed(x,weights=rep(1,length(x)))
  res$range1090[res$t==i]       = 
    unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
  res$rangeratio1090[res$t==i]  =  
    unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
  res$variance[res$t==i]       =      
    var(x)
  res$MAD[res$t==i]            =  
    mad(x)
  res$cv[res$t==i]            =  
    mean(x)/max(x)
  
  
}
resdata=res
t(round(res,5))


## How do the inequalities in the data compare to the estimated surface?
res=data.table(t=1:4,
               gini=NA,madmed=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA,cv=NA)

for(i in 1:4){
  x=as.vector(fits[[i]])
  #x=plogis(x)
  
  res$gini[res$t==i]           =  
    IID(x,alpha=1,beta=1)
  res$madmed[res$t==i]         =     
    madmed(x,weights=rep(1,length(x)))
  res$range1090[res$t==i]       = 
    unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
  res$rangeratio1090[res$t==i]  =  
    unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
  res$variance[res$t==i]       =      
    var(x)
  res$MAD[res$t==i]            =  
    mad(x)
  res$cv[res$t==i]            =  
    mean(x)/max(x)
  
  
}
resest=res
t(round(res,5))

resest$source   = 'Estimate'
resdata$source  = 'Observed Data'
restrue$source  = 'Truth'

res=rbind(resest,resdata,restrue)
res$c=res$t

ggplot(res,aes(x=c,y=gini,colour=source))+
  geom_point(size=40,alpha=.7,shape=95)+theme_bw()


res2=subset(res,source!='Observed Data')
grid.arrange(
ggplot(res2,aes(x=c,y=madmed,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ggplot(res2,aes(x=c,y=range1090,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ggplot(res2,aes(x=c,y=rangeratio1090,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ggplot(res2,aes(x=c,y=variance,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ggplot(res2,aes(x=c,y=MAD,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ggplot(res2,aes(x=c,y=cv,colour=source))+
  geom_point(size=15,alpha=.7,shape=95)+theme_bw() + theme(legend.position="none"),
ncol=3)

## CHECK OUT CROSSSECTIONS
library(scales)

par(mfrow=c(2,1))
# look at fun cross-sections
plot(x=1,y=1,xlim=c(1,102),ylim=c(0,1),xlab='y',ylab='5q0',main='truth')
for(i in 1:102){
 # plot((as.matrix(fit)[,i]))
  lines((as.matrix(newr))[,i],col=alpha("black",(30+i)/150),lwd=1.5)
  
}
# look at fun cross-sections
plot(x=1,y=1,xlim=c(1,102),ylim=c(0,1),xlab='y',ylab='5q0',main='fit')
for(i in 1:102){
  # plot((as.matrix(fit)[,i]))
  lines((as.matrix(fit))[,i],col=alpha("black",(30+i)/150),lwd=1.5)
}

  par(mfrow=c(2,1))
  # look at fun cross-sections
  plot(x=1,y=1,xlim=c(1,102),ylim=c(0,1),xlab='x',ylab='5q0',main='truth')
  for(i in 1:102){
    # plot((as.matrix(fit)[,i]))
    lines((as.matrix(newr))[i,],col=alpha("black",(30+i)/150),lwd=1.5)
    
  }
  
  # look at fun cross-sections
  plot(x=1,y=1,xlim=c(1,102),ylim=c(0,1),xlab='x',ylab='5q0',main='fit')
  for(i in 1:102){
    # plot((as.matrix(fit)[,i]))
    lines((as.matrix(fit))[i,],col=alpha("black",(30+i)/150),lwd=1.5)
    
  }

  
  
  
  
  
  
####################################################  

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  

