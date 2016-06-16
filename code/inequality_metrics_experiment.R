 # Test Inequality Metrics on Fake data



rm(list=ls())

root="C:/Users/royburst/Documents/tmb_inla_sim"
setwd(root)
source('code/utils.R')

library(seegMBG)
library(boot)
library(raster)


l=51 # matrix length
countries=16 # must be a square




##################################
## SIMULATE SOME DATA

#  estimated u5m probability surface (logit space)
set.seed(12344)
r = makeRandomSurface(sd=2,l=l,scale=.2,offset=-5)
r = calc(r, fun = inv.logit)
plot(r)

# TODO need more draws


#  population raster
set.seed(12344)
p = makeRandomSurface(sd=1,l=l,scale=.05)
p = calc(p, fun = exp)*10
p = insertRaster(p,data.frame(rpois(l**2,(as.vector(p)))))
names(p)='population'


#  country ID raster surface
cntrysOrig = raster(matrix(seq(1:countries),nrow=sqrt(countries),ncol=sqrt(countries)))
cntrys     = resample(cntrysOrig,r,method='ngb')  
cntryPoly  = rasterToPolygons(cntrysOrig)

onecountry = resample(raster(matrix(1,nrow=1,ncol=1)),r,method='ngb') 


# population weighting
pop_wt = makePopulationWeightsRaster(cntrys,p)
pop_wtall = makePopulationWeightsRaster(onecountry,p)

# we must not count areas with no pop as areas with no mort, so make NA
#pop_wt[pop_wt==0]=NA
#pop_wtall[pop_wtall==0]=NA


# pop weights that just 1 or 0 (o for zero pop)
nopopdrop_wt = pop_wt #!=0


############################ 
# CALCULATE INEQUALITY MEASURES
res = data.frame(

  

  
  RGR_nonwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                     weights = NULL,
                     group   = as.vector(extract(cntrys,cellIdx(r))),
                     fun     = RGR,hi=.9,lo=.1),
  range_nonwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                            weights = NULL,
                            group   = as.vector(extract(cntrys,cellIdx(r))),
                            fun     = range,hi=.9,lo=.1),

  madmed_nonwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                             weights = NULL,
                             group   = as.vector(extract(cntrys,cellIdx(r))),
                             fun     = madmed),
  cv_nonwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                         weights = NULL,
                         group   = as.vector(extract(cntrys,cellIdx(r))),
                         fun     = coefvar),

  gini_nonwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                           weights = NULL,
                           group   = as.vector(extract(cntrys,cellIdx(r))),
                           fun     = gini),
  
  
  RGR_popwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                       weights = as.vector(extract(p,cellIdx(r))),
                       group   = as.vector(extract(cntrys,cellIdx(r))),
                       fun     = RGR,hi=.9,lo=.1),
  
 
  
  range_popwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                         weights = as.vector(extract(p,cellIdx(r))),
                         group   = as.vector(extract(cntrys,cellIdx(r))),
                            fun     = range,hi=.9,lo=.1),
  

  
  madmed_popwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                             weights = as.vector(extract(p,cellIdx(r))),
                             group   = as.vector(extract(cntrys,cellIdx(r))),
                             fun     = madmed),
  
  
  cv_popwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                         weights = as.vector(extract(pop_wt,cellIdx(r))),
                         group   = as.vector(extract(cntrys,cellIdx(r))),
                         fun     = coefvar),
  
  gini_popwt_sim = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                           weights = as.vector(extract(pop_wt,cellIdx(r))),
                           group   = as.vector(extract(cntrys,cellIdx(r))),
                           fun     = gini)
)


pdf('output/sandbox/inequality_simulation3.pdf',width=15,height=15)
  ############################ 
  # PLOT SIMULATED DATA
  colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))
  par(mfrow=c(2,2))
  plot(r,main='mortality',col=colz,legend=T,xaxt='n',yaxt='n') 
  lines(cntryPoly)
  text(coordinates(cntryPoly), label=cntryPoly@data[,1])
  plot(p,main='population',legend=T,xaxt='n',yaxt='n') 
  lines(cntryPoly)
  text(coordinates(cntryPoly), label=cntryPoly@data[,1])
  
  plot(r*pop_wt,main='pop-weighted mortality',col=colz,legend=T,xaxt='n',yaxt='n') 
  lines(cntryPoly)
  text(coordinates(cntryPoly), label=cntryPoly@data[,1])
  
  
  ############################
  # country histograms non-weighted
  mx=max(as.matrix(extract(r,cellIdx(r))))
  mn=min(as.matrix(extract(r,cellIdx(r))))
  
  par(mfrow=c(sqrt(countries),sqrt(countries)))
  for(i in unique(as.vector(extract(cntrys,cellIdx(r)))))
    hist(as.matrix(extract(r,cellIdx(r)))[as.vector(extract(cntrys,cellIdx(r)))==i],
         main=i,xlab="",xlim=c(mn,mx),breaks=round(l**2/(countries*30),0))

  title("Non Weighted MR by Country",outer=TRUE)
  
  
  ############################
  # country histograms weighted cells
  mx=max(as.matrix(extract(r,cellIdx(r)))*as.matrix(extract(pop_wt,cellIdx(r))))
  mn=min(as.matrix(extract(r,cellIdx(r)))*as.matrix(extract(pop_wt,cellIdx(r))))
  
  par(mfrow=c(sqrt(countries),sqrt(countries)))
  for(i in unique(as.vector(extract(cntrys,cellIdx(r)))))
    hist(as.matrix(extract(r,cellIdx(r)))*as.matrix(extract(pop_wt,cellIdx(r)))[as.vector(extract(cntrys,cellIdx(r)))==i],
         main=i,xlab="",xlim=c(mn,mx),breaks=round(l**2/(countries*30),0))
  title("Pop-weighted MR by Country",outer=TRUE)
  
  
  
  ############################
  # plot table of results
  library(gridExtra)
  library(grid)
  
  result = res[order(-res$gini_popwt_sim),]
  mytheme <- gridExtra::ttheme_default(
    core = list(fg_params=list(cex = 1.0)),
    colhead = list(fg_params=list(cex = 1.0)),
    rowhead = list(fg_params=list(cex = 1.0)))
  
  plot(tableGrob(round(result,2),rows=rownames(result)),theme=ttheme_minimal())
  
  
  ###############################
  #### Heatmap ranks of each metric
  require(lattice)
  
  resrank = apply(result,2,rank)
  plot(levelplot((resrank),xlab='Country',ylab='Inequality Metric',main='Inequality Rank'))


  
  
  
  ################# 
  ## Look at sensitivity of range and range ratio to different quantile cutoffs. 
  
  his = rev(75:95)/100
  los = 1-his
  
  
  sensres = data.table()
  for(i in 1:length(los)){
    hi=his[i]
    lo=los[i]
    ratio = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                            weights = as.vector(extract(nopopdrop_wt,cellIdx(r))),
                            group   = as.vector(extract(cntrys,cellIdx(r))),
                            fun     = RGR,hi=hi,lo=lo)
    absolute = condSim(vals    = as.matrix(extract(r,cellIdx(r))),
                              weights = as.vector(extract(nopopdrop_wt,cellIdx(r))),
                              group   = as.vector(extract(cntrys,cellIdx(r))),
                              fun     = range,hi=hi,lo=lo)

    tmp=data.table(country=row.names(ratio),hi=hi,lo=lo,ratio=ratio,absolute=absolute)
  
    sensres=rbind(sensres,tmp)
  
  }
  sensres$ratio.V1[sensres$ratio.V1==Inf]=NA
  sensres$country=as.numeric(sensres$country)
  
  # plot as two line graphs
  require(ggplot2)
  labels = subset(sensres,hi==.75)
  
  r=ggplot(sensres,aes(x=hi,y=ratio.V1,colour=factor(country)))+
    geom_line(size=2,aes(alpha=.5))+
    geom_text(data=labels,aes(label=country,x=.75,y=ratio.V1),colour='black',size=4)+
    theme_bw()+
    scale_x_continuous(breaks = c(.75,.8,.85,.9,.95),labels = c('25%-75%','20%-80%','15%-85%','10%-90%','5%-95%'))+
    xlab('Quantile Range')+
    scale_y_log10(breaks=c(2:10,15,seq(20,100,10)))+
    ylab('Ratio (hi/low) (relative measure)')+
    ggtitle('Threshold sensitivity in RATIO metric') +
    theme(legend.position="none")
  
  
  a=ggplot(sensres,aes(x=hi,y=absolute.V1,group=country,colour=factor(country)))+
    geom_line(size=2,aes(alpha=.5))+
    geom_text(data=labels,aes(label=country,x=.75,y=absolute.V1),colour='black',size=4)+
    theme_bw()+
    scale_x_continuous(breaks = c(.75,.8,.85,.9,.95),labels = c('25%-75%','20%-80%','15%-85%','10%-90%','5%-95%'))+
    xlab('Quantile Range')+
    ylab('Difference (hi-low) (absolute measure)') + theme(legend.position="none") +
    ggtitle('Threshold sensitivity in ABSOLUTE metric') 
    
  
  
  
  grid.arrange(r,a,ncol=2)
  
  


dev.off()












##########################
# plots for slides
require(weights)
colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))
colz2=(colorRampPalette(c('#004e66','#e1eef6','#fcbe32',  '#ff5f2e'))(255))
colz3=(colorRampPalette(c('#00dffc', '#008c9e', '#005f6b', '#343838'))(255))

plot(r,main='Simulated Mortality Surface',col=colz,legend=T,xaxt='n',yaxt='n') 
#plot(1-r,main='Simulated survival Surface',col=colz,legend=T,xaxt='n',yaxt='n') 


plot(p,main='Simulated Population Surface',col=colz2,legend=T,xaxt='n',yaxt='n') 

raster::plot(r*pop_wtall*l**2,main='Simulated Pop-Wt Mortality',col=colz,legend=T,xaxt='n',yaxt='n',colNA='gray') 

#plot(r*pop_wt,main='pop-weighted mortality',col=colz,legend=T,xaxt='n',yaxt='n') 
par(mfrow=c(2,2))
nonmissing= !is.na(as.vector(pop_wtall))
hist(as.vector(r), main='Mortality Pixels',breaks=100,xlab='Simulated MR')


hist(logit(as.vector(r)), main='Logit Mortality Pixels',breaks=100,xlab='Simulated MR')

wtd.hist(as.vector(r)[nonmissing], 
         weight=as.vector(pop_wtall*l**2)[nonmissing],
         main='Mortality Pixels (Weighted)',
         breaks=100,
         xlab='Simulated Pop-weighted MR')


wtd.hist(logit(as.vector(r))[nonmissing], 
         weight=as.vector(pop_wtall*l**2)[nonmissing],
         main='Logit Mortality Pixels (Weighted)',
         breaks=100,
         xlab='Simulated Pop-weighted MR')








# lorenz
require(lawstat)
require(ineq)


par(mfrow=c(1,2))
plot(Lc(as.vector(r)[nonmissing],
        n=as.vector(pop_wtall)[nonmissing]),
     main="Pop Weighted Lorenz Curve",
     ylab="Cumulative Proportion of Mortality",
     xlab="Cumulative Proportion of Population")
plot(Lc(as.vector(r)),
     main="Unweighted Lorenz Curve",
     ylab="Cumulative Proportion of Mortality",
     xlab="Cumulative Proportion of Population")

gini(as.vector(r)[nonmissing],
     weight=as.vector(pop_wtall)[nonmissing])


gini(as.vector(r),weights=rep(1,l**2))


# survival

par(mfrow=c(1,2))
plot(Lc(as.vector(1-r)[nonmissing],
        n=as.vector(pop_wtall)[nonmissing]),
     main="Pop Weighted Lorenz Curve",
     ylab="Cumulative Proportion of Survival Risk",
     xlab="Cumulative Proportion of Population")
plot(Lc(as.vector(1-r)),
     main="Unweighted Lorenz Curve",
     ylab="Cumulative Proportion of Survival Risk",
     xlab="Cumulative Proportion of Population")

gini(as.vector(1-r)[nonmissing],
     weights=as.vector(pop_wtall)[nonmissing])

gini(as.vector(1-r),weights=rep(1,l**2))





# Expected Survival time
x= as.vector(r)
m = -(log(1-x)/5)
S = (1/m)-((exp(-5*m))/m)

hist(S,breaks=100)


par(mfrow=c(1,2))
plot(Lc(S[nonmissing],
        n=as.vector(pop_wtall)[nonmissing]),
     main="Pop Weighted Lorenz Curve",
     ylab="Cumulative Proportion of Expected Survival Time",
     xlab="Cumulative Proportion of Population")
plot(Lc(S),
     main="Unweighted Lorenz Curve",
     ylab="Cumulative Proportion of Expected Survival Time",
     xlab="Cumulative Proportion of Population")

gini(S[nonmissing],
     weights=as.vector(pop_wtall)[nonmissing])

gini(S,weights=rep(1,l**2))







##############
# Do a sensitivity analysis on changing alpha and beta on IID and IMD Metrics

i=1
x=as.vector(r)
res=data.table(alpha=rep(seq(1,4,.3),11),beta=rep(seq(0,1,.1),each=11),IMD=NA,IID=NA)
pb <- txtProgressBar(min = 1, max = nrow(res), style = 3)

for(a in unique(res$alpha)){
  for(b in unique(res$beta)){
    setTxtProgressBar(pb, i); i=i+1
    res$IID[res$alpha==a&res$beta==b]=IID(x,alpha=a,beta=b)
    res$IMD[res$alpha==a&res$beta==b]=IMD(x,alpha=a,beta=b)
  }
}
close(pb)

# or use outer
rIID=rasterFromXYZ(cbind(res$beta,res$alpha,res$IID))
rIMD=rasterFromXYZ(cbind(res$beta,res$alpha,res$IMD))
# plot them
colz=rev(colorRampPalette(c('#f9320c','#f9c00c', '#00b9f1'))(255))
  
par(mfrow=c(1,2))
plot(rIID,col=colz,main="IID",ylab=expression(alpha),xlab=expression(beta))
plot(rIMD,col=colz,main="IMD",ylab=expression(alpha),xlab=expression(beta),xlim=c(-.05,1.05))






###############
## plot r with changing variance and range
#  estimated u5m probability surface (logit space)
set.seed(12344)
r = 
r = calc(r, fun = inv.logit)
plot(r)

ranges=c(0.01,0.05,0.1,0.2,0.5)
vars  =c(0.25,0.50,1.0,2.0,3.0)
colz=rev(colorRampPalette(c('#383A3F','#E71D36', '#2EC4B6', '#EFFFE9'))(255))
par(mfrow=c(5,5),mar=c(0,0,0,0))

for(i in 1:5){
  for(j in 1:5){
    set.seed(12344)
    plot(calc(makeRandomSurface(sd=vars[i],l=l,scale=ranges[j],offset=-5),fun=inv.logit),
         col=colz,zlim=c(0,0.6),
         legend=F,
         xaxt='n',
         yaxt='n') 
}}
# also, do histograms


dev.off()

plot(calc(makeRandomSurface(sd=vars[i],l=l,scale=ranges[j],offset=-5),fun=inv.logit),
     col=colz,zlim=c(0,0.6),
     legend=T,
     xaxt='n',
     yaxt='n')


par(mfrow=c(5,5),mar=c(0,0,0,0))

for(i in 1:5){
  for(j in 1:5){
    set.seed(12344)
    hist(as.vector(calc(makeRandomSurface(sd=vars[i],l=l,scale=ranges[j],offset=-5),fun=inv.logit)),
         main="",
         xlim=c(0,.6),
         breaks=25,     xaxt='n',
         yaxt='n')
  }}





# see how IMD and IID change as SD and range change
ranges=c(0.01,0.05,0.1,0.2,0.5)
vars  =c(0.25,0.50,1.0,2.0,3.0)

pb <- txtProgressBar(min = 1, max = 25*25, style = 3)
z=0
ress=data.table()
for(i in 1:5){
  for(j in 1:5){
    
    
    set.seed(12344)
    
    res=data.table(alpha=rep(seq(1,3,.5),5),
                   beta=rep(seq(0,1,.25),each=5),
                   IMD=NA,IID=NA,range=ranges[j],sd=vars[i])
    
    x=as.vector(calc(makeRandomSurface(sd=vars[i],
                                       l=l,
                                       scale=ranges[j],
                                       offset=-5),fun=inv.logit))
    for(a in unique(res$alpha)){
      for(b in unique(res$beta)){
        setTxtProgressBar(pb, z); z=z+1
        res$IID[res$alpha==a&res$beta==b]=IID(x,alpha=a,beta=b)
        res$IMD[res$alpha==a&res$beta==b]=IMD(x,alpha=a,beta=b)
      }
    }
    ress=rbind(ress,res)
  }
}


# plot them
colz=rev(colorRampPalette(c('#f9320c','#f9c00c', '#00b9f1'))(255))

par(mfrow=c(5,5))
for(i in 1:5){
  for(j in 1:5){
    res = subset(ress,sd==vars[i]&range==ranges[j])
    rIID=rasterFromXYZ(cbind(res$beta,res$alpha,res$IID))
    #rIMD=rasterFromXYZ(cbind(res$beta,res$alpha,res$IMD))
    plot(rIID,col=colz,legend=F,xaxt='n', yaxt='n',zlim=c(min(ress$IID),max(ress$IID))) 

    #  plot(rIMD,col=colz,legend=F,xaxt='n', yaxt='n',zlim=c(min(ress$IMD),max(ress$IMD))) 
      
  }
}

dev.off()
plot(rIMD,col=colz,legend=T,xaxt='n', yaxt='n',zlim=c(min(ress$IMD),max(ress$IMD))) 





# Simulate 20 countries with the same SD and Range
#  estimated u5m probability surface (logit space)
#  estimated u5m probability surface (logit space)
set.seed(12344)
r=calc(makeRandomSurface(sd=2,l=l,scale=.2,offset=-5), fun = inv.logit)
st=stack(r)

for(i in 1:9)
  st= addLayer(st,calc(makeRandomSurface(sd=2,l=l,scale=.2,offset=-5), fun = inv.logit))
names(st)=c('Squareistan',paste('Country',2:10))


plot(st,
     col=colz,zlim=c(0,0.4),
     legend=F,
     xaxt='n',
     yaxt='n', nc = 4, nr = 2) 


# Gini
res=data.table(country=1:10,gini=NA)

for(i in 1:10)
 res$gini[res$country==i]=IID(x=as.vector(st[[i]]),alpha=1,beta=1)



#### Gini as a function of SD
sds = seq(0,5,.25)
set.seed(12345)
res=data.table(sd=rep(sds,10),c=rep(1:10,length(sds)),gini=NA)
for(i in 1:length(sds)){
  message(sds[i])
  for(c in 1:10){
    message(c)
    res$gini[res$sd==sds[i]&res$c==c]=
      IID(x=as.vector(calc(makeRandomSurface(sd=sds[i],l=l,scale=.2,offset=-5), fun = inv.logit)),
          alpha=1,beta=1)
  }
}

# plot with line and bars
ress=res[,list(mean=mean(gini),
               max=max(gini),
               min=min(gini)),
         by=sd]

require(ggplot2)
ggplot(ress,aes(x=sd,y=mean,ymin=min,ymax=max))+
  geom_point()+
  geom_errorbar()+ylab('gini')+xlab('SD')+
  theme_bw()




#########
# which inequality metrics are least sensititive to model specification?
## SD

sds = seq(0,5,.25)
set.seed(12345)
res=data.table(sd=rep(sds,10),c=rep(1:10,length(sds)),
               gini=NA,cv=NA,madmed=NA,range=NA,rangeratio=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA)

for(i in 1:length(sds)){
  message(sds[i])
  for(c in 1:10){
    message(c)
    x=as.vector(calc(makeRandomSurface(sd=sds[i],l=l,scale=.2,offset=-5), fun = inv.logit))
                     
    res$gini[res$sd==sds[i]&res$c==c]           =  
      IID(x,alpha=1,beta=1)
    res$madmed[res$sd==sds[i]&res$c==c]         =     
      madmed(x,weights=rep(1,length(x)))
    res$range[res$sd==sds[i]&res$c==c]          =  
      max(x)-mean(x)
    res$rangeratio[res$sd==sds[i]&res$c==c]     =   
      max(x)/mean(x)
    res$range1090[res$sd==sds[i]&res$c==c]      = 
      unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
    res$rangeratio1090[res$sd==sds[i]&res$c==c] =  
      unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
    res$variance[res$sd==sds[i]&res$c==c]       =      
      var(x)
    res$MAD[res$sd==sds[i]&res$c==c]            =  
      mad(x)
    res$cv[res$sd==sds[i]&res$c==c]            =  
      mean(x)/max(x)
      
  }
}

ress=res[, lapply(.SD, mean, na.rm=TRUE), by=sd ]
ress=reshape(ress,idvar='sd',direction='long',
             varying=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
             times=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
             v.names='metric')

ggplot(ress,aes(x=sd,y=metric))+
  geom_line()+theme_bw()+
  facet_wrap(~time,scales='free')


########################
## RANGE
ranges = c(0.01, 0.05 ,seq(.1,1,.1))
set.seed(12345)
res=data.table(sd=rep(ranges,10),c=rep(1:10,length(ranges)),
               gini=NA,cv=NA,madmed=NA,range=NA,rangeratio=NA,range1090=NA,rangeratio1090=NA,
               variance=NA,MAD=NA)

for(i in 1:length(ranges)){
  message(ranges[i])
  for(c in 1:10){
    message(c)
    x=as.vector(calc(makeRandomSurface(sd=2,l=l,scale=ranges[i],offset=-5), fun = inv.logit))
    
    res$gini[res$sd==ranges[i]&res$c==c]           =  
      IID(x,alpha=1,beta=1)
    res$madmed[res$sd==ranges[i]&res$c==c]         =     
      madmed(x,weights=rep(1,length(x)))
    res$range[res$sd==ranges[i]&res$c==c]          =  
      max(x)-mean(x)
    res$rangeratio[res$sd==ranges[i]&res$c==c]     =   
      max(x)/mean(x)
    res$range1090[res$sd==ranges[i]&res$c==c]      = 
      unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
    res$rangeratio1090[res$sd==ranges[i]&res$c==c] =  
      unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
    res$variance[res$sd==ranges[i]&res$c==c]       =      
      var(x)
    res$MAD[res$sd==ranges[i]&res$c==c]            =  
      mad(x)
    res$cv[res$sd==ranges[i]&res$c==c]            =  
      mean(x)/max(x)
    
  }
}

ress=res[, lapply(.SD, mean, na.rm=TRUE), by=sd ]
ress=reshape(ress,idvar='sd',direction='long',
             varying=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
             times=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
             v.names='metric')

ggplot(ress,aes(x=sd,y=metric))+
  geom_line()+theme_bw()+xlab('RANGE')+
  facet_wrap(~time,scales='free')

  
  
  
  ######################################
  ## INTERCEPT
  sds = seq(-10,0,1)
  set.seed(12345)
  res=data.table(sd=rep(sds,10),c=rep(1:10,length(sds)),
                 gini=NA,cv=NA,madmed=NA,range=NA,rangeratio=NA,range1090=NA,rangeratio1090=NA,
                 variance=NA,MAD=NA)
  
  for(i in 1:length(sds)){
    message(sds[i])
    for(c in 1:10){
      message(c)
      x=as.vector(calc(makeRandomSurface(sd=2,l=l,scale=.2,offset=sds[i]), fun = inv.logit))
      
      res$gini[res$sd==sds[i]&res$c==c]           =  
        IID(x,alpha=1,beta=1)
      res$madmed[res$sd==sds[i]&res$c==c]         =     
        madmed(x,weights=rep(1,length(x)))
      res$range[res$sd==sds[i]&res$c==c]          =  
        max(x)-mean(x)
      res$rangeratio[res$sd==sds[i]&res$c==c]     =   
        max(x)/mean(x)
      res$range1090[res$sd==sds[i]&res$c==c]      = 
        unname(quantile(x,p=.9,na.rm=T)-quantile(x,p=.1,na.rm=T))
      res$rangeratio1090[res$sd==sds[i]&res$c==c] =  
        unname(quantile(x,p=.9,na.rm=T)/quantile(x,p=.1,na.rm=T))
      res$variance[res$sd==sds[i]&res$c==c]       =      
        var(x)
      res$MAD[res$sd==sds[i]&res$c==c]            =  
        mad(x)
      res$cv[res$sd==sds[i]&res$c==c]            =  
        mean(x)/max(x)
      
    }
  }
  
  ress=res[, lapply(.SD, mean, na.rm=TRUE), by=sd ]
  ress=reshape(ress,idvar='sd',direction='long',
               varying=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
               times=c('MAD','gini','cv','madmed','range','rangeratio','range1090','rangeratio1090','variance'),
               v.names='metric')
  
  ggplot(ress,aes(x=sd,y=metric))+
    geom_line()+theme_bw()+xlab('Intercept')+
  facet_wrap(~time,scales='free')
  
  
  
  
  
  
  
  
###########################
### Does our model attenuate GINI over time? 
  ## compare data gini with surface gini
  
### simulate observations over 3 periods
set.seed(12344)
  r2 = stack(makeRandomSurface(sd=2,l=l,scale=.2,offset=-5))
  r2 = stack(r2,0.5*r2[[1]]+makeRandomSurface(sd=0.5,l=l,scale=.2,offset=-2))
  r2 = stack(r2,0.5*r2[[2]]+makeRandomSurface(sd=0.2,l=l,scale=.2,offset=-3))
  
  r2 = calc(r2, fun = inv.logit)
  plot(r2)
  
  coo = cbind(x=runif(100,0,1),y=runif(100,0,1))
  d=data.table(x=rep(coo[,1],3),
               y=rep(coo[,2],3),
               per=rep(c('1','2','3'),each=100),
               died=NA)
  d$died[d$per=='1']=rbinom(n=100,prob=extract(r2[[1]],coo),size=1000)
  d$died[d$per=='2']=rbinom(n=100,prob=extract(r2[[2]],coo),size=1000)
  d$died[d$per=='3']=rbinom(n=100,prob=extract(r2[[3]],coo),size=1000)
               
  # plot simulated data
  ggplot(d,aes(x=x,y=y,colour=per))+
    scale_colour_manual(values=c('#E71D36','#2EC4B6', '#011638'))+
    geom_jitter(aes(size=died),shape=21,stroke=1.7,alpha=.5)+ #,width=.01,height=.01)+
    theme_bw()+ 
    scale_size(range = c(0, 20))

  par(mfrow=c(3,1))
  hist(d$died[d$per==1],col='#E71D36',main="period 1",xlab='# died',xlim=c(0,175),breaks=20)
  hist(d$died[d$per==2],col='#2EC4B6',main="period 2",xlab='# died',xlim=c(0,175),breaks=20)
  hist(d$died[d$per==3],col='#011638',main="period 3",xlab='# died',xlim=c(0,175),breaks=20)
  
  IID(d$died[d$per==1])
  IID(d$died[d$per==2])
  IID(d$died[d$per==3])
  
  
  fit=fittmb2(dt=d)
  plot(fit$mean_ras)
  
  
  
  