# RB

# 1. Simulate a surface.
# 2. Simulate survey points#
# 3. Resample points in polygons and then back to points
# 4. Return prediction using INLA for point data with (UC)
# 5. Return prediction surface for polygon data (with UC)
# 6. Look at differences in UC and predictive validity

rm(list=ls())

root="C:/Users/royburst/Documents/code/tmb_inla_sim"
setwd(root)
source('code/utils.R')
colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))


library(seegMBG)
library(boot)
library(raster)
library(ggplot2)
library(sp)

###################################################################
###################################################################
## 1 and 2
# data simulation
simobj= mortsim(nu         = 2  ,
                betas      = c(-5,0)      ,
                scale      = .2            ,
                Sigma2     = c(3,1,.5,.25),
                rho        = 0.9          ,
                l          = 151           ,
                n_clusters = 75          ,
                n_periods  = 4            ,
                mean.exposure.months = 100,
                extent = c(0,1,0,1)       ,
                ncovariates = 1           ,
                seed   = 1234455           ,
                returnall=TRUE            )

# just take 1 year
plot(simobj$r.true.mr)
tr = simobj$r.true.mr[[4]] # raster of truth
d  = simobj$d[period==4,]
mn = min(as.vector(tr))
mx = max(as.vector(tr))

# plot truth
plot(tr,zlim=c(mn,mx),col=colz,legend=T,xaxt='n',yaxt='n', main='True mortality rates in an area')


# plot samples
d$observed_mr=as.numeric(d$mr)
d$children_sampled = d$exposures
ggplot(d,aes(x=x,y=y))+
  geom_jitter(aes(size=children_sampled,color=observed_mr),shape=21,stroke=5,alpha=.85)+
  scale_colour_gradient2(high='#E71D36',mid='#2EC4B6', low='#EFFFE9',limits=c(mn,mx))+
  theme_bw()+
  scale_size(range = c(0, 15)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())+
  ggtitle('Sampling Locations')

###################################################################
###################################################################
## 3
#make polygon
polycoords =  data.frame(rbind( c(0.20,0.05),
                                c(0.15,0.30),
                                c(0.19,0.90),
                                c(0.40,0.90),
                                c(0.50,0.30),

                                c(0.55,0.00),
                                c(0.35,0.00),
                                c(0.20,0.05)))
p = Polygon(polycoords)
p = Polygons(list(p),1)
p = SpatialPolygons(list(p))

# identify if points are in polygon
d$inpoly = !is.na(over(SpatialPoints(cbind(d$x,d$y)),p))

# aggregate only those points
mr=d[,list(polydeaths=sum(deaths),polyexp =sum(exposures)), by = inpoly]
mr[,mr:=polydeaths/polyexp,]
polycoords[,3]=mr$mr[mr$inpoly]
names(polycoords)=c('x','y','observed_mr')



# save d1 as a copy of d and d2 as just those not in poly
d1 = d[,c('x','y','exposures','children_sampled','deaths','inpoly','observed_mr','X1'),with=FALSE]
d1$ period = 1
d2 = d1[inpoly==FALSE,]

# resample polygons as 10 points equally weighted
samples = 20
resamp = data.table(coordinates(spsample(p,20,type='random')),
               children_sampled=mr$polyexp[mr$inpoly],
               exposures=mr$polyexp[mr$inpoly],
               deaths=mr$polydeaths[mr$inpoly],
               observed_mr = mr$mr[mr$inpoly],
               inpoly = TRUE,
               period = 1)

resamp[,children_sampled:=children_sampled*1/samples,]
resamp[,exposures:=exposures*1/samples,]
resamp[,deaths:=deaths*1/samples,]
resamp$X1=extract(simobj$cov.raster,cbind(resamp$x,resamp$y))

d2 = rbind(d2,resamp)


# plot samples with polygon
ggplot(d2,aes(x=x,y=y))+
  geom_jitter(aes(size=children_sampled,color=observed_mr),shape=21,stroke=4,alpha=.9)+
  scale_colour_gradient2(high='#E71D36',mid='#2EC4B6', low='#EFFFE9',limits=c(mn,mx))+
  theme_bw()+
  scale_size(range = c(0, 15)) +
  theme(axis.ticks = element_blank(), axis.text = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())+
  ggtitle('Sampling Locations') +
  geom_polygon(data=polycoords,aes(fill=observed_mr),alpha=.2)+
  scale_fill_gradient2(high='#E71D36',mid='#2EC4B6', low='#EFFFE9',limits=c(mn,mx),guide=FALSE)





###################################################################
###################################################################
## 4 and 5 run inla models
mod1=fitinla2(dt=d1,template=simobj$template)

mod2=fitinla2(dt=d2,template=simobj$template)


# plot truth
pdf('./plots/truthversuspoly.pdf',width=12,height=4.5)
par(mfrow=c(1,3))
plot(tr  ,zlim=c(mn,mx),col=colz,legend=T,xaxt='n',yaxt='n', main='Truth')
plot(mod1,zlim=c(mn,mx),col=colz,legend=F,xaxt='n',yaxt='n', main='All Points');lines(p)
plot(mod2,zlim=c(mn,mx),col=colz,legend=F,xaxt='n',yaxt='n', main='With Poly');lines(p)
dev.off()
