 # Test Inequality Metrics on Fake data



rm(list=ls())

root="C:/Users/royburst/Google Drive/spatial_ecology/hw/project/"
setwd(root)
source('code/utils.R')

library(seegMBG)
library(boot)
library(raster)

#  Simulate surface to test this stuff on
r = makeRandomCovariate(sd=10,l=100,ext=F)

#  Simulate population raster
p = makeRandomCovariate(sd=2,l=100,ext=F)
p = calc(p, fun = exp)*10
p = insertRaster(p,data.frame(rpois(100*100,(as.vector(p)))))
names(p)=population


#  Simulate country ID raster surface
cntrys = crop(resample(raster(outer(seq(1,2, l = 2), 
                      seq(3,4, l = 2))),r,method='ngb'),r)   
                         





