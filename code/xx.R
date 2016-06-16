# Use draw level raked outputs to estimate some inequality metrics:
# RGE, CV, Gini, Theil


## ~~~~~~~~
rm(list=ls())

library(raster)
library(seegMBG)
library(seegSDM)

source('code/functions.R')


## ~~~~~~~
# bring in data we need ( shoul have already run combine_predictions.R)
message('Loading in data')

# raw predictions
load('output/child_raked_raw_predictions.RData')
#load('output/infant_raked_raw_predictions.RData')
#load('output/neonatal_raked_raw_predictions.RData')

# ad0 mapping
ad0 <-raster(paste0('data/clean/shapefiles/ad0_raster.grd'))




## ~~~~~~~~
# write functions to pull various inequality metrics
message('Writing and loading inequality functions')

# abs range over mean (lowlimit and highlimit are quantiles)

# range ratio (lowlimit and highlimit are quantiles)
# make sure to only give it places with a population
RGR<-function(x,weights=rep(1, length = length(x))){
  weights=rep(1, length = length(x)) # no actual weights, just tricking condSim
  unname(quantile(x*weights,p=.9,na.rm=T)/quantile(x*weights,p=.1,na.rm=T))
}


# Coefficient of Variation already exists as cv()

# gini exists already as gini() from reldist:
gini=function (x, weights = rep(1, length = length(x))) 
{
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

# Theil's T statistic: come back to this




## ~~~~~~~~~~~~~~
# run condSim() with our functions of inequality
name='child'

name_to_gini  <- sprintf('%s_raked_gini_sim', name)
name_to_rge   <- sprintf('%s_raked_rge_sim', name)
name_to_rgr   <- sprintf('%s_raked_rgr_sim', name)
name_to_cv    <- sprintf('%s_raked_cv_sim', name)
name_to_theil <- sprintf('%s_raked_theil_sim', name)


# get object (draws)
x <- get('raked_5q0')

# get NA mask
good_cells <- which(!is.na(x[, 1]))


message('gini')
gini_sim <- condSim(x[good_cells, ],
                    weights = pop_wt_all_ad0[good_cells],
                    group = ad0_code_all[good_cells],
                    fun = gini)

message('range (10% and 90%) ratio')

# maybe instead of weight, look only in populated cells..
rgr_sim<- condSim(x[good_cells, ][pop_cell>10,],
                  weights = NULL,
                  group   = ad0_code_all[good_cells][pop_cell>10],
                  fun     = RGR)

# assign summary stats to object

assign(name_to_gini, gini_sim)
assign(name_to_rgr , rgr_sim)

  


 


## ~~~~~~~~~~~~~~
# Visualize the metrics

