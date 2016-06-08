
# set up the true parameters
# TODO:
# set up full framework for common comparisons
# run many models on the cluster. 
# run each 100 times at each level, and plot the mean value
# wrapper function for running simulation, fitting models, and returning results

# things to modulate about the simulation
# number of covariates
# number of periods
# relative predictive strength of covariates to latent field
# probability of death
# number of nodes
# mean exposures per cluster

# things to graph as a function of the above, comparing inla to tmb
# runtime
# parameter error
# mean predictive error
# root mean squared error
# objective (max L )


### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit Model using INLA
rm(list=ls())

require(boot)
library(doParallel)
require(snow)
library(data.table)

#root=paste0(ifelse(Sys.info()['sysname']=='Windows','J:','/home/j'),
#            '/temp/geospatial/tmb_inla_sim')
root="C:/Users/royburst/Google Drive/spatial_ecology/hw/project/"

setwd(root)
source("code/utils.r")



# named list of vars to change and their values
todolist = list(#'exposures'=seq(50,1000,25),
  #'numperiods'=2:10,
  #'nc'=1:10,
  'numclusters'=seq(50,1000,50)) #,
#'betavalues'=seq(-3,3,.1),
#'int'=seq(-3,3,.1),
#'ratio'=seq(-2,0,.05))


sims=3

# loop through list else do this
# maybe parallelize it as well
for(j in 1:length(todolist)){
  #for(j in 2:2) {
  
  var=names(todolist)[j]
  values=todolist[[var]]
  message(var)
  
  cl <- snow::makeCluster(ifelse(Sys.info()['sysname']=='Windows',4,40)) # use 60 threads if on cluster
  print(cl)
  
  snow::clusterExport(cl,fcnlist)
  
  pvmetrics=data.table() # start results table
  for(simulation in 1:sims){
    message(simulation)
    
    # in parallel, run full simulation for the changing params of interest
    if(var=='nc')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(nc=x))   ,TRUE)
    if(var=='exposures')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(exposures=x))   ,TRUE)
    if(var=='numclusters')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(numclusters=x))   ,TRUE)
    if(var=='numperiods')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(numperiods=x))   ,TRUE)
    if(var=='rho')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(rho=x))   ,TRUE)
    if(var=='int')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(betas=c(x,-.01)))     ,TRUE) 
    if(var=='betavalues')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(betas=c(-.1,x)))       ,TRUE)
    if(var=='ratio')
      x = try(snow::parLapply(cl=cl, x= values,function(x) simNfit(betas=c(x,x)))       ,TRUE)
    
    
    # clean up reporting table ( if x didnt crash)
    if(!inherits(x, "try-error")){
      
      for(i in 1:length(values)){
        tmp=data.frame(x[[i]]$pv)
        tmp[,'metric']=rownames(tmp)
        tmp$value=values[i]
        if(var=="ratio")
          tmp$value=as.numeric(round(as.numeric(x[[i]]$l.s.ratio),2))
        tmp=data.table(tmp)
        
        pvmetrics=rbind(pvmetrics,
                        tmp,
                        data.table(cbind(x[[i]]$runtimes,metric='runtime',value=values[i])),
                        data.table(cbind(x[[i]]$ll,metric='ll',value=values[i])))
        
        # param differences
        tmp=x[[i]]$params$diff
        tmp$metric <- row.names(tmp)  
        tmp$metric[grep('X',tmp$metric)]='X'
        
        tmp$value=values[i]
        if(var=="ratio")
          tmp$value=as.numeric(round(as.numeric(x[[i]]$l.s.ratio),1))
        tmp=data.table(tmp)
        
        pvmetrics=rbind(pvmetrics,tmp)
        
        
        
      } 
    }
    
  }
  
  
  snow::stopCluster(cl)
  
  pvmetrics$sim=1
  pvmetrics = pvmetrics[,list(inla=mean(as.numeric(inla)),
                              tmb =mean(as.numeric(tmb)),
                              sims=sum(as.numeric(sim))),
                        by=list(value,metric)]
  
  # clean up pesky inla.model files that get created. 
  # DO NO DO OVER VPN. CLUSTER IS BEST
  setwd(root)
  unlink(dir(c(getwd(),paste0(getwd(),'/code/')),pattern="inla.model",full.names=T),recursive=T,force=T)
  
  # save a CSV for graphing
  write.csv(pvmetrics,
            paste0('output/data/',var,'.csv'))
  
  
}

#










if(TRUE==FALSE)
  
{
  colz=rev(colorRampPalette(c('#E71D36', '#2EC4B6', '#EFFFE9'))(255))
  
  # sim examples
  sim=mortsim(betas=c(0,0),
              nc=1,
              rho=.85,
              n_periods = 4,
              n_clusters = 30)
  raster::plot(sim$r.true.mr,col=colz,legend=F,xaxt='n',yaxt='n')
  
  sim=mortsim(betas=c(-.5,-3,-1),
              nc=2,
              rho=.85,
              n_periods = 6)
  raster::plot(sim$r.true.mr,col=colz,legend=F,xaxt='n',yaxt='n')
  
  
  # simulate clusters
  par(pch=1)
  plot(sim$r.true.mr[[5]],col=(colz),legend=F,xaxt='n',yaxt='n')
  points(sim$d$x[sim$d$period==5],
         sim$d$y[sim$d$period==5],
         cex=-30*log10(sim$d$mr[sim$d$period==5]))
  
  # true fit comparison
  
  sim=mortsim(betas=c(-3,-3,-1),nc=2,mean.exposure.months = 100,n_clusters=50)
  mod.inla=fitinla(simobj=sim)
  mod.tmb= fittmb(simobj=sim)
  
  p=1
  z=as.vector(cbind(as.matrix(sim$r.true.mr[[p]]),
                    as.matrix(mod.inla$mean_ras[[p]]),
                    as.matrix(mod.tmb$mean_ras[[p]])))
  mn=min(z);mx=max(z)
  
  
  setwd(root)
  pdf('output/sandbox/temp3.pdf',width=9,height=10)
  par(mfrow=c(2,2))
  plot(sim$r.true.mr[[p]]     ,col=colz,xaxt='n',yaxt='n', main='TRUTH',zlim=c(mn,mx))
  points(sim$d$x[sim$d$period==p],
         sim$d$y[sim$d$period==p],
         cex=-1*log10(sim$d$mr[sim$d$period==p]))
  plot(mod.inla$mean_ras[[p]],col=colz,xaxt='n',yaxt='n',main='MEAN ESTIMATE (INLA)',zlim=c(mn,mx))
  plot(mod.tmb$mean_ras[[p]] ,col=colz,xaxt='n',yaxt='n',main='MEAN ESTIMATE (TMB)',zlim=c(mn,mx))
  dev.off()
  
  
  
  # compare outputs
  # make reporting plots for this.. 
  p=1
  
  z=as.vector(cbind(as.matrix(sim$r.true.mr[[p]]),
                    as.matrix(mod.inla$mean_ras[[p]]),
                    as.matrix(mod.tmb$mean_ras[[p]])))
  mn=min(z);mx=max(z)
  
  par(mfrow=c(2,2))
  plot(sim$r.true.mr[[p]]    ,zlim=c(mn,mx),main='TRUTH')
  plot(mod.inla$mean_ras[[p]],zlim=c(mn,mx),main='MEAN ESTIMATE (INLA)')
  plot(mod.tmb$mean_ras[[p]] ,zlim=c(mn,mx),main='MEAN ESTIMATE (TMB)')
  
  
  # print runtimes
  print(mod.inla$model.runtime)
  print(mod.tmb$model.runtime)
  
  # print covariates (with their coefficients)
  plot(sim$cov.raster)
  
  
  # compare epsilons
  plot(sim$S.raster,zlim=c(min(as.matrix(sim$S.raster)),max(as.matrix(sim$S.raster))),main='TRUE Latent Field')
  plot(mod.tmb$s.raster,zlim=c(min(as.matrix(mod.tmb$s.raster)),max(as.matrix(mod.tmb$s.raster))),main='TMB Latent Field')
  plot(mod.inla$s.raster,zlim=c(min(as.matrix(mod.inla$s.raster)),max(as.matrix(mod.inla$s.raster))),main='INLA Latent Field')
  # tmb is systematically underestimating the RF
  
  
  
  ###
  # Function to return graphs with model object inputs, options to choose which graphs get returned
  
  
}















