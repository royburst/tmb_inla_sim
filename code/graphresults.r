# bring in results table, plot comparison over time


# setup
require(ggplot2); theme_set(theme_bw())
#root=paste0(ifelse(Sys.info()['sysname']=='Windows','J:','/home/j'),
#            '/temp/geospatial/tmb_inla_sim')
root="C:/Users/royburst/Google Drive/spatial_ecology/hw/project/"
setwd(root)
source("code/utils.r")


# search for results and plot them
xs=gsub('.csv','',dir('output/data'))

for(x in xs){
  message(x)
  
  d=fread(paste0('output/data/',x,'.csv'))
  d$value = as.numeric(d$value)
  
  ys=unique(d$metric)
  sims=min(d$sims)
  
  if(!dir.exists(paste0('output/figures/',x)))
    dir.create(paste0('output/figures/',x))
  
  for(y in ys){
    message(y)
    png(paste0('output/figures/',x,'/',y,'.png'),width=1000,height=1000)
    
    print(ggplot(d[d$metric==y,],aes(x=value,y=inla))+
            geom_line(color='red')+
            geom_line(aes(y=tmb),colour='blue') +
            ylab(y)+
            xlab(x)+
            ggtitle(paste0(sims,' Simulations per level. (TMB=blue,INLA=red)')))
    
    dev.off()
  }
}

