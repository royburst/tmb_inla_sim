############
# look at different ineq metrics on real ad0 africa data..




################
# Load data
  
  setwd('J:/temp/geospatial/U5M_africa')
  
  # load miscellaneous project functions
  source('code/functions.R')
  
  # start logging the run
  #startLog() 
  
  # set global variables for final model run
  defaultOptions(resolution = 5,        # raster resolution
                 country = '',       # MOZ. if there is a 3 letter ISO here, the analysis will be run only on that country
                 location = 'seattle',  # location for final run 
                 cores = 30,            # number of cores to use
                 full = FALSE,           # whether to do a full run
                 age_bins = c(1,2,3,4,6),        # which age bins to run models for
                 spacetime = TRUE,      # whether to do space-time inference
                 sbh_variance = TRUE,   # whether to incorporate variance in SBH estimate
                 champs_experiment = FALSE, # If true, will save outputs elsewhere
                 ce_holdout = FALSE,    # Only applicable if ce is true, will hold out data
                 validation = F,    # if true, will run 5-fold cross-validation instead of full moel
                 start = Sys.time())    # start time
  
  
  # data
  load('output/child_raked_raw_predictions.RData')
  
  # country IDs and population weights
  
    # get samples of national, population-weighted mortality rates in each fraction.
    
    # admin codes, matched up with names
    # template raster
    mask <- brick('data/clean/covs_transformed.grd') * 0
    pop<-brick('data/raw/covariates/new_20160421/pop_stack.tif')
    cell_idx <- cellIdx(mask)
    
    ad0 <- raster(paste0('data/clean/shapefiles/ad0_raster',getOption('country'),'.grd'))
    ad0_cell <- extract(ad0, cell_idx)
    sad0 <- shapefile(paste0('data/clean/shapefiles/africa_ad0.shp'))@data[c('name','gaul_code')]
    sad0$name[sad0$gaul_code==66]="Cote dIvoire"
    ad0_code <- match(ad0_cell,sad0$gaul_code)
    ad0_code <- sad0$name[ad0_code]
    # cell populations
    pop_cell <- extract(pop[[2]], cell_idx) # using 2005 for now, later split these up
  
    # set to zero if the population or admin code is NA
    pop_cell[is.na(pop_cell)] <- 0
    pop_cell[is.na(ad0_code)] <- 0  #
  
    # changed to cell by RB (6APR2016) to make sure 0 pop adm2s get in there
    # get population totals in all cells
    pop_totals_ad0 <- tapply(pop_cell, ad0_code, sum)
  
    # find those with 0 population andremove them from both totals and ad0_code
    rem_ad0 <- names(pop_totals_ad0)[pop_totals_ad0 == 0]
    ad0_code[ad0_code %in% rem_ad0] <- NA
    pop_totals_ad0 <- pop_totals_ad0[pop_totals_ad0 > 0]
    pop_totals_ad0_cell <- as.vector(pop_totals_ad0)[match(ad0_code,
                                                           names(pop_totals_ad0))]
    # # for all zero population, set to 1 to avoid divide-by zero errors
    pop_totals_ad0_cell[pop_totals_ad0_cell == 0] <- 1
    # for all cells with some population, but unknown country, set to 0/1
    disowned_ad0 <- which(is.na(pop_totals_ad0_cell))
    pop_totals_ad0_cell[disowned_ad0] <- 1
    pop_cell[disowned_ad0] <- 0
  
    # get national population weights for each cell
    pop_wt_ad0 <- pop_cell / pop_totals_ad0_cell
    
    # make sure these sum to one or zero
    wt_sum_ad0 <- tapply(pop_wt_ad0, ad0_code, sum)
    stopifnot(all.equal(wt_sum_ad0, round(wt_sum_ad0)))
  
    
    # replicate for multiple years 
    prs=4
    pop_wt_all_ad0 <- rep(pop_wt_ad0, prs)
    pop_cell_all <- rep(pop_cell, prs)
    periods <- rep(c(2000, 2005, 2010, 2015), each = length(ad0_code))
    periods[is.na(ad0_code)] <- ''
    ad0_code_all <- paste(ad0_code, periods, sep = '_')
    ad0_code_all[ad0_code_all == 'NA'] <- NA
  
    
    
    
  # calculate inequality metrics:
    
    root="C:/Users/royburst/Documents/tmb_inla_sim"
    setwd(root)
    source('code/utils.R')
    
    x = cbind(rowMeans(raked_5q0))
    
    good_cells <- which(!is.na(x[, 1]))

  
    # CALCULATE INEQUALITY MEASURES
    res = data.frame(
      gini_popwt_ad0 = condSim(vals=    cbind(x[good_cells, ]),
                               weights = pop_cell_all[good_cells],
                               group   = ad0_code_all[good_cells],
                               fun     = gini),

      RGR_popwt_sim = condSim(vals=     cbind(x[good_cells, ]),
                              weights = rep(pop_cell,4)[good_cells],
                              group =   ad0_code_all[good_cells],
                              fun     = RGR,hi=.9,lo=.1),
      
      range_popwt_sim = condSim(vals=    cbind(x[good_cells, ]),
                                weights = rep(pop_cell,4)[good_cells],
                                group = ad0_code_all[good_cells],
                                fun     = range,hi=.9,lo=.1),
      
      cv_popwt_sim = condSim(vals=    cbind(x[good_cells, ]),
                             weights = rep(pop_cell,4)[good_cells],
                             group = ad0_code_all[good_cells],
                              fun     = coefvar)# ,
     
     # madmed_popwt_sim = condSim(vals    = cbind(x[good_cells, ]),
     #                            weights = pop_cell[good_cells],
     #                            group   = ad0_code_all[good_cells],
     #                            fun     = madmed)
    )
    
  # plot them

    plot(res)
  
    
    
    
    
    
    
    
    #
    setwd('J:/temp/geospatial/U5M_africa')
    

    
    # load packages
    library(RColorBrewer)
    library(plyr)
    library(ggplot2)
    library(gridExtra)
    library(raster)
    # source misc functions
    source('code/functions.R')
    library(seegMBG)
    # ~~~~~~~~~~~
    # load data
    
    # iso3 lookup tables
    lookup <- read.csv('data/raw/gbd2015_national_estimates/iso3_lookup.csv',
                       stringsAsFactors = FALSE)
    lookup2 <- read.csv('data/clean/ad0_iso_code_lookup.csv',
                        stringsAsFactors = FALSE)
    
   
    # save csvs of adm 0, 1 , 2

    gini_5q0=data.frame(iso3=splitGeoNames(res)$iso3,
                        year=splitGeoNames(res)$year,
                        med=splitGeoNames(res)$gini_popwt_ad0)
    RGR_5q0=data.frame(iso3=splitGeoNames(res)$iso3,
                        year=splitGeoNames(res)$year,
                        med=splitGeoNames(res)$RGR_popwt_sim)
    RNG_5q0=data.frame(iso3=splitGeoNames(res)$iso3,
                        year=splitGeoNames(res)$year,
                        med=splitGeoNames(res)$range_popwt_sim)
    CV_5q0=data.frame(iso3=splitGeoNames(res)$iso3,
                        year=splitGeoNames(res)$year,
                        med=splitGeoNames(res)$cv_popwt_sim)
    

    years <- c(2000, 2005, 2010, 2015)
    pop <- brick('data/raw/covariates/new_20160421/pop_stack.tif')[[2]]
    cell_idx <- cellIdx(pop)
    ad0 <- raster(paste0('data/clean/shapefiles/ad0_raster',getOption('country'),'.grd'))
    ad0_cell <- extract(ad0, cell_idx)
    sad0 <- shapefile(paste0('data/clean/shapefiles/africa_ad0.shp'))@data[c('name','gaul_code')]
    sad0$name[sad0$gaul_code==66]="Cote dIvoire"
    ad0_code <- match(ad0_cell,sad0$gaul_code)
    ad0_code <- sad0$name[ad0_code]
    pops<-as.vector(pop)[cell_idx]
    pops<-aggregate(pops~ad0_code,FUN=sum)
    
    
    
    addregionalaverages<-function(d,p,ci=F){
      if(ci==F)d$upper=d$lower=NULL
      
      
      store=d
      d<-merge(d,p,by.x='iso3',by.y='ad0_code',all.x=T)
      d$region=NA
      d$region[d$iso3%in%c('Tunisia','Egypt','Morocco','Sudan')]='North'
      d$region[d$iso3%in%c('Burundi', 'Comoros', 'Djibouti', 'Eritrea' ,
                           'Ethiopia', 'Kenya', 'Madagascar', 'Malawi' ,
                           'Mozambique', 'Rwanda', 'Somalia', 'United Republic of Tanzania' ,
                           'Uganda', 'Zambia', 'South Sudan')]="East"
      d$region[d$iso3%in%c('Botswana' ,'Lesotho' ,'Namibia' ,'South Africa','Swaziland' ,'Zimbabwe')]="South" 
      d$region[d$iso3%in%c('Benin','Burkina Faso','Cameroon' ,'Cape Verde' ,'Chad',
                           'Cote dIvoire','Gambia','Ghana','Guinea','Guinea-Bissau',
                           'Liberia','Mali','Mauritania','Niger','Nigeria','Sao Tome & Principe',
                           'Senegal','Sierra Leone', 'Togo')]='West'
      d$region[d$iso3%in%c('Angola','Central African Republic','Congo','Democratic Republic of the Congo',
                           'Equatorial Guinea','Gabon')]='Central'
      
      # population weights
      rp<-aggregate(pops~region+year,data=d,FUN=sum)
      d<-merge(d,rp,by=c('region','year'))        
      d$weight = d$pops.x/d$pops.y
      d$med = d$med*d$weight
      
      if(ci){
        d$upper = d$upper*d$weight
        d$lower = d$lower*d$weight
        
        # final average
        tmp<-aggregate(cbind(med,upper,lower)~region+year,data=d,FUN=sum)
      } else {
        tmp<-aggregate(cbind(med)~region+year,data=d,FUN=sum)
        
        
        
      }
      
      
      tmp$iso3=tmp$region;tmp$region=NULL
      
      return(rbind(store,tmp))
    }
    
    
    gini_5q0=addregionalaverages(gini_5q0,p=pops)
    RGR_5q0=addregionalaverages(RGR_5q0,p=pops)
    RNG_5q0=addregionalaverages(RNG_5q0,p=pops)
    CV_5q0=addregionalaverages(CV_5q0,p=pops)
    
    
    
    
    # keep only requid countries
    countries <- sort(unique(c(gini_5q0$iso3)))
    
    countries = countries[!countries%in%c('Algeria','Libya','Tunisia','Western Sahara',"Ma'tan al-Sarra")]
    
    n_ctry <- length(countries)
    
    set <- brewer.pal(8, 'Dark2')

    gini_5q0$yearnum <- match(gini_5q0$year, years)
    RGR_5q0$yearnum <- match(RGR_5q0$year, years)
    RNG_5q0$yearnum <- match(RNG_5q0$year, years)
    CV_5q0$yearnum <- match(CV_5q0$year, years)
    
    
    ############################################################################################################
    #nick gs gini plot
    
    

    ctry_highlight <- ''
    regions <- c('North', 'West', 'East', 'Central', 'South')
    
    # sort countries
    countries <- countries[!(countries %in% regions)]
    
    n_ctry <- length(countries)
    
    # get distinctive colours for each
    set <- brewer.pal(8, 'Dark2')[-c(5, 7)]
    
   
    
    
    
    # ~~~~~~~~~~~
    # plot population-weighted Gini coefficients for each country
    
    # define countries to highlight
    region_col <- set[1:length(regions)]
    
    
    
    
    
    ## ~~~~
    # RB remake of gini plot using ggplot

    for (rate in c('gini_5q0', 'RGR_5q0', 'RNG_5q0','CV_5q0')) {
      
      # get human-friendly name
      name <- switch(rate,
                     `gini_5q0` = 'GINI',
                     `RGR_5q0` = '9010 Ratio',
                     `RNG_5q0` = '9010 Diff',
                     `CV_5q0` = 'Coef. of Variation')
      
      pdf(sprintf('./output/figures/mortality_gini_%s4.pdf',
                  name),
          width = 14,
          height = 8,
          pointsize = 18)
      
          gini=get(rate)
          
            
          gini$Region = gini$iso3
  
          gini$year=as.numeric(gini$year)
  
      gini.sub  <- subset(gini, iso3%in%c('North','East','South','Central','West')) 
      gini.cnt  <- subset(gini,!iso3%in%c('North','East','South','Central','West')) 
      gini.dum = gini.sub
      gini.dum$med=NaN
      cols=c('#f9c00c', '#00b9f1', '#7200da', '#f9320c', '#090707')
      
      plot(
      ggplot(data=gini,aes(x=year,y=med,group=iso3))+
        geom_path(data=gini.cnt,size=2,alpha=.4,colour='grey',lineend="round")+
        geom_path(data=gini.sub,size=4,alpha=.7,aes(colour=Region),lineend="round")+
        geom_path(data=gini.dum,size=4,alpha=1,aes(colour=Region),lineend="round",na.rm=T)+
        guides(colour = guide_legend(reverse=TRUE))+
        theme_bw()+
        theme(axis.line = element_line(colour = "black"),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              legend.key = element_blank())+
        scale_colour_manual(values=cols)+
        ylab(name)
      )
      dev.off()
    
    }
    
    
    
    
    
    
    
    
 ### Check out random cross-sections of africa.. 
    ad0 <- raster(paste0('data/clean/shapefiles/ad0_raster',getOption('country'),'.grd'))
    r=raked_5q0[1:(nrow(raked_5q0)/4),]
    
    x=insertRaster(mask,cbind(r[,1:20]))
    x=raster::as.array(x)
    m=rowMeans(cbind(x[700:900,1000,1],x[700:900,1000,10],
                     x[700:900,1000,2],x[700:900,1000,11],
                     x[700:900,1000,3],x[700:900,1000,12],
                     x[700:900,1000,4],x[700:900,1000,13],
                     x[700:900,1000,5],x[700:900,1000,14],
                     x[700:900,1000,6],x[700:900,1000,15],
                     x[700:900,1000,7],x[700:900,1000,16],
                     x[700:900,1000,8],x[700:900,1000,17],
                     x[700:900,1000,9],x[700:900,1000,18],x[700:900,1000,19],x[700:900,1000,20]))
    require(scales)
    plot(x[700:900,1000,1],type='l',col=alpha("black",.5),ylim=c(.1,.25))
    for(i in 2:20)
    lines(x[700:900,1000,i],type='l',col=alpha("black",.5))
    lines(m,type='l',col='red',lwd=5)
    

    

    
    
    
    
    
    