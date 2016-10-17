# May 2016


# functions
library(data.table)





## USEFUL OBJECTS
#pi
pi=3.141593

# List of all functions in here (to pass to cluster in parallelizing)
fcnlist = list('mortsim','simNfit','fitinla','fittmb','makeRandomCovariate',
               'pv','rasterFromXYZT','compareparams')





#######
## Helper function for turning an xyzt table into a raster
rasterFromXYZT <- function(table,
                           z,t){
  require(data.table)
  n_periods = length(unique(table[,t]))
  table$t=table[,t]
  res=  stack(rasterFromXYZ(as.matrix(table[t==1,c('x','y',z),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      res=addLayer(res, rasterFromXYZ(as.matrix(table[t==r,c('x','y',z),with=F])))
  return(res)
}




#######
## A function to make a random covariate
makeRandomCovariate <- function(extent = c(0,1,0,1),
                                mean   = 0,
                                sd=.1,
                                l=51,
                                scale=1,
                                ext=T){

  require(RandomFields)
  require(raster)

  if(ext) {
    extension = abs(extent[2]-extent[1])/3
  } else { extension=0 }

  cmean = 500
  csd  = 500

  while(!(cmean>(mean-1)&cmean<(mean+1))&csd>sd){
    cov.raster<- raster(outer(seq(0,1, l = l),
                              seq(0,1, l = l),
                              FUN = function(x,y) rnorm(1,0,1)*x^(rnorm(1,2,5))-
                                1*extension*(abs(y))^(runif(1,1,4)/runif(1,1,2))),
                        xmn=extent[1]-extension, xmx=extent[2]+extension,
                        ymn=extent[3]-extension, ymx=extent[4]+extension)
    cmean=(mean(as.matrix(cov.raster)))
    csd  =(sd(as.matrix(cov.raster)))

    if(cmean==Inf|csd==Inf|is.na(cmean)|is.na(csd)) {
      cmean = 500
      csd = 500
    }
  }

  RMmodel = RMmatern(nu    = 1,
                     var   = sd,
                     scale = scale)

  z=RFsimulate(model=RMmodel,
               x=coordinates(cov.raster)[,1],
               y=coordinates(cov.raster)[,2])@data[,1]

  cov.raster=cov.raster+rasterFromXYZ(data.frame(x=coordinates(cov.raster)[,1],
                                                 y=coordinates(cov.raster)[,2],
                                                 z=z))
  return(cov.raster)

}












#######
## simulate mortality data

mortsim <- function(
  nu         = 2            ,  #  Matern smoothness parameter (alpha in INLA speak)
  betas      = c(-3,-1,1,1) ,  #  Intercept coef and covariate coef For Linear predictors
  scale      = 2            ,  #  Matern scale eparameter
  Sigma2     = (1) ^ 2      ,  #  Variance (Nugget)
  rho        = 0.9          ,  #  AR1 term
  l          = 51           ,  #  Matrix Length
  n_clusters = 100          ,  #  number of clusters sampled ]
  n_periods  = 4            ,  #  number of periods (1 = no spacetime)
  mean.exposure.months = 100,  #  mean exposure months per cluster
  extent = c(0,1,0,1)       ,  #  xmin,xmax,ymin,ymax
  ncovariates = 3           ,  #  how many covariates to include?
  seed   = NULL             ,
  returnall=TRUE            ) {


  require(fields)
  require(RandomFields)
  #require(ggplot2); theme_set(theme_bw())
  #require(geoR)
  require(data.table)
  require(raster)


  if(length(betas)!=ncovariates+1)
    stop('length(betas)!=ncovariates+1')

  # seed
  if(!is.null(seed)) set.seed(seed)


  # Generate a covariate field, something simple that changes linearly with location
  # cov. raster has a deeper extent..
  cov.raster=stack(makeRandomCovariate())
  if(ncovariates>1)
    for(i in 2:ncovariates)
      cov.raster <- addLayer(cov.raster,makeRandomCovariate())
  names(cov.raster) <- paste0('X',1:ncovariates)

  #plot(cov.raster,main="Theoretical Covariate")
  ## TODO: Incorporate more realistic covariates

  # make a l by l matrix (raster to sample from), exposures, and covariate values
  template<- raster(outer(seq(extent[1],extent[2], l = l),
                          seq(extent[3],extent[4], l = l)))
  samplespace = data.table('x'=rep(coordinates(template)[,1],n_periods),
                           'y'=rep(coordinates(template)[,2],n_periods),
                           't'=rep(seq(1:n_periods),each=l*l))

  # sample frame
  d=data.table( "x"=runif(n_clusters*n_periods, min=extent[1],max=extent[2]),
                "y"=runif(n_clusters*n_periods, min=extent[3],max=extent[4]),
                "exposures"=round(abs(rnorm(n=n_clusters*n_periods,mean=mean.exposure.months,sd=mean.exposure.months/5))), # exposure months
                "period"=rep(1:n_periods,each=n_clusters),
                'int'=1)

  ## TODO: use pop data to get clumpier (and more realistic) sampling locations
  datanames='int'
  for(cov in 1:ncovariates){
    d[,paste0('X',cov)]=raster::extract(cov.raster[[cov]],cbind(d$x,d$y)) # extract covariate value at sampling locations
    datanames=c(datanames,paste0('X',cov))
  }

  if(!is.null(seed)) set.seed(seed)
  # Define the matern object describing the underlying spatial model
  # can have varying underlying SD
  if(length(Sigma2)==1){
    RMmodel = RMmatern(nu    = nu,
                       var   = Sigma2,
                       scale = scale)
    # A different way of doing this (similar to kronecker?)
    # JT calls this the 'innovations' approach, which should be identical to making kronecker covariance matrix
    # Simulate S
    Epsilon1 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      Epsilon1[,t] = RFsimulate(model=RMmodel,
                                x=samplespace$x[samplespace$t==t],
                                y=samplespace$y[samplespace$t==t])@data[,1]
    }

    # Rho
    Epsilon2 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      if(t==1) Epsilon2[,t] = Epsilon1[,t]
      if(t>=2) Epsilon2[,t] = rho * Epsilon1[,t-1] + Epsilon1[,t]
    }

  } else {
    stopifnot(n_periods==length(Sigma2))

    RMmodel = list()
    for(i in 1:n_periods){
      RMmodel[[i]]= RMmatern(nu    = nu,
                             var   = Sigma2[[i]],
                             scale = scale)
    }
    # A different way of doing this (similar to kronecker?)
    # JT calls this the 'innovations' approach, which should be identical to making kronecker covariance matrix
    # Simulate S
    Epsilon1 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      Epsilon1[,t] = RFsimulate(model=RMmodel[[t]],
                                x=samplespace$x[samplespace$t==t],
                                y=samplespace$y[samplespace$t==t])@data[,1]
    }

    # Rho
    Epsilon2 = array(NA, dim=c(l**2,n_periods))
    for(t in 1:n_periods){
      if(t==1) Epsilon2[,t] = Epsilon1[,t]
      if(t>=2) Epsilon2[,t] = rho * Epsilon1[,t-1] + Epsilon1[,t]
    }

  }



  # as vector
  S=as.vector(Epsilon2)


  # Make S into an n_periods rasterStack to extract samples layer

  samplespace$S=S
  r.samplespace=  rasterFromXYZT(samplespace,'S','t')
  names(r.samplespace)=paste0('period',1:n_periods)

  # plot(r.samplespace)
  # sample 'clusters' or sample sites from the sample space
  d$S=NA
  for(t in 1:n_periods)
    d$S[d$period==t] =
    raster::extract(r.samplespace[[t]],cbind(d$x[d$period==t],d$y[d$period==t]))


  ## Calculate linear predictor
  d$linpred = as.matrix(d[,datanames,with=FALSE]) %*% betas
  d$p = plogis(d$linpred+d$S) # logit link

  # Simulate deaths
  d$deaths <- rbinom(n_clusters*n_periods,size=d$exposures, prob=d$p)
  d$mr = d$death/d$exposure

  # save also a raster of true P on the surface for model comparison
  true.mr=samplespace

  samplespace$int=1
  samplespace=cbind(samplespace,raster::extract(cov.raster,cbind(samplespace$x,samplespace$y)))
  samplespace$P  =  plogis(as.matrix(samplespace[,datanames,with=FALSE]) %*% betas+samplespace$S)


  r.true.mr=  stack(rasterFromXYZ(as.matrix(samplespace[t==1,c('x','y','P'),with=F])))
  if(n_periods>1)
    for(r in 2:n_periods)
      r.true.mr=addLayer(r.true.mr,
                         rasterFromXYZ(as.matrix(samplespace[t==r,c('x','y','P'),with=F])))


  names(r.true.mr)<-paste0('period',1:n_periods)

  if(returnall==T){
    return(list(d=d,
                r.true.mr=r.true.mr,
                RMmodel=RMmodel,
                cov.raster=cov.raster,
                S.raster=r.samplespace,
                fullsamplespace=samplespace,
                template=template,
                fe.names=datanames[-1],
                l.s.ratio=median(abs(as.vector(as.matrix(samplespace[,datanames,with=FALSE]) %*% betas))/abs(samplespace$S))))
  } else { return(d) }


}


















############################
############################
## Fit in INLA


fitinla <- function(simobj=mortsim() # takes a full list object from mortsim
){
  # req inla to run, obviously
  require(INLA)
  require(data.table)

  # parse out some useful objects
  dt=simobj$d
  nperiod = unique(dt$period)

  # Meshes
  message('Building Meshes')
  mesh_s = inla.mesh.2d(
    loc=cbind(dt$x,dt$y),
    max.edge=c(0.2,0.2),
    cutoff=0.05)



  nodes=mesh_s$n
  # A is the nxm 'projector matrix' to project the process at the mesh nodes to locations
  message('Building Projector Matrix (A)')

  pcoords = cbind(x=dt$x, y=dt$y)

  coords_periods <- do.call(rbind,
                            replicate(length(nperiod),
                                      pcoords,
                                      simplify = FALSE))

  groups_periods <- rep(nperiod,
                        each = nrow(pcoords))


  # use inla helper functions to project the spatial effect.
  # want the same mesh for all Time periods to match TMB way.
  A <- inla.spde.make.A(
    mesh = mesh_s,
    loc = as.matrix(dt[, c('x', 'y'),with=F]),
    group = dt$period
  )


  # construct an SPDE model with a Matern kernel
  message('Constructing an SPDE model with a Matern kernel')
  spde <- inla.spde2.matern(mesh = mesh_s,
                            alpha = 2)
  space = inla.spde.make.index("space",
                               n.spde = spde$n.spde,
                               n.group = length(nperiod))

  # make a design matrix
  design_matrix <- data.frame(int = 1,dt[,simobj$fe.names,with=F])


  ####
  # create a data stack as a prep for inla to estimate, includes all info so far that inla needs
  message('Making Stack Object')
  stack.obs=inla.stack(
    tag='est',
    data=list(died=dt$deaths), # response
    A=list(A,1), # proj matrix, not sure what the 1 is for
    effects=list(
      space,
      design_matrix)
  )


  # make the formula
  # f_null <- died ~ -1 + int
  # f_lin  <-  reformulate(paste(simobj$fe.names,collapse = ' + '))
  # f_space <- ~  f(space,
  #                 model = spde,
  #                 group = space.group,
  #                 control.group = list(model = 'ar1'))
  # formula <- (f_null + f_lin + f_space)
  Sys.sleep(runif(1,0,4)) # to prevent some cluster issues
  formula <-
    formula(paste0('died ~ -1+int+',
                   (paste(simobj$fe.names,collapse = ' + ')),
                   ' + f(space,
                   model = spde,
                   group = space.group,
                   control.group = list(model = \'ar1\'))'
                   ))
  ####
  # Fit Model
  message('Fitting Model')
  start_time = Sys.time()
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           link = 1,
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec.intercept = 1,
                                       mean.intercept = 0),
                  control.compute = list(dic = TRUE,
                                         cpo = TRUE,
                                         config = TRUE),
                  control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                  family = 'binomial',
                  num.threads = 1,
                  Ntrials = dt$exposures,
                  verbose = TRUE,
                  keep = TRUE)
  model.runtime = (Sys.time() - start_time)

  # Predict back the surface.
  message('Predicting full surface')

  # sample from posterior over latents
  n_draws=100
  draws <- inla.posterior.sample(n_draws, res_fit)

  # get parameter names
  par_names <- rownames(draws[[1]]$latent)

  # index to spatial field and linear coefficient samples
  s_idx <- grep('^s.*', par_names)
  l_idx <- match(sprintf('%s.1', res_fit$names.fixed),
                 par_names)


  # get samples as matrices
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  rownames(pred_l) <- res_fit$names.fixed

  # replicate coordinates and years
  sspc=simobj$fullsamplespace
  coords <- cbind(x=sspc$x,y=sspc$y)

  groups_periods <- sspc$t

  # Projector matrix
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = coords,
    group = groups_periods
  )


  # get samples of s for all coo locations
  s <- A.pred %*% pred_s
  s <- as.matrix(s)


  # predict out linear effects
  vals=as.matrix(cbind(int=1,sspc[,simobj$fe.name,with=F]))
  l <- vals %*% pred_l


  # n_draws of estimates at each location
  y_hat = plogis(l+s)
  sspc$mean_pred_mr = rowMeans(plogis(l+s))

  # make a latent field raster
  res2 = data.table(coords,
                    epsilon=rowMeans(s),
                    t=rep(1:length(nperiod),each=nrow(sspc)/length(nperiod)))
  names(res2)=c('x','y','epsilon','t')
  epsilon_ras=rasterFromXYZT(res2,"epsilon","t")



  return(list(
    inla.object = res_fit,
    fullsamplespace = sspc,
    draws = y_hat,
    nodes=nodes,
    model.runtime=model.runtime,
    mean_ras=rasterFromXYZT(sspc,"mean_pred_mr","t"),
    s.raster=epsilon_ras,
    nodes=mesh_s$n
  ))

}








############################
############################
## Fit in INLA


fitinla2 <- function(dt=d,template=newr # takes a full list object from mortsim
){
  # req inla to run, obviously
  require(INLA)
  require(data.table)

  # parse out some useful objects

  # Meshes
  message('Building Meshes')
  mesh_s = inla.mesh.2d(
    loc=cbind(dt$x,dt$y),
    max.edge=c(0.2,0.2),
    cutoff=0.05)



  nodes=mesh_s$n
  # A is the nxm 'projector matrix' to project the process at the mesh nodes to locations
  message('Building Projector Matrix (A)')

  pcoords = cbind(x=dt$x, y=dt$y)



  # use inla helper functions to project the spatial effect.
  # want the same mesh for all Time periods to match TMB way.
  A <- inla.spde.make.A(
    mesh = mesh_s,
    loc = as.matrix(dt[, c('x', 'y'),with=F])
  )


  # construct an SPDE model with a Matern kernel
  message('Constructing an SPDE model with a Matern kernel')
  spde <- inla.spde2.matern(mesh = mesh_s,
                            alpha = 2)
  space = inla.spde.make.index("space",
                               n.spde = spde$n.spde)

  # make a design matrix
  design_matrix <- data.frame(int = rep(1,nrow(dt)))


  ####
  # create a data stack as a prep for inla to estimate, includes all info so far that inla needs
  message('Making Stack Object')
  stack.obs=inla.stack(
    tag='est',
    data=list(died=dt$deaths), # response
    A=list(A,1), # proj matrix, not sure what the 1 is for
    effects=list(
      space,
      design_matrix)
  )


  # make the formula
  # f_null <- died ~ -1 + int
  # f_lin  <-  reformulate(paste(simobj$fe.names,collapse = ' + '))
  # f_space <- ~  f(space,
  #                 model = spde,
  #                 group = space.group,
  #                 control.group = list(model = 'ar1'))
  # formula <- (f_null + f_lin + f_space)
  Sys.sleep(runif(1,0,4)) # to prevent some cluster issues
  formula <-
    formula(paste0('died ~ -1+int+ f(space,model = spde)'
                   ))
  ####
  # Fit Model
  message('Fitting Model')
  start_time = Sys.time()
  res_fit <- inla(formula,
                  data = inla.stack.data(stack.obs),
                  control.predictor = list(A = inla.stack.A(stack.obs),
                                           link = 1,
                                           compute = FALSE),
                  control.fixed = list(expand.factor.strategy = 'inla',
                                       prec.intercept = 1,
                                       mean.intercept = 0),
                  control.compute = list(dic = TRUE,
                                         cpo = TRUE,
                                         config = TRUE),
                  control.inla = list(int.strategy = 'eb', h = 1e-3, tolerance = 1e-6),
                  family = 'binomial',
                  num.threads = 1,
                  Ntrials = dt$exposures,
                  verbose = TRUE,
                  keep = TRUE)
  model.runtime = (Sys.time() - start_time)

  # Predict back the surface.
  message('Predicting full surface')

  # sample from posterior over latents
  n_draws=100
  draws <- inla.posterior.sample(n_draws, res_fit)

  # get parameter names
  par_names <- rownames(draws[[1]]$latent)

  # index to spatial field and linear coefficient samples
  s_idx <- grep('^s.*', par_names)
  l_idx <- match(sprintf('%s.1', res_fit$names.fixed),
                 par_names)


  # get samples as matrices
  pred_s <- sapply(draws, function (x) x$latent[s_idx])
  pred_l <- sapply(draws, function (x) x$latent[l_idx])
  if(is.null(dim(pred_l))) pred_l=t(pred_l)
  rownames(pred_l) <- res_fit$names.fixed

  # replicate coordinates and years


  coords <- data.table(x= coordinates(template)[,1],y= coordinates(template)[,2])


  # Projector matrix
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = as.matrix(coords)
  )


  # get samples of s for all coo locations
  s <- A.pred %*% pred_s
  s <- as.matrix(s)


  # predict out linear effects
  vals=as.matrix(cbind(int=rep(1,nrow(coords))))
  l <- vals %*% pred_l


  # n_draws of estimates at each location
  y_hat = plogis(l+s)
  coords$mean_pred_mr = rowMeans(plogis(l+s))

  # make a latent field raster
  # res2 = data.table(coords,
  #                   epsilon=rowMeans(s),
  #                   t=rep(1:length(nperiod),each=nrow(sspc)/length(nperiod)))
  # names(res2)=c('x','y','epsilon','t')
  # epsilon_ras=rasterFromXYZT(res2,"epsilon","t")
  #
  #

  mean_ras=rasterFromXYZ(cbind(coords$x,coords$y,coords$mean_pred_mr))

  return(mean_ras)

}









###########
## FIT IN TMB


fittmb <- function(simobj=mortsim(),
                   silent=FALSE,
                   wd='code/',
                   Version = "spde"){

  require(INLA)
  require(TMB)
  require(RandomFields)
  # require(TMBdebug)
  require(data.table)

  message('setup')
  dt = simobj[["d"]]

  dt$id=1:nrow(dt)
  coords = cbind(dt$x,dt$y)
  nperiod = length(unique(dt$period))

  # use same mesh per time point
  mesh_s = inla.mesh.2d(
    loc=coords,
    max.edge=c(0.2,0.2),
    cutoff=0.05)

  nodes=mesh_s$n
  # Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
  spde = inla.spde2.matern( mesh_s )



  # pull covariate(s) at knots
  covs <- raster::extract(simobj$cov.raster,cbind(mesh_s$loc[,1],mesh_s$loc[,2]))



  ###################
  # Parameter estimation
  ###################

  #####  Version 0 -- Sweep upstream to downstream through time


  # make sure I am in the code directory
  if(length(grep('code',getwd()))==0)
    setwd(wd)
  #
  # # need to recompile for different systems, drop compiled files from previoulsy run os if needed
  # ifelse(Sys.info()['sysname']=='Windows',
  #        try(file.remove(c(paste0(Version,c('.so')))),TRUE),   # if running on Windows, drop .so and .o files
  #        try(file.remove(c(paste0(Version,c('.dll')))),TRUE))  # if running on unix, drop .dll and .o files
  #
  # if((file.exists(c(paste0(Version,c('.dll')))))+(file.exists(c(paste0(Version,c('.so')))))!=1) # may need to also delete the .o file
  #   try(file.remove(c(paste0(Version,c('.o')))),TRUE)

  # Compile
  message('compiling likelihood function')
  TMB::compile( paste0(Version,".cpp") )
  dyn.load( dynlib(Version) )

  # Data to pass to TMB
  X_xp = cbind( 1, covs)

  Data = list( n_i=nrow(dt),                   # Total number of observations
               n_x=mesh_s$n,                   # Number of vertices in SPDE mesh
               n_t=nperiod,                    # Number of periods
               n_p=ncol(X_xp),                 # Number of columns in covariate matrix X
               x_s=mesh_s$idx$loc-1,           # Association of each cluster with a given vertex in SPDE mesh
               c_i=dt$deaths,                  # Number of observed deaths in the cluster (N+ in binomial likelihood)
               Exp_i=dt$exposures,             # Number of observed exposures in the cluster (N in binomial likelihood)
               s_i=dt$id-1,                    # no site specific effect in my model (different sampling locations in time)
               t_i=dt$period-1,                # Sample period ( starting at zero )
               X_xp=X_xp,                      # Covariate design matrix
               G0=spde$param.inla$M0,          # SPDE sparse matrix
               G1=spde$param.inla$M1,          # SPDE sparse matrix
               G2=spde$param.inla$M2)          # SPDE sparse matrix



  # staring values for parameters
  Parameters = list(alpha   =  rep(0,ncol(X_xp)),        # FE parameter alpha
                    log_tau_E=1.0,                # log inverse of tau  (Epsilon)
                    #                  log_tau_O=1.0,               # log inverse of tau (SP)
                    log_kappa=0.0,	              # Matern Range parameter
                    rho=0.5,                      # Autocorrelation term
                    epsilon=matrix(1,ncol=nperiod,nrow=mesh_s$n),  # GP
                    sp=matrix(rnorm(mesh_s$n)))         # RE for mesh points



  # which parameters are random
  Random = c("epsilon",'sp')



  # Make object
  message('making object')
  obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)
  if(silent) obj$env$beSilent()

  # Run optimizer
  message('running optimizer')
  start_time = Sys.time()
  opt0 = nlminb(start       =    obj$par,
                objective   =    obj$fn,
                gradient    =    obj$gr,
                lower       =    c(rep(-20,sum(names(obj$par)=='alpha')),rep(-10,2),-0.999),
                upper       =    c(rep(20 ,sum(names(obj$par)=='alpha')),rep(10,2),0.999),
                control     =    list(eval.max=1e4, iter.max=1e4, trace=0))
  model.runtime = (Sys.time() - start_time)
  opt0[["final_gradient"]] = obj$gr( opt0$par )

  # Get standard errors
  message('getting standard errors')
  Report0 = obj$report()
  SD0 = try( sdreport(obj) )


  #### Prediction
  message('making predictions')
  # get back a surface (use epsilon, then interpolate somehow? )
  epsilon=SD0$par.random[names(SD0$par.random)=="epsilon"] # should be mesh_s$n*nperiods long (make sure its in the right order)
  # later  figure out how to get draws of this, for now, must be mean.. see SD0$diag.cov.random (mesh_s$n*5 long.. )

  # get surface to project on to
  pcoords = cbind(x=simobj$fullsamplespace$x, y=simobj$fullsamplespace$y)
  groups_periods <- simobj$fullsamplespace$t



  # use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)



  # values of S at each cell (long by nperiods)
  cell_s <- as.matrix(A.pred %*% epsilon)


  # fixed effects
  npars <- sum(names(opt0$par)=='alpha')

  # extract cell values  from covariates
  vals <- extract(simobj$cov.raster, pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
  vals <- (cbind(int = 1, vals))

  cell_l <- vals %*% opt0$par[1:npars]
  cell_l = rep(cell_l,nperiod) # since there are no time varying components

  # add together linear and st components
  pred <- cell_l + cell_s
  pred <- plogis(as.vector(pred))

  # make them into time bins
  len = length(pred)/nperiod
  res = data.table(pcoords,
                   pred,
                   t=rep(1:nperiod,each=len))
  mean_ras=rasterFromXYZT(res,"pred","t")

  # make a latent field raster
  res2 = data.table(pcoords,
                    epsilon=cell_s,
                    t=rep(1:nperiod,each=len))
  names(res2)=c('x','y','epsilon','t')
  epsilon_ras=rasterFromXYZT(res2,"epsilon","t")


  return(list(
    mean_ras=mean_ras,
    obj=obj,
    opt=opt0,
    pred=pred,
    res=res,
    SD=SD0,
    nodes=nodes,
    model.runtime=model.runtime,
    s.raster=epsilon_ras,
    nodes=mesh_s$n


  ))


}




###########
## FIT IN TMB for Ineq testing


fittmb2 <- function(dt=d,
                   silent=FALSE,
                   wd='code/',
                   Version = "spde",rastertemplate=r2[[1]]){

  require(INLA)
  require(TMB)
  require(RandomFields)
  # require(TMBdebug)
  require(data.table)

  dt$exposures=1000
  dt$period=as.numeric(dt$per)
  dt$deaths=dt$died
  dt$id=1:nrow(dt)

  coords = cbind(dt$x,dt$y)
  nperiod = length(unique(dt$period))

  # use same mesh per time point
  mesh_s = inla.mesh.2d(
    loc=coords,
    max.edge=c(0.2,0.2),
    cutoff=0.05)

  nodes=mesh_s$n
  # Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
  spde = inla.spde2.matern( mesh_s )

  ###################
  # Parameter estimation
  ###################

  #####  Version 0 -- Sweep upstream to downstream through time


  # make sure I am in the code directory
  if(length(grep('code',getwd()))==0)
    setwd(wd)
  #
  # # need to recompile for different systems, drop compiled files from previoulsy run os if needed
  # ifelse(Sys.info()['sysname']=='Windows',
  #        try(file.remove(c(paste0(Version,c('.so')))),TRUE),   # if running on Windows, drop .so and .o files
  #        try(file.remove(c(paste0(Version,c('.dll')))),TRUE))  # if running on unix, drop .dll and .o files
  #
  # if((file.exists(c(paste0(Version,c('.dll')))))+(file.exists(c(paste0(Version,c('.so')))))!=1) # may need to also delete the .o file
  #   try(file.remove(c(paste0(Version,c('.o')))),TRUE)

  # Compile
  Version='spde'
  message('compiling likelihood function')
  TMB::compile( paste0(Version,".cpp") )
  dyn.load( dynlib(Version) )

  # Data to pass to TMB
  X_xp = cbind( int=rep(1,length(mesh_s$loc[,1])))


  Data = list( n_i=nrow(dt),                   # Total number of observations
               n_x=mesh_s$n,                   # Number of vertices in SPDE mesh
               n_t=nperiod,                    # Number of periods
               n_p=ncol(X_xp),                 # Number of columns in covariate matrix X
               x_s=mesh_s$idx$loc-1,           # Association of each cluster with a given vertex in SPDE mesh
               c_i=dt$died,                  # Number of observed deaths in the cluster (N+ in binomial likelihood)
               Exp_i=dt$exposures,             # Number of observed exposures in the cluster (N in binomial likelihood)
               s_i=dt$id-1,                    # no site specific effect in my model (different sampling locations in time)
               t_i=dt$period-1,                # Sample period ( starting at zero )
               X_xp=X_xp,                      # Covariate design matrix
               G0=spde$param.inla$M0,          # SPDE sparse matrix
               G1=spde$param.inla$M1,          # SPDE sparse matrix
               G2=spde$param.inla$M2)          # SPDE sparse matrix



  # staring values for parameters
  Parameters = list(alpha   =  rep(0,ncol(X_xp)),        # FE parameter alpha
                    log_tau_E=1.0,                # log inverse of tau  (Epsilon)
                    #                  log_tau_O=1.0,               # log inverse of tau (SP)
                    log_kappa=0.0,	              # Matern Range parameter
                    rho=0.5,                      # Autocorrelation term
                    epsilon=matrix(1,ncol=nperiod,nrow=mesh_s$n),  # GP
                    sp=matrix(rnorm(mesh_s$n)))         # RE for mesh points



  # which parameters are random
  Random = c("epsilon",'sp')



  # Make object
  message('making object')
  obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, hessian=FALSE, DLL=Version)
  if(silent) obj$env$beSilent()

  # Run optimizer
  message('running optimizer')
  start_time = Sys.time()
  opt0 = nlminb(start       =    obj$par,
                objective   =    obj$fn,
                gradient    =    obj$gr,
                lower       =    c(rep(-20,sum(names(obj$par)=='alpha')),rep(-10,2),-0.999),
                upper       =    c(rep(20 ,sum(names(obj$par)=='alpha')),rep(10,2),0.999),
                control     =    list(eval.max=1e4, iter.max=1e4, trace=0))
  model.runtime = (Sys.time() - start_time)
  opt0[["final_gradient"]] = obj$gr( opt0$par )

  # Get standard errors
  message('getting standard errors')
  Report0 = obj$report()
  SD0 = try( sdreport(obj) )


  #### Prediction
  message('making predictions')
  # get back a surface (use epsilon, then interpolate somehow? )
  epsilon=SD0$par.random[names(SD0$par.random)=="epsilon"] # should be mesh_s$n*nperiods long (make sure its in the right order)
  # later  figure out how to get draws of this, for now, must be mean.. see SD0$diag.cov.random (mesh_s$n*5 long.. )
  l=51
  samplespace = data.table('x'=rep(coordinates(rastertemplate)[,1],nperiod),
                           'y'=rep(coordinates(rastertemplate)[,2],nperiod),
                           't'=rep(seq(1:nperiod),each=l*l))

  # get surface to project on to
              pcoords <- cbind(samplespace$x,samplespace$y)


  groups_periods <- samplespace$t



  # use inla helper functions to project the spatial effect.
  A.pred <- inla.spde.make.A(
    mesh = mesh_s,
    loc = pcoords,
    group = groups_periods)



  # values of S at each cell (long by nperiods)
  cell_s <- as.matrix(A.pred %*% epsilon)


  # fixed effects
  npars <- sum(names(opt0$par)=='alpha')

  # extract cell values  from covariates
  #vals <- extract(simobj$cov.raster, pcoords[1:(nrow(simobj$fullsamplespace)/nperiod),])
  vals <- (cbind(int = rep(1, l*l)))

  cell_l <- vals %*% opt0$par[1:npars]
  cell_l = rep(cell_l,nperiod) # since there are no time varying components

  # add together linear and st components
  pred <- cell_l + cell_s
  pred <- plogis(as.vector(pred))

  # make them into time bins
  len = length(pred)/nperiod
  res = data.table(pcoords,
                   pred,
                   t=rep(1:nperiod,each=len))
  names(res)<-c('x','y','pred','t')
  mean_ras=rasterFromXYZT(res,"pred","t")

  # make a latent field raster
  res2 = data.table(pcoords,
                    epsilon=cell_s,
                    t=rep(1:nperiod,each=len))
  names(res2)=c('x','y','epsilon','t')
  epsilon_ras=rasterFromXYZT(res2,"epsilon","t")


  return(list(
    mean_ras=mean_ras,
    obj=obj,
    opt=opt0,
    pred=pred,
    res=res,
    SD=SD0,
    nodes=nodes,
    model.runtime=model.runtime,
    s.raster=epsilon_ras,
    nodes=mesh_s$n


  ))


}



###
# Wrapper function for full sim and fit

simNfit <- function(  nc=2,
                      betas=c(-1,-.1), # the second will be repeated nc times
                      log_kappa=log(1),
                      SigmaE=.5,
                      rho=.90,
                      exposures = 100,
                      numclusters = 100,
                      numperiods = 4,
                      returnmodelobjects=FALSE ) {

  # rep betas nc times
  betas = c(betas[1],rep(betas[2],nc))

  # simulation
  message('Simulating Data')
  sim=mortsim(
    betas      = betas,
    rho        = rho,
    scale      = exp(log_kappa),
    Sigma2     = SigmaE ** 2,
    n_clusters = numclusters,
    n_periods  = numperiods,
    mean.exposure.months = exposures,
    ncovariates = nc)




  # fit in INLA and TMB
  message('RUNNING INLA')
  mod.inla = fitinla(simobj=sim)
  message('RUNNING TMB')
  mod.tmb  = fittmb (simobj=sim)


  # ll
  message('WRAPPING UP')
  ll=cbind(inla=mod.inla$inla.object$mlik[1,1],
           tmb=-1* mod.tmb$opt$objective)
  row.names(ll)='LL'


  # runtimes
  runtimes=cbind(inla=as.numeric(mod.inla$model.runtime,units='secs'),
                 tmb=as.numeric(mod.tmb $model.runtime,units='secs'))
  row.names(runtimes)= 'runtimes'

  # others
  params=compareparams(tbetas=betas,tSigmaE=SigmaE,tlog_kappa=log_kappa,trho=rho,ncs=nc,
                       inla.obj=mod.inla,
                       tmb.obj=mod.tmb)
  params$diff = data.frame(inla=params$inla.params[,'mean']-params$tru.params,
                           tmb =params$tmb.params[,'mean']-params$tru.params)


  pvm=pv(simobj=sim,
         inla.obj=mod.inla,
         tmb.obj=mod.tmb)

  # objects for models, in case
  inla.object=ifelse(returnmodelobjects,mod.inla$inla.object,0)
  tmb.object =ifelse(returnmodelobjects,mod.tmb$opt,0)
  sim.object =ifelse(returnmodelobjects,sim,0)


  return(list(
    ## things I am curious about over the two models
    params=params,
    pvm=pvm,
    ll=ll,
    runtimes=runtimes,
    # add 95% coverage ( should worsen with worsening Ps)

    ## things I change
    l.s.ratio=sim$l.s.ratio,
    numcovs  = nc,
    numpers  = numperiods,
    numnodes = mod.tmb$nodes,
    meanexpos=exposures,
    mean.death.prob = mean(sim$fullsamplespace$P),

    # objects just in case
    inla.object=inla.object,
    tmb.object =tmb.object,
    sim.object =sim.object)

  )
}








#######
## A function to print Predictive Validity Metrics
pv  = function(simobj=sim,
               inla.obj=mod.inla,
               tmb.obj=mod.tmb){
  require(boot)
  # rmse
  inla.rmse = sqrt(mean((logit(simobj$fullsamplespace$P)-logit(inla.obj$fullsamplespace$mean_pred_mr))**2))
  tmb.rmse = sqrt(mean((logit(simobj$fullsamplespace$P)-logit(tmb.obj$pred))**2))

  # prediction error
  inla.sum.pe = sum(logit(simobj$fullsamplespace$P)-logit(inla.obj$fullsamplespace$mean_pred_mr))
  tmb.sum.pe  = sum(logit(simobj$fullsamplespace$P)-logit(tmb.obj$pred))
  inla.mean.pe = mean(logit(simobj$fullsamplespace$P)-logit(inla.obj$fullsamplespace$mean_pred_mr))
  tmb.mean.pe  = mean(logit(simobj$fullsamplespace$P)-logit(tmb.obj$pred))

  res =   as.matrix(cbind(inla=c(inla.rmse,inla.sum.pe,inla.mean.pe),
                          tmb= c(tmb.rmse,tmb.sum.pe,tmb.mean.pe)))
  row.names(res) = c('rmse','sum of pe','mean of pe')
  return(res)
}



#####
# return a list comparing param estimates
compareparams<-function(inla.obj=mod.inla,
                        tmb.obj=mod.tmb,
                        tbetas=betas,tSigmaE=SigmaE,tlog_kappa=log_kappa,trho=rho,ncs=nc){
  # compare parameters
  # extract parameters from inla
  inla.params = as.matrix(inla.obj$inla.object$summary.fixed[,1:2])
  inla.params = rbind(inla.params,
                      inla.obj$inla.object$internal.summary.hyperpar[,1:2])

  # extract parameters from tmb
  tmb.params = as.matrix(cbind(mean=tmb.obj$SD$par.fixed,
                               sd=sqrt(diag(tmb.obj$SD$cov.fixed))),ncol=2)

  # consistent naming for both
  rownames(inla.params) <- rownames(tmb.params) <-
    c(row.names(inla.obj$inla.object$summary.fixed)[1:(ncs+1)],
      names(tmb.obj$SD$par.fixed)[-(1:(ncs+1))])

  # pull truth
  tru.params = as.matrix(c(tbetas,1/(log(tSigmaE^2)),tlog_kappa,trho))
  row.names(tru.params)=rownames(inla.params)


  return(list(tru.params  = tru.params,
              inla.params = inla.params,
              tmb.params  = tmb.params
  ))
}









#####
## Make population weights
## ~~~~~~~~~~~~~~
makePopulationWeightsRaster <- function(adm,   # raster of adm names or codes
                                        pop){  # population raster

  cell_idx <- cellIdx(brick(adm)[[1]] * 0)
  ad_code  <- extract(adm, cell_idx)
  pop_cell <- extract(pop, cell_idx)


  pop_cell[is.na(pop_cell)] <- 0
  pop_cell[is.na(ad_code)] <- 0

  pop_totals_ad <- tapply(pop_cell, ad_code, sum)

  pop_totals_ad_cell <- as.vector(pop_totals_ad)[match(ad_code, names(pop_totals_ad))]
  pop_totals_ad_cell[pop_totals_ad_cell == 0] <- .1

  disowned_ad <- which(is.na(pop_totals_ad_cell))
  pop_totals_ad_cell[disowned_ad] <- .1

  pop_cell[disowned_ad] <- 0
  pop_wt_ad <- pop_cell / pop_totals_ad_cell
  wt_sum_ad <- tapply(pop_wt_ad, ad_code, sum)
  stopifnot(all.equal(wt_sum_ad, round(wt_sum_ad)))

  return(insertRaster(adm,pop_wt_ad))

}







##################################
##  FUNCTIONS TO GET INEQUALITY METRICS

# gini: exists already as gini() from reldist:
gini=function (x, weights = pop_wt) { # will use apply, so think of doing this for each country
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}


# range ratio (lowlimit and highlimit are quantiles)
# make sure to only give it places with a population
RGR<-function(x,weights,hi=.9,lo=.1) {
  require(reldist)
  unname(reldist::wtd.quantile(x,weight=weights,q=hi,na.rm=T)/
           reldist::wtd.quantile(x,weight=weights,q=lo,na.rm=T))
}



# range (lowlimit and highlimit are quantiles)
# make sure to only give it places with a population
range<-function(x,weights,hi=.9,lo=.1) {
  require(reldist)
  unname(reldist::wtd.quantile(x,weight=weights,q=hi,na.rm=T)-
    reldist::wtd.quantile(x,weight=weights,q=lo,na.rm=T))

}



# range (lowlimit and highlimit are quantiles)
# make sure to only give it places with a population
varr<-function(x,weights) {
  require(Hmisc)
  wtd.var(x,weights=weights,na.rm=T)

}

# coefficient of variation
coefvar<-function(x,weights){
  require(SDMTools)
  wt.sd(x,weights)/wt.mean(x,weights)
}

# coefficient of variation (MAD/MEDIAN)
madmed<-function(x,weights){
  require(matrixStats)
  require(reldist)
  weightedMad(x,w=weights)/
    reldist::wtd.quantile(x,weight=weights,q=.5,na.rm=T)

}


# Generalized code for inter-individual difference metrics, a=1 and b=1 should return gini
# NOT CURRENTLY WORKING WITH WREIGHTS - GIVES DIFF ANWER THAN GINI FUNCTION
IID = function(x,pop=NULL,alpha=1,beta=1){ # weights should be population counts per pixel
  if(is.null(pop)) pop = rep(1,length(x))
  z=x
  xx = outer(z,z,FUN='-')
  numerator = sum(abs(xx)^alpha)
  denom=2*length(z)^2*mean(x*(pop))^beta # need to exlude zero pop from length here?
  return(numerator/denom)
}
# unweighted the above returns the correct gini

# Generalized code for inter-meandifference metrics,
# NOT CURRENTLY WORKING WITH WEIGHTS -
IMD = function(x,alpha=2,beta=0){
  numerator = sum(abs(x-mean(x))^alpha)
  denom=length(x)*(mean(x)^beta) # need to exlude zero pop from length here?
  return(numerator/denom)
}
# unweighted it returns the correct variance
# not correct CV (but maybe a mistake in the paper.. see page 48).


# conditional simulation function
condSim=function (vals, weights = NULL, group = NULL, fun = NULL, hi=NULL,lo=NULL, ...)
{
  ncell <- nrow(vals)
  ndraw <- ncol(vals)
  fun_string <- deparse(substitute(fun))
  if (is.null(weights)) {
    weights <- rep(1, ncell)
  }
  else {
    if (length(weights) != ncell) {
      stop(sprintf("number of elements in weights (%i) not equal to number of cells in vals (%i)",
                   length(weights), ncell))
    }
  }
  if (is.null(group)) {
    group <- rep(1, length(weights))
  }
  else {
    if (length(group) != ncell) {
      stop(sprintf("number of elements in group (%i) not equal to number of cells in vals (%i)",
                   length(group), ncell))
    }
  }
  levels <- unique(na.omit(group))
  nlevel <- length(levels)
  ans <- matrix(NA, ncol = ndraw, nrow = nlevel)
  rownames(ans) <- levels
  for (lvl in 1:nlevel) {
    idx <- which(group == levels[lvl])
    if (is.null(fun)) {
      if (all(dim(t(vals[idx, ])) == c(1, ndraw))) {
        ans[lvl, ] <- weights[idx] %*% t(vals[idx, ])
      }
      else {
        ans[lvl, ] <- weights[idx] %*% vals[idx, ]
      }
    }
    else {
      if(fun_string %in% list('RGR','range')){
        ans[lvl, ] <- apply(as.matrix(vals[idx, ]), 2, fun, weights = weights[idx], hi=hi,lo=lo,
                            ...)
      } else {
        ans[lvl, ] <- apply(as.matrix(vals[idx, ]), 2, fun, weights = weights[idx],
                          ...)
      }
    }
  }
  if (nlevel == 1)
    ans <- as.vector(ans)
  return(ans)
}



## A function to make a random surface
makeRandomSurface <- function(extent = c(0,1,0,1),
                              sd=.1,
                              l=51,
                              scale=1,
                              offset=0){

  require(RandomFields)
  require(raster)

  cov.raster<- raster(outer(seq(0,1, l = l),
                            seq(0,1, l = l),
                            FUN = function(x,y) x*y*0),
                      xmn=extent[1], xmx=extent[2],
                      ymn=extent[3], ymx=extent[4])


  RMmodel = RMmatern(nu    = 1,
                     var   = sd,
                     scale = scale)

  z=RFsimulate(model=RMmodel,
               x=coordinates(cov.raster)[,1],
               y=coordinates(cov.raster)[,2])@data[,1]

  cov.raster=cov.raster+rasterFromXYZ(data.frame(x=coordinates(cov.raster)[,1],
                                                 y=coordinates(cov.raster)[,2],
                                                 z=z))+offset
  return(cov.raster)

}
