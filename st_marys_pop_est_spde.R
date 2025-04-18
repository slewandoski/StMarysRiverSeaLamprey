library(sf)
library(ggplot2)
library(RTMB)
library(plyr)
library(fmesher)
library(glmmTMB)

##catch data
##read in catch data
catch=read.csv('all_data_through_2024/2024-10-24-CatchALL.csv')
lengths=read.csv('all_data_through_2024/2024-10-24-LengthALL.csv')

#calc p from Bergstedt and Genovese (2003)
select_intercept=2.2429
select_slope= -0.0164
lengths$p=plogis(select_intercept + lengths$length*select_slope)
lengths$adj_1=1/lengths$p

#calc curently used to adjust catch
var1=0.0229
var2=-1.732
lengths$adj_2=1 + exp(var1*lengths$length + var2)
lengths$p2=1/lengths$adj_2

mean_p_by_sample=ddply(lengths,.(sampid), summarize, mean_p=mean(p))
catch$mean_p=mean_p_by_sample$mean_p[match(catch$sampid,mean_p_by_sample$sampid)]

##how to model detection probability for samples with 0 catch? simplest approach here
catch$mean_p[is.na(catch$mean_p)]=mean(mean_p_by_sample$mean_p)

#type 3 habitat does not get surveyed -remove the effort
catch$sample_area_m2=ifelse(catch$hab.type==3,0,2.44)
catch$year_f=as.factor(catch$year)

#period -1 is pre control sampling. Several years early in the time series were sampled both before and after control
catch_no_pre=catch[catch$period!= -1,]

river_pop_est_data=catch_no_pre[catch_no_pre$sample_area_m2>0,]
river_pop_est_data=river_pop_est_data[river_pop_est_data$region!=6,]
river_pop_est_data=river_pop_est_data[is.na(river_pop_est_data$longitude)==F,]

##plot catch verses habitat and depth
mean_by_habitat_type=ddply(river_pop_est_data, .(depth,hab.type),summarize, 
                           mean_n=mean(sl.larv.n),n=length(sl.larv.n),
                           sd=sd(sl.larv.n))
mean_by_habitat_type_with_data=mean_by_habitat_type[mean_by_habitat_type$n>100,]
ggplot() + geom_point(data=mean_by_habitat_type_with_data,aes(x=depth,y=mean_n,color=as.factor(hab.type)))


##add upper/lower designation
river_pop_est_data$upper=ifelse(river_pop_est_data$region %in% c(1,2,5),1,0)
river_pop_est_data$year_num=as.numeric(river_pop_est_data$year_f)
river_pop_est_data$upper_f=as.factor(river_pop_est_data$upper)
river_pop_est_data$region_f=as.factor(river_pop_est_data$region)
river_pop_est_data$hab.type_f=as.factor(river_pop_est_data$hab.type)

##habitat + region model
glm.hab=glm(sl.larv.n ~ 1 + year_f + region_f + hab.type_f + depth + I(depth^2), 
            offset=log(sample_area_m2) + log(mean_p),family='poisson',data=river_pop_est_data)

glm.hab.nb=glmmTMB(sl.larv.n ~ 1 + year_f + region_f + hab.type_f + depth + I(depth^2), 
                   offset=log(sample_area_m2) + log(mean_p),family='nbinom2',data=river_pop_est_data)

river_pop_est_data_sf=st_as_sf(river_pop_est_data,coords=c('longitude','latitude'),crs=4326)
river_pop_est_data_utm=st_transform(river_pop_est_data_sf,crs = 26916)
#river_pop_est_data_utm=river_pop_est_data_utm[river_pop_est_data_utm$year_num %in% c(1,2,3,4,5,6,7,9,10,11),]
#river_pop_est_data_utm$year_num=as.numeric(as.factor(river_pop_est_data_utm$year_num))
#problem years: 8, 
## spde methods
#--------------------------------------------------------------------------------
# set up spde approximation
#--------------------------------------------------------------------------------
mesh <- fm_mesh_2d(river_pop_est_data_utm, refine = TRUE, cutoff = 500)
mesh$n # mesh nodes --> random effects we must estimate, so coarse to start
plot(mesh, main = "study area with mesh", 
     xlab = "Easting", ylab = "Northing")
points(st_coordinates(river_pop_est_data_utm), cex = 0.1, pch = 1, col = "steelblue4")

# create matrices in fmesher / INLA
spde <- fm_fem(mesh, order = 2)

# create projection matrix from vertices to sample locations
A_is <- fm_evaluator(mesh, loc = river_pop_est_data_utm)$proj$A


#--------------------------------------------------------------------------------
# fit spatial GLMM using SPDE approximation via RTMB
#--------------------------------------------------------------------------------
## design matrix
river_pop_est_data_utm=river_pop_est_data_utm[order(river_pop_est_data_utm$year_num),]
regionID=as.factor(river_pop_est_data_utm$region)
yearID=as.factor(river_pop_est_data_utm$year_num)
habitatID=as.factor(river_pop_est_data_utm$hab.type)
d <- model.matrix( ~ -1 + yearID + regionID + habitatID + river_pop_est_data_utm$depth + I(river_pop_est_data_utm$depth^2))

#ragged array indexing for linking data to random effects
samplesPerYear=as.data.frame(table(river_pop_est_data_utm$year_num))$Freq

data <- list(
  "c_i" = river_pop_est_data_utm$sl.larv.n, 
  "A_is" = A_is, 
  "M0" = spde$c0, "M1" = spde$g1, "M2" = spde$g2,
  "d"=d,
  "survey_area"=river_pop_est_data_utm$sample_area_m2,
  "deepwater.ef_p"=river_pop_est_data_utm$mean_p,
  "A2yr"=cumsum(samplesPerYear),
  "nperyear"=samplesPerYear,
  "n_years"=max(river_pop_est_data_utm$year_num),
  "n_data"=nrow(river_pop_est_data_utm)
)

par <- list(
  "ln_tau" = 0, "ln_kappa" = 0,
  "omega_s" = matrix(0,data$n_years,nrow(spde$c0)),
  "beta_depth"=0,"beta_depth2"=0,"beta_habitat"=0, "beta_region"=rep(0,max(river_pop_est_data_utm$region-1)),
  "beta_year"=rep(0,max(river_pop_est_data_utm$year_num))
)


sm_f<-function(par){
  getAll(data, par, warn = FALSE)
  
  c_i <- OBS(c_i)
  sigE <- 1 / sqrt(4 * pi * exp(2 * ln_tau) * exp(2 * ln_kappa))
  Q <- exp(4 * ln_kappa) * M0 + 2 * exp(2 * ln_kappa) * M1 + M2
  kappa=exp(ln_kappa)
  tau=exp(ln_tau)
 # Q <- tau^2*(kappa^4 * M0 + 2 * kappa^2 * M1 + M2)        ## GMRF prior
  jnll <- 0
  jnll <- jnll - sum(dgmrf(omega_s, 0.0, Q, TRUE, scale = 1 / exp(ln_tau)))
  
  omega_i=matrix(0,1,n_data)
  pos=1
  for(i in 1:n_years){
  omega_i[1,pos:A2yr[i]] <- as.vector(A_is[pos:A2yr[i],] %*% omega_s[i,])
  pos=pos+nperyear[i]
  }
  
  paramM=matrix(c(beta_year,beta_region,beta_habitat, beta_depth, beta_depth2),ncol=1, 
                nrow=length(c(beta_year,beta_region,beta_habitat, beta_depth, beta_depth2)))
  
  predDensity_fixed_eff <-d%*%paramM
  
  expected_catch=exp(predDensity_fixed_eff + as.vector(omega_i) + log(survey_area)+log(deepwater.ef_p))
  
  jnll <- jnll - sum(dpois(c_i, expected_catch,TRUE))
  range <- sqrt(8) / exp(ln_kappa)
  ADREPORT(range)
  REPORT(sigE)
  REPORT(expected_catch)
  jnll
}

sm_f(par)


map=list(#beta_region=factor(c(NA,NA,NA,NA)),
        # beta_year=factor(c(1,rep(NA,28))),
         ln_theta=factor(NA))


obj.sp <- RTMB::MakeADFun(sm_f, par,random=c('omega_s'))
obj.sp$gr()

start_time=Sys.time()
opt.sp <- nlminb(obj.sp$par, obj.sp$fn, obj.sp$gr,control = list(rel.tol=1e-5, eval.max=350,iter.max=350))
#a 0 for convergence indicates successful convergence
opt.sp

run.time=Sys.time()-start_time
run.time

invert_time_start=Sys.time()
sdr=sdreport(obj.sp)
sdr

invert.time=Sys.time()-invert_time_start
invert.time

##goodness of fit
expected_catch=obj.sp$report()$expected_catch
sim_data=rnbinom(length(expected_catch),mu=expected_catch,size=exp(0))
sim_obs_plot=data.frame(id=1:length(sim_data),
                        sim=sim_data[order(sim_data)],
                        obs=river_pop_est_data$sl.larv.n[order(river_pop_est_data$sl.larv.n)])

ggplot() + geom_line(data=sim_obs_plot,aes(y=id,x=sim),color='blue') +
  geom_line(data=sim_obs_plot,aes(y=id,x=obs),color='black') +
  scale_y_continuous("cummulative freq") + 
  scale_x_continuous("catch")


#--------------------------------------------------------------------------------
# approximate overdispersion
# phi = sum( (y_i - mu_i)^2 / mu_i ) / (n - p)
# where:
#   y_i    = observed count
#   mu_i   = expected count (from model)
#   n      = number of observations
#   p      = total number of estimated parameters (fixed + random effects)
# note we could ADREPORT this but I am being lazy tonight
#--------------------------------------------------------------------------------

# compute Pearson residuals
resid_pearson <- (data$c_i - expected_catch) / sqrt(expected_catch)

# overdispersion estimate (phi)
phi <- sum(resid_pearson^2) / (length(data$c_i) - mesh$n - 1)
print(phi) # values >> 1 are bad

##compare parameter estimates -glm, glmm, glmm + spatial cor
year_coef_plot=data.frame(model=rep(c('sp glmm','glm'),each=data$n_year),
                          year=rep(1:data$n_year,2),
                          year_coef=c(opt.sp$par[which(names(opt.sp$par)=="beta_year")],
                                      glm.hab.nb$fit$par[1:data$n_year]))

ggplot() + geom_line(data=year_coef_plot,aes(x=year,y=year_coef,color=model)) + 
  geom_point(data=year_coef_plot,aes(x=year,y=year_coef,color=model))

