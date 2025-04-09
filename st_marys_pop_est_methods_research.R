library(sf)
library(ggplot2)
library(RTMB)
library(plyr)
library(fmesher)
library(glmmTMB)

matern <- function (u, phi, kappa) 
{
  if (is.vector(u)) 
    names(u) <- NULL
  if (is.matrix(u)) 
    dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
                                                    gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
                 1)
  if (!RTMB:::ad_context()) ## FIXME
    uphi[u > 600 * phi] <- 0
  return(uphi)
}

##spatial boundaries
regions=st_read('region_polys/Regions_2013.shp')

st_marys_shoreline=st_read("AOC_MI_StMarys_2022/AOC_MI_StMarys_2022.shp")

st_marys_shoreline_f=st_transform(st_marys_shoreline,st_crs(regions))



ggplot() + geom_sf(data=st_marys_shoreline_f,fill='blue') +
  geom_sf(data=regions,colour='orange',fill=NA) + 
  geom_sf_label(data=regions,aes(label=Region))


st_marys_shoreline_region_list=list()
for(i in 1:nrow(regions)){
  river_section=st_intersection(st_marys_shoreline_f,regions[i,])
  river_section$Region=regions$Region[i]
  river_section$area_ha=as.numeric(st_area(river_section)/10000)
  st_marys_shoreline_region_list[[i]]=river_section
}

st_marys_shoreline_region=do.call(rbind,st_marys_shoreline_region_list)

ggplot() + geom_sf(data=st_marys_shoreline_region,fill='blue') +
  geom_sf(data=regions,colour='orange',fill=NA) + 
  geom_sf_label(data=regions,aes(label=Region))

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

##plot catch verses habitat and depth
mean_by_habitat_type=ddply(river_pop_est_data, .(depth,hab.type),summarize, 
                           mean_n=mean(sl.larv.n),n=length(sl.larv.n),
                           sd=sd(sl.larv.n))
mean_by_habitat_type_with_data=mean_by_habitat_type[mean_by_habitat_type$n>100,
                                                    ]
ggplot() + geom_point(data=mean_by_habitat_type_with_data,aes(x=depth,y=mean_n,color=as.factor(hab.type)))

##add spatial region/grid assigment to catch data
st_marys_shoreline_region_ordered=st_marys_shoreline_region[order(st_marys_shoreline_region$Region),]

##grid spatial correlation
spatial_corr_grid=st_as_sf(st_make_grid(st_marys_shoreline_region,cellsize = c(250,250))) #dimension of grid in meters

spatial_corr_grid_in_river=spatial_corr_grid[st_marys_shoreline_region,]

##remove parts of spatial grids that overlap with plot polys
##not needed when a small grid size is used

#spatial_corr_grid_in_river_ll=st_transform(spatial_corr_grid_in_river,4326)
#
#
#plot_polys_old_utm=st_transform(plot_polys_old,st_crs(spatial_corr_grid_in_river))
#spatial_corr_grid_in_river_less_old_polys_list=list()
#
#erase_poly=st_union(st_make_valid(st_combine(plot_polys_old_utm)))
#
#spatial_corr_grid_in_river_less_old_polys=st_difference(spatial_corr_grid_in_river,erase_poly)
#
#spatial_corr_grid_in_river_less_old_polys_f=cbind(spatial_corr_grid_in_river_less_old_polys,as.data.frame(plot_polys_old_utm)[0,1:(length(plot_polys_old_utm)-1)][1:nrow(spatial_corr_grid_in_river_less_old_polys),])
#
#st_geometry(spatial_corr_grid_in_river_less_old_polys_f)<-'geometry'

#spatial_corr_grid_in_river=rbind(spatial_corr_grid_in_river_less_old_polys_f,plot_polys_old_utm)

spatial_corr_grid_in_river$grid_id=1:nrow(spatial_corr_grid_in_river)

ggplot() + geom_sf(data=spatial_corr_grid_in_river) + geom_sf(data=st_marys_shoreline_region)

distmat_km=list(distmat_region=matrix(as.numeric(st_distance(st_marys_shoreline_region_ordered))/1000,
                                      nrow(st_marys_shoreline_region_ordered),
                                      nrow(st_marys_shoreline_region_ordered)),
                distmat_grid=matrix(as.numeric(st_distance(st_centroid(spatial_corr_grid_in_river)))/1000,
                                    nrow(spatial_corr_grid_in_river),
                                    nrow(spatial_corr_grid_in_river)))

distmat_km$distmat_region_scale=distmat_km$distmat_region /max(distmat_km$distmat_region)

distmat_km$distmat_grid_scale=distmat_km$distmat_grid/max(distmat_km$distmat_grid)

##add grid_id to catch data
river_pop_est_data_sf=st_as_sf(river_pop_est_data,coords=c('longitude','latitude'),crs=4326)
river_pop_est_data_utm=st_transform(river_pop_est_data_sf,st_crs(spatial_corr_grid_in_river))

regions_ll=st_transform(regions,4326)
spatial_corr_grid_in_river_ll=st_transform(spatial_corr_grid_in_river,4326)

river_pop_est_data$region_sp_id=regions$Region[as.numeric(st_intersects(river_pop_est_data_utm,regions))]
river_pop_est_data$grid_id=spatial_corr_grid_in_river_ll$grid_id[as.numeric(st_intersects(river_pop_est_data_sf,spatial_corr_grid_in_river_ll))]

##drop catch data outside of region and/or grid polys
river_pop_est_data_sp_filter=river_pop_est_data[is.na(river_pop_est_data$grid_id)==F,]
river_pop_est_data_sp_filter=river_pop_est_data_sp_filter[is.na(river_pop_est_data_sp_filter$region_sp_id)==F,]

river_pop_est_data_sp_filter$year_num=as.numeric(river_pop_est_data_sp_filter$year_f)
river_pop_est_data_sp_filter$region_f=as.factor(river_pop_est_data_sp_filter$region_sp_id)

##id grids greater than max cor distance away from a sample each year
##thought about using this to map sp_dev values to fix at zero to help with convergence/speed
#max_cor_dist=15
#survey_years=as.numeric(levels(as.factor(river_pop_est_data_sp_filter$year_num)))
#effortbygrid=ddply(river_pop_est_data_sp_filter,.(year_num,grid_id), summarize, n_survey=length(grid_id))
#
#sp_dev_2_maxcordist=matrix(0,nrow=max(river_pop_est_data_sp_filter$year_num), ncol=max(spatial_corr_grid_in_river$grid_id))
#for(i in 1:length(survey_years)){
#  grids_surveyed=effortbygrid$grid_id[effortbygrid$year_num==survey_years[i]]
#  grid_min_dist_from_survey=apply(distmat_km$distmat_grid[,grids_surveyed],1,min)
#  grid_further_than_max_cor_dist=which(grid_min_dist_from_survey>max_cor_dist)
#  
#  ##assign a 1 for grids that are > max_cor_dist from a survey location in a given year
#  sp_dev_2_maxcordist[i,grid_further_than_max_cor_dist]=1
#}


##add upper/lower designation
river_pop_est_data_sp_filter$upper=ifelse(river_pop_est_data_sp_filter$region_sp_id %in% c(1,2,5),1,0)
river_pop_est_data_sp_filter$year_num=as.numeric(river_pop_est_data_sp_filter$year_f)
river_pop_est_data_sp_filter$upper_f=as.factor(river_pop_est_data_sp_filter$upper)
river_pop_est_data_sp_filter$hab.type_f=as.factor(river_pop_est_data_sp_filter$hab.type)

##habitat + region model
glm.hab=glm(sl.larv.n ~ 1 + year_f + upper_f + hab.type_f + depth + I(depth^2), 
         offset=log(sample_area_m2) + log(mean_p),family='poisson',data=river_pop_est_data_sp_filter)

glm.hab.nb=glmmTMB(sl.larv.n ~ 1 + year_f + upper_f + hab.type_f + depth + I(depth^2), 
            offset=log(sample_area_m2) + log(mean_p),family='nbinom2',data=river_pop_est_data_sp_filter)

##habitat + region model with grid_id ranef
#poisson doesn't converge
glmm.hab=glmmTMB(sl.larv.n ~ 1 + year_f + upper_f + hab.type_f + depth + I(depth^2) + (1|grid_id), 
                    offset=log(sample_area_m2) + log(mean_p),family='poisson',data=river_pop_est_data_sp_filter)

glmm.hab.nb=glmmTMB(sl.larv.n ~ 1 + year_f + upper_f + hab.type_f + depth + I(depth^2) + (1|grid_id), 
            offset=log(sample_area_m2) + log(mean_p),family='nbinom2',data=river_pop_est_data_sp_filter)

## spde methods
river_points=st_transform(st_centroid(spatial_corr_grid_in_river),4326)
loc <-st_coordinates(river_points)  # Spatial coordinates
bnd1 <- fmesher::fm_nonconvex_hull(loc, convex=0.05)
bnd2 <- fmesher::fm_nonconvex_hull(loc, convex=0.25)
#
mesh <- fmesher::fm_mesh_2d(
  loc=loc,
  boundary=list(bnd1,bnd2),
  min.angle=21,
  max.edge=c(0.05, 0.2),
  cutoff=0.0005,
  plot.delay=0.5
)
mesh$crs<-st_crs(river_points)
#plot(mesh)
spde=fmesher::fm_fem(mesh)

grid2upper=spatial_corr_grid_in_river[st_marys_shoreline_region_ordered[c(1,2,5),],]
spatial_corr_grid_in_river$river_region=ifelse(spatial_corr_grid_in_river$grid_id %in% grid2upper$grid_id,2,1)

ggplot() + geom_sf(data=spatial_corr_grid_in_river,aes(color=as.factor(river_region))) +
   geom_sf(data=st_marys_shoreline_region,fill=NA) +
   geom_sf_label(data=st_marys_shoreline_region,aes(label=Region))


sm_autocorrelated_error_f<-function(pars){
  getAll(pars,rtmb_sm_data)
  "[<-"<-ADoverload("[<-")
  
  nll=0
  
  D_log=matrix(0,nrow=n_years,ncol=n_grid)

 # L=diag(rep(1,n_grid))
  phi=exp(phi_log)
  kappa=exp(kappa_log)
#  cov=matern(distmat,phi,kappa)
 # cov=exp(-distmat/kappa)
#  L[lower.tri(L, diag=F)]=cov[lower.tri(cov, diag=F)]
#  R <- cov2cor(L%*%t(L))  # Correlation matrix of X (guarantied positive definite)
 Q <- phi^2*(kappa^4 * c0 + 2 * kappa^2 * g1 + g2)        ## GMRF prior
  nll <- nll - sum(dgmrf(sp_dev, 0, Q,log=TRUE))       ## Negative log likelihood
 # nll <- nll - sum(dmvnorm(sp_dev,0,R,TRUE)) 
#  nll <- nll - sum(dnorm(sp_dev,0,1,TRUE)) 
  #single variance for whole river
 # sp_dev_tau=sp_dev*exp(sp_sd_log)

 # L2=diag(rep(1,n_region))
#  r_cov=exp(r_cov_log)
  #cov=exp(-distmat/kappa)
#  L2[lower.tri(L2, diag=F)]=r_cov #single covariance parameter
  
#  R2 <- cov2cor(L2%*%t(L2))  # Correlation matrix of X (guarantied positive definite)
  
 # nll <- nll - sum(dnorm(r_dev,0,1,TRUE)) #sum(dmvnorm(r_dev,0,R2,TRUE)) 
  
 # r_dev_tau=r_dev*exp(r_sd_log)

#  D_log[1,]=D0_log[grid2region] + sp_dev_tau[1,] + r_dev_tau[1,grid2region]
  
 # sp_ar1_phi=plogis(sp_ar1_phi_logit)
#  sp_ar1=matrix(0, nrow=n_years,ncol=nrow(c0))
#  sp_ar1[1,]=sp_dev[1,]
#  
#  for(i in 2:n_years){
#   sp_ar1[i,]=sp_ar1_phi*sp_ar1[i-1,] + sp_dev[i,]
#  }
  

  
  paramM=matrix(c(D_year,D_region, habitatTypeBeta, depthBeta, depth_2Beta),ncol=1, 
                nrow=length(c(D_year,D_region, habitatTypeBeta, depthBeta, depth_2Beta)))
  
  predDensity_fixed_eff <-d%*%paramM
  
  expected_catch_log=rep(0,n_data)
  for(i in 1:n_data){
    expected_catch_log[i]=predDensity_fixed_eff[i]+log(survey_area[i])+log(deepwater.ef_p[i]) + sp_dev[year_num[i],meshidxloc[grid_id[i]]] #survey area in m2, change to meshidxloc for spde
  }
  expected_catch=exp(expected_catch_log)
  
  log_var_minus_mu_catch = 2* expected_catch_log - logtheta;
  
  nll <- nll - sum(dnbinom_robust(x=catch, log_mu=expected_catch_log, log_var_minus_mu=log_var_minus_mu_catch,log=TRUE))
  
  #nll <- nll - sum(dpois(catch,expected_catch,TRUE))
  
REPORT(expected_catch)
REPORT(sp_dev)
#REPORT(sp_ar1)
#ADREPORT(sp_ar1)
#ADREPORT(expected_catch)

  return(nll)
}

##spatial correlation random effects
rtmb_sm_data=list(
  region=river_pop_est_data_sp_filter$upper, #upper = 1 lower =0
  distmat=distmat_km$distmat_grid_scale,
  grid_id=river_pop_est_data_sp_filter$grid_id,
  year_num=river_pop_est_data_sp_filter$year_num,
  catch=river_pop_est_data_sp_filter$sl.larv.n,
  deepwater.ef_p=river_pop_est_data_sp_filter$mean_p,
  survey_area=river_pop_est_data_sp_filter$sample_area_m2,
  grid2region=spatial_corr_grid_in_river$river_region,
  n_grid=max(max(spatial_corr_grid_in_river$grid_id)),
  n_years=max(river_pop_est_data_sp_filter$year_num),
  n_region=2,
  n_data=nrow(river_pop_est_data_sp_filter),
  hab.type=river_pop_est_data_sp_filter$hab.type,
  depth=river_pop_est_data_sp_filter$depth,
  c0=spde$c0,
  g1=spde$g1,
  g2=spde$g2,
  meshidxloc=mesh$idx$loc
)

## design matrix
regionID=as.factor(rtmb_sm_data$region)
yearID=as.factor(rtmb_sm_data$year_num)
habitatID=as.factor(rtmb_sm_data$hab.type)
d <- model.matrix( ~ -1 + yearID + regionID + habitatID + rtmb_sm_data$depth + I(rtmb_sm_data$depth^2))
rtmb_sm_data$d=d

pars=list(sp_dev=matrix(0,max(rtmb_sm_data$year_num),nrow(rtmb_sm_data$c0)),
          D_region=rep(0,max(rtmb_sm_data$region)),
          D_year=rep(0,max(rtmb_sm_data$year_num)),
          habitatTypeBeta=0,
          depthBeta=0,
          depth_2Beta=0,
          phi_log=log(0.2),
          kappa_log=log(0.2),
       #   ln_tauO=log(0.2),
          logtheta=log(0.2))
        # sp_ar1_phi_logit=0)
          


sm_autocorrelated_error_f(pars)

#sp_dev_map= sp_dev_2_maxcordist
#sp_dev_map[which(sp_dev_map==0)]=(1:length(which(sp_dev_map==0)))+1
#sp_dev_map[sp_dev_map==1]<-NA
#sp_dev_map=factor(sp_dev_map)
#
#map=list(sp_dev=sp_dev_map)


obj.sp <- RTMB::MakeADFun(sm_autocorrelated_error_f, pars,random=c('sp_dev'))

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
sim_data_mean=obj.sp$report()$expected_catch
sim_data=rnbinom(length(sim_data_mean),mu=sim_data_mean,size=exp(opt.sp$par["logtheta"]))
sim_obs_plot=data.frame(id=1:length(sim_data),
sim=sim_data[order(sim_data)],
obs=river_pop_est_data_sp_filter$sl.larv.n[order(river_pop_est_data_sp_filter$sl.larv.n)])

ggplot() + geom_line(data=sim_obs_plot,aes(y=id,x=sim),color='blue') +
geom_line(data=sim_obs_plot,aes(y=id,x=obs),color='black') +
  scale_y_continuous("cummulative freq") + 
  scale_x_continuous("catch")

##compare parameter estimates -glm, glmm, glmm + spatial cor
  year_coef_plot=data.frame(model=rep(c('sp glmm','glm','glmm'),each=rtmb_sm_data$n_year),
    year=rep(1:rtmb_sm_data$n_year,3),
  year_coef=c(opt.sp$par[which(names(opt.sp$par)=="D_year")],
  glm.hab.nb$fit$par[1:rtmb_sm_data$n_year],
  glmm.hab.nb$fit$par[1:rtmb_sm_data$n_year]))

  ggplot() + geom_line(data=year_coef_plot,aes(x=year,y=year_coef,color=model)) + 
    geom_point(data=year_coef_plot,aes(x=year,y=year_coef,color=model))
  
  
  hab_coef_plot=data.frame(model=rep(c('sp glmm','glm','glmm'),each=4),
                            year=rep(c( "upper river","habitat","depth","depth_2"),3),
                            hab_coef=c(opt.sp$par[c("D_region","habitatTypeBeta","depthBeta","depth_2Beta")],
                                        glm.hab.nb$fit$par[(rtmb_sm_data$n_year+1):(rtmb_sm_data$n_year+4)],
                                        glmm.hab.nb$fit$par[(rtmb_sm_data$n_year+1):(rtmb_sm_data$n_year+4)]))
  
  ggplot() +  geom_jitter(data=hab_coef_plot,aes(x=year,y=hab_coef,color=model),width=0.3)
  
