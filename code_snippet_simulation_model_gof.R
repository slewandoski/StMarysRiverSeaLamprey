library(plyr)
library(ggplot2)
##goodness of fit code snippet
gof_sims=500

sim_data_mean=exp(obj$report()$ln_pred_i)
est_theta=exp(opt$par["ln_theta"])
gof_list=list()
for(g in 1:gof_sims){
  
  sim_data=rnbinom(length(sim_data_mean),mu=sim_data_mean,size=est_theta)
  gof_list[[g]]=data.frame(sim_id=g,
                           id=1:length(sim_data),
                           sim=sim_data[order(sim_data)])
  
  #ggplot() + geom_line(data=sim_obs_plot,aes(y=id,x=sim),color='blue') +
  #  geom_line(data=sim_obs_plot,aes(y=id,x=obs),color='black') +
  #  scale_y_continuous("cummulative freq") + 
  #  scale_x_continuous("catch")
}

gof_df=do.call(rbind,gof_list)
sim_summary=ddply(gof_df, .(id), summarize, mean=mean(sim),lcl=quantile(sim,0.025),ucl=quantile(sim,0.975))

sim_gof_plot_data=data.frame(id=1:length(sim_data),
                             sim_mean=sim_summary$mean,
                             sim_lcl=sim_summary$lcl,
                             sim_ucl=sim_summary$ucl,
                             obs=data$c_i[order(data$c_i)])
sim_gof_plot_data$cum_freq=sim_gof_plot_data$id/length(sim_gof_plot_data$id)

ggplot() + geom_ribbon(data=sim_gof_plot_data,aes(y=cum_freq,xmax=sim_ucl,xmin=sim_lcl),color='orange',alpha=0.3) +
  geom_line(data=sim_gof_plot_data,aes(y=cum_freq,x=sim_mean),color='orange') +
  geom_line(data=sim_gof_plot_data,aes(y=cum_freq,x=obs),color='black') +
  scale_y_continuous("cummulative freq") + 
  scale_x_continuous("catch")

ggplot() + geom_ribbon(data=sim_gof_plot_data,aes(y=cum_freq,xmax=sim_ucl,xmin=sim_lcl),color='orange',alpha=0.3) +
  geom_line(data=sim_gof_plot_data,aes(y=cum_freq,x=sim_mean),color='orange') +
  geom_line(data=sim_gof_plot_data,aes(y=cum_freq,x=obs),color='black') +
  scale_y_continuous("cummulative freq") + 
  scale_x_continuous("catch") + scale_y_continuous(limits=c(0.88,1.0))


##find the point of the curve with maximum distance between cummulative frequencies
obs_sim_diff=sim_gof_plot_data$obs-sim_gof_plot_data$sim_mean
sim_gof_plot_data[which.max(abs(obs_sim_diff)),]
