library(ggplot2)
library(patchwork)

setwd('~/Maple_Ridge_Inference/Obs_Data/')


## works for each dataset, just change these two lines below
obs_avg_10kb_20kindep <- read.csv('Maple_Ridge_10kbgenesfix_X_A_obs_data_ABCRF.csv') 
file_prefix = 'Maple_Ridge_10kgenesfix'


### Mean Values of Each Summary Statistic #### 
#create dataframes for statistics with X chrom and A value 

pi_df = obs_avg_10kb_20kindep[, which(grepl('pi_m', colnames(obs_avg_10kb_20kindep)))]
tw_df = obs_avg_10kb_20kindep[, which(grepl('theta_m', colnames(obs_avg_10kb_20kindep)))]
r2_df = obs_avg_10kb_20kindep[, which(grepl('r2_m', colnames(obs_avg_10kb_20kindep)))]
r2_var_df = obs_avg_10kb_20kindep[, which(grepl('r2_var_m', colnames(obs_avg_10kb_20kindep)))]
TD_df = obs_avg_10kb_20kindep[, which(grepl('TD_m', colnames(obs_avg_10kb_20kindep)))]

colnames(pi_df) <- c("Autosome", "X_Chr")
colnames(tw_df) <- c("Autosome", "X_Chr")
colnames(r2_df) <- c("Autosome", "X_Chr")
colnames(r2_var_df) <- c("Autosome", "X_Chr")
colnames(TD_df) <- c("Autosome", "X_Chr")

mean_list = list(pi_df, tw_df, r2_df, r2_var_df, TD_df)


stat_list = list(NA)

labs = c('pi', "Watterson's theta", 'r2', 'r2 Variance', "Tajima's D")

for(i in c(1:5)){
  stat_list[[i]] <- data.frame(matrix(nrow = 2, ncol=3))
  colnames(stat_list[[i]]) <- c('Chr', 'Val', 'Label')
  stat_list[[i]]$Chr <- colnames(mean_list[[i]])
  stat_list[[i]]$Val <- unlist(mean_list[[i]][1,], use.names = F)
  stat_list[[i]]$Label <- labs[i]
  
}

plot_list = list(NA)

for (i in c(1:length(stat_list))){
  plot_list[[i]] <- ggplot(stat_list[[i]], aes(x=Chr, y=Val))+
    geom_point()+
    theme_bw()+
    ggtitle(paste("Maple Ridge", stat_list[[i]]$Label[1], sep=' '))+
    ylab(paste("mean",stat_list[[i]]$Label[1]))
}

mean_out_plots <- (plot_list[[1]] | plot_list[[2]] | plot_list[[5]] ) / (plot_list[[3]] | plot_list[[4]])

ggsave(paste(file_prefix,'mean_obs_stats.png', sep="_"), mean_out_plots, width = 9, height = 5, units = 'in')

### Variances of Each Summary Statsitic #####


pi_df = obs_avg_10kb_20kindep[, which(grepl('pi_v', colnames(obs_avg_10kb_20kindep)))]
tw_df = obs_avg_10kb_20kindep[, which(grepl('theta_v', colnames(obs_avg_10kb_20kindep)))]
r2_df = data.frame(matrix(nrow = nrow(tw_df), ncol=ncol(tw_df)))
colnames(r2_df) <- c("Ar2_v", "Xr2_v")
r2_df$Ar2_v <- obs_avg_10kb_20kindep$Ar2_v
r2_df$Xr2_v <- obs_avg_10kb_20kindep$Xr2_v
r2_var_df = obs_avg_10kb_20kindep[, which(grepl('r2_var_v', colnames(obs_avg_10kb_20kindep)))]
TD_df = obs_avg_10kb_20kindep[, which(grepl('TD_v', colnames(obs_avg_10kb_20kindep)))]

colnames(pi_df) <- c("Autosome", "X_Chr")
colnames(tw_df) <- c("Autosome", "X_Chr")
colnames(r2_df) <- c("Autosome", "X_Chr")
colnames(r2_var_df) <- c("Autosome", "X_Chr")
colnames(TD_df) <- c("Autosome", "X_Chr")

var_list = list(pi_df, tw_df, r2_df, r2_var_df, TD_df)


stat_list = list(NA)

labs = c('pi', "Watterson's theta", 'r2', 'r2 Variance', "Tajimas D")

for(i in c(1:5)){
  stat_list[[i]] <- data.frame(matrix(nrow = 2, ncol=3))
  colnames(stat_list[[i]]) <- c('Chr', 'Val', 'Label')
  stat_list[[i]]$Chr <- colnames(var_list[[i]])
  stat_list[[i]]$Val <- unlist(var_list[[i]][1,], use.names = F)
  stat_list[[i]]$Label <- labs[i]
  
}

plot_list = list(NA)

for (i in c(1:length(stat_list))){
  plot_list[[i]] <- ggplot(stat_list[[i]], aes(x=Chr, y=Val))+
    geom_point()+
    theme_bw()+
    ggtitle(paste("Maple Ridge", stat_list[[i]]$Label[1], sep=' '))+
    ylab(paste("variance",stat_list[[i]]$Label[1]))
}

var_out_plots <- (plot_list[[1]] | plot_list[[2]] | plot_list[[5]] ) / (plot_list[[3]] | plot_list[[4]])

ggsave(paste(file_prefix,'var_obs_stats.png', sep="_"), var_out_plots, width = 9, height = 5, units = 'in')

