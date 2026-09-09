library(ggplot2)
library(patchwork)
library(reshape2)


setwd('~/X_Chr_Var_Calling/MR_windows/')


autosomes = c('NC_056008.1', 'NC_056009.1', 'NC_056010.1', 'NC_056011.1', 'NC_056012.1', 'NC_056013.1', 'NC_056014.1', 'NC_056015.1', 'NC_056016.1', 'NC_056017.1', 'NC_056018.1', 'NC_056019.1', 'NC_056020.1', 'NC_056021.1', 'NC_056022.1', 'NC_056023.1', 'NC_056024.1', 'NC_056025.1', 'NC_056026.1', 'NC_056027.1', 'NC_056028.1', 'NC_056029.1', 'NC_056030.1')

x_chr = 'NC_056031.1'

autosome_list = list(NA)

for (i in c(1:length(autosomes))) {
  filename = paste('MR_10kgenes_',autosomes[i],'_div_windows_1000accessible_20k_indep.csv', sep='')
  
  auto_file = read.csv(filename)
  
  autosome_list[[i]] <- auto_file 
  autosome_list[[i]]$Chromosome <- autosomes[i]
}

Autosome_df <- autosome_list[[1]]

for (i in c(2:length(autosomes))) {
  Autosome_df <- rbind(Autosome_df, autosome_list[[i]])
}

Autosome_df <- Autosome_df[sample.int(nrow(Autosome_df), size=2500), ]

Autosome_avg <-  data.frame(matrix(ncol = ncol(Autosome_df)))
colnames(Autosome_avg) <- colnames(Autosome_df)

for (j in c(1:ncol(Autosome_df))) {
  Autosome_avg[1,j] <- mean(Autosome_df[,j], na.rm=T)
}

Autosome_var <- data.frame(matrix(ncol = ncol(Autosome_df)))
colnames(Autosome_var) <- colnames(Autosome_df)

for (j in c(1:ncol(Autosome_df))) {
  Autosome_var[1,j] <- var(Autosome_df[,j], na.rm=T)
}



# Calculate the proportion of variants for each frequecy bin 
total_vars_A <- sum(Autosome_df$var_counts)

Prop_SFS <- Autosome_avg[,which(grepl("Number_vars_bin", colnames(Autosome_avg)))]
A_bins <- c(0:(ncol(Prop_SFS)-1))


for (bin in A_bins){
  col_bin <- paste('Number_vars_bin', bin, sep="_")
  Prop_SFS[1, col_bin] <- sum(Autosome_df[,col_bin], na.rm = T) / total_vars_A
}


X_df <- read.csv('MR_10kgenes_NC_056031.1_div_windows_1000accessible_20k_indep.csv')

X_avg <- data.frame(matrix(ncol=ncol(X_df)))
colnames(X_avg) <- colnames(X_df)

for (j in c(1:ncol(X_df))){
  X_avg[1,j] <- mean(X_df[,j], na.rm =T)
}

X_var <- data.frame(matrix(ncol=ncol(X_df)))
colnames(X_var) <- colnames(X_df)

for (j in c(1:ncol(X_df))){
  X_var[1,j] <- var(X_df[,j], na.rm =T)
}
# Calculate the proportion of variants for each frequecy bin 
total_vars_X <- sum(X_df$var_counts)
X_Prop_SFS <- X_avg[,which(grepl("Number_vars_bin", colnames(X_avg)))]
X_bins <- c(0:(ncol(X_Prop_SFS)-1))
for (bin in X_bins){
  col_bin <- paste('Number_vars_bin', bin, sep="_")
  X_Prop_SFS[1, col_bin] <- sum(X_df[,col_bin], na.rm = T) / total_vars_X
}



## format the dataset to match the ABC simbank data ####
ABC_obs_colnames <- c('name','Api_m','Api_v','Atheta_m','Atheta_v','Ar2_m','Ar2_v','Ar2_var_m','Ar2_var_v','ATD_m','ATD_v','A_Prop_0','A_Prop_1','A_Prop_2','A_Prop_3','A_Prop_4','A_Prop_5','A_Prop_6','A_Prop_7','A_Prop_8','A_Prop_9','A_Prop_10', 'A_Prop_11', 'A_Prop_12', 'A_Prop_13', 'A_Prop_14', 'A_Prop_15', 'A_Prop_16', 'A_Prop_17','Xpi_m','Xpi_v','Xtheta_m','Xtheta_v','Xr2_m','Xr2_v','Xr2_var_m','Xr2_var_v','XTD_m','XTD_v','X_Prop_0','X_Prop_1','X_Prop_2','X_Prop_3','X_Prop_4','X_Prop_5','X_Prop_6','X_Prop_7','X_Prop_8','X_Prop_9','X_Prop_10','X_Prop_11', 'X_Prop_12')


ABC_obs_data <- data.frame(matrix(nrow=1, ncol = length(ABC_obs_colnames)))
colnames(ABC_obs_data) <- ABC_obs_colnames

ABC_obs_data$name <- "Maple Ridge"

# populate the dataframe with means 
ABC_obs_data$Api_m <- Autosome_avg$pi
ABC_obs_data$Atheta_m <- Autosome_avg$Watterson_theta
ABC_obs_data$ATD_m <- Autosome_avg$Tajimas_D
ABC_obs_data$Ar2_m <- Autosome_avg$mean_r2
ABC_obs_data$Ar2_var_m <- Autosome_avg$var_r2

ABC_obs_data$Xpi_m <- X_avg$pi
ABC_obs_data$Xtheta_m <- X_avg$Watterson_theta
ABC_obs_data$XTD_m <- X_avg$Tajimas_D
ABC_obs_data$Xr2_m <- X_avg$mean_r2
ABC_obs_data$Xr2_var_m <- X_avg$var_r2

# populate the dataframe with variances
ABC_obs_data$Api_v <- Autosome_var$pi
ABC_obs_data$Atheta_v <- Autosome_var$Watterson_theta
ABC_obs_data$ATD_v <- Autosome_var$Tajimas_D
ABC_obs_data$Ar2_v <- Autosome_var$mean_r2
ABC_obs_data$Ar2_var_v <- Autosome_var$var_r2

ABC_obs_data$Xpi_v <- X_var$pi
ABC_obs_data$Xtheta_v <- X_var$Watterson_theta
ABC_obs_data$XTD_v <- X_var$Tajimas_D
ABC_obs_data$Xr2_v <- X_var$mean_r2
ABC_obs_data$Xr2_var_v <- X_var$var_r2

# populate Autosome proportion of variants in each SFS bins 
for (bin in A_bins){
  bin_name_P <- paste('Number_vars_bin', bin, sep="_")
  bin_name_ABC <- paste("A_Prop", bin, sep="_")
  ABC_obs_data[,bin_name_ABC] <- Prop_SFS[,bin_name_P]
}

for (bin in X_bins){
  bin_name_P <- paste('Number_vars_bin', bin, sep="_")
  bin_name_ABC <- paste("X_Prop", bin, sep="_")
  ABC_obs_data[,bin_name_ABC] <- Prop_SFS[,bin_name_P]
}

# Write the data to a csv, to be used later
write.csv(ABC_obs_data, file = "../../Maple_Ridge_Inference/ABC_RF/Maple_Ridge_10kbgenes_20kindep_X_A_obs_data_ABCRF.csv", row.names = F)


## Visualize some of the stasitics
# 
# # Change stat to reflect the pattern of variation of interest 
# stat = 'Watterson_theta'
# 
# stat_df <- data.frame(matrix(nrow=1, ncol=2))
# colnames(stat_df) <- c("Autosome", "X_Chromosome")
# 
# stat_df$Autosome <- Autosome_avg[,stat]
# stat_df$X_Chromosome <- X_avg[,stat]
# 
# stat_plot <- melt(stat_df, variable.name = "Chromosome", value.name = stat)
# 
# 
# theta_W <- ggplot(stat_plot, aes(x=Chromosome, y = Watterson_theta))+
#   geom_point(size=4)+
#   theme_bw()+
#   ylab(stat)+
#   ggtitle(paste("Maple Ridge ",stat," X/A Comparison", sep=''))
# 
# 
# (pi + theta_W + TD) / (mean_r2 + var_r2) + plot_annotation(title = "Maple Ridge Mean Stats")
# 


## Create SFS for dadi input #### 

# autosome sfs 
Auto_df_n_vars <- Autosome_df[,which(grepl("Number_vars", colnames(Autosome_df)))]
auto_fsfs <- data.frame(matrix(data=NA,nrow=1, ncol = ncol(Auto_df_n_vars)))
colnames(auto_fsfs) <- colnames(Auto_df_n_vars)

for(j in c(1:ncol(auto_fsfs))){
  auto_fsfs[1,j] <- sum(Auto_df_n_vars[,j], na.rm = TRUE)
}

write.csv(auto_fsfs, file = '../../Maple_Ridge_Inference/Obs_Data/MR_10kbgenes_20kindep_Autosome_FSFS.csv', row.names = F)

# x sfs 
X_df_n_vars <- X_df[,which(grepl('Number_vars', colnames(X_df)))]
x_fsfs <- data.frame(matrix(data=NA,nrow=1,ncol=ncol(X_df_n_vars)))
colnames(x_fsfs) <- colnames(X_df_n_vars)

for(j in c(1:ncol(x_fsfs))){
  x_fsfs[1,j] <- sum(X_df_n_vars[,j], na.rm = TRUE)
}

## same X fsfs for the downsampled A and regular. Just trying to decrease the delta between X and A data
write.csv(x_fsfs, file='../../Maple_Ridge_Inference/Obs_Data/MR_10kbgenes_20kindep_X_FSFS.csv', row.names = F)
#write.csv(x_fsfs, file='../../Maple_Ridge_Inference/Obs_Data/MR_10kbgenes_20kindep_downsampledA_X_FSFS.csv', row.names = F)
# 
# ### After initial exploration, sample fewer Autosomal windows to make the ratio of sequence on the X:A a little larger 
# Auto_df_down <- Autosome_df[sample.int(nrow(Autosome_df), size=1500), ]
# Auto_down_n_var <- Auto_df_down[,which(grepl("Number_vars", colnames(Auto_df_down)))]
# auto_fsfs_down <- data.frame(matrix(data=NA,nrow=1, ncol = ncol(Auto_down_n_var)))
# colnames(auto_fsfs_down) <- colnames(Auto_down_n_var)
# 
# for (j in c(1:ncol(auto_fsfs_down))){
#   auto_fsfs_down[1,j] <- sum(Auto_down_n_var[,j], na.rm=T)
# }
#write.csv(auto_fsfs_down, file='../../Maple_Ridge_Inference/Obs_Data/MR_10kbgenes_20kindep_downsampledA_Autosome_FSFS.csv', row.names = F)

# ratio of X sequence to Auto sequence
X_bases <- sum(X_df$n_bases, na.rm = T)
#Auto_down_bases <- sum(Auto_df_down$n_bases, na.rm = T)
Auto_full_bases <- sum(Autosome_df$n_bases, na.rm = T)



## create csv of list of windows for autosome and X chromosome to export and filter VCF with 
A_windows = Autosome_df[,c('Chromosome', 'window_start', 'window_end')]
X_windows = X_df[,c("window_start", "window_end")]
X_windows$Chromosome = 'NC_056031.1'

write.csv(A_windows, file='Maple_Ridge_Autosome_ind_windows.csv', row.names = F)
write.csv(X_windows, file='Maple_Ridge_X_ind_windows.csv', row.names = F)






