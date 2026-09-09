setwd('/mnt/payseurgenomes/will/Maple_Ridge_Inference/simbank')

library(abcrf)

popname <- 'Maple_Ridge'
filename <- 'Maple_Ridge_a2_Xloci154'

start_time <- Sys.time()

## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

## read in simbank using notation from abcrf documentation
data <- read.csv('noNA_simbank_Maple_Ridge_alpha2_Xloci154_v5_no_phase.csv')
# read in observed data
obs_data <- read.csv("../obs_data/Maple_Ridge_X_A_obs_data_ABCRF.csv")
print("Read in Data")

## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

# ## for troubleshooting ##
# setwd("~/ABC_Dadi_Combined/ABC_Results/")
# data <- read.csv("dadi_terms_simbank_pilot.csv")
# obs_data <- read.csv('pop_red_pf0.5_alpha1_abc_obs_data_v5.csv')


data$name = popname
data <- data[,which(!grepl("Prop_0", colnames(data)))]
simbank_name <- "dadi_terms_Maple_Ridge"


## ABC-RF Models ####
model.BSR <- readRDS(paste("../models/BSR_model",simbank_name,"no_phase.rds", sep='_'))
model.NC_NA <- readRDS(paste("../models/NC_NA_model",simbank_name,"no_phase.rds", sep='_'))
model.NM_NA <- readRDS(paste("../models/NM_NA_model",simbank_name,"no_phase.rds", sep='_'))
model.t_curr <- readRDS(paste("../models/t_curr_model",simbank_name,"no_phase.rds", sep='_'))
model.t_main <-  readRDS(paste("../models/t_main_model",simbank_name,"no_phase.rds", sep='_'))

### Training Data ####
data.BSR <- data.frame(BSR = data[, "beta"], data[, 2:50])
data.NC_NA <- data.frame(NC_NA = data[, 'NC_NA_ratio'], data[,2:50])
data.NM_NA <- data.frame(NM_NA = data[, 'Nmain_NA_ratio'], data[,2:50])
data.t_curr <- data.frame(t_curr = data[,'t_curr_NA_ratio'], data[,2:50])
data.t_main <- data.frame(t_main = data[,'t_main_NA_ratio'], data[,2:50])



## Output the posterior distribution summaries #####

# reorder obs_data to match colnames of simbank data
obs_data <- obs_data[,names(data[, 1:50])]
obs_data$name = popname


print("Beginning Predictions ")


## Iterate over simulated observed data to predict posterior distribution ## 


pred.obsBSR <- predict(object = model.BSR, obs = obs_data, training = data.BSR, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
pred.obst_curr <- predict(object = model.t_curr, obs = obs_data, training = data.t_curr, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
pred.obst_main <- predict(object = model.t_main, obs = obs_data, training = data.t_main, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
pred.obsNC_NA <- predict(object = model.NC_NA, obs = obs_data, training = data.NC_NA, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
pred.obsNM_NA <- predict(object = model.NM_NA, obs = obs_data, training = data.NM_NA, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)


BSR_posterior_pred <- as.data.frame(pred.obsBSR)
t_curr_posterior_pred <- as.data.frame(pred.obst_curr)
t_main_posterior_pred <- as.data.frame(pred.obst_main)
NC_NA_posterior_pred <- as.data.frame(pred.obsNC_NA)
NM_NA_posterior_pred <- as.data.frame(pred.obsNM_NA)


## Save datasets of predictions #### 
write.csv(BSR_posterior_pred, file = paste("../posterior_ests/BSR_posterior_pred", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(t_curr_posterior_pred, file = paste("../posterior_ests/t_curr_posterior_pred_", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(t_main_posterior_pred, file = paste("../posterior_ests/t_main_posterior_pred_", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(NC_NA_posterior_pred, file = paste("../posterior_ests/NC_NA_posterior_pred", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(NM_NA_posterior_pred, file = paste("../posterior_ests/NM_NA_posterior_pred", filename, "no_phase.csv", sep = '_'), row.names = F)


### Compile predicted value and simulated values into one dataframe for output
abc_perform <- data.frame(matrix(nrow = nrow(BSR_posterior_pred), ncol=5))
colnames(abc_perform) <- c("pred_beta", "pred_t_curr", "pred_NC_NA", 'pred_t_main', 'pred_NM_NA')

## popululate output dataframe with expected value of posterior 
abc_perform$pred_beta <- BSR_posterior_pred$median
abc_perform$pred_t_curr <- t_curr_posterior_pred$median
abc_perform$pred_t_main <- t_main_posterior_pred$median
abc_perform$pred_NC_NA <- NC_NA_posterior_pred$median
abc_perform$pred_NM_NA <- NM_NA_posterior_pred$median

write.csv(abc_perform, file=paste("../posterior_ests/abcrf_", filename, "_no_phase_posterior_medians.csv", sep=""), row.names = F)


# Vis approx of posterior distributions for each parameter 

png(filename = paste("../plots/",filename,"_BSR_est_posterior_dist.png",sep=""), width = 1400, height = 1200, res = 200)
densityPlot(object = model.BSR, obs = obs_data, training= data.BSR)
dev.off()

png(filename = paste("../plots/",filename,"_NM_NA_est_posterior_dist.png",sep=""), width = 1400, height = 1200, res = 200)
densityPlot(object = model.NM_NA, obs=obs_data, training=data.NM_NA)
dev.off()

png(filename = paste("../plots/",filename,"_NC_NA_est_posterior_dist.png",sep=""), width = 1400, height = 1200, res = 200)
densityPlot(object = model.NC_NA, obs = obs_data, training = data.NC_NA)
dev.off()

png(filename = paste("../plots/",filename,"_t_curr_est_posterior_dist.png",sep=""), width = 1400, height = 1200, res = 200)
densityPlot(object = model.t_curr, obs = obs_data, training = data.t_curr)
dev.off()

png(filename = paste("../plots/",filename,"_t_main_est_posterior_dist.png",sep=""), width = 1400, height = 1200, res = 200)
densityPlot(object = model.t_main, obs = obs_data, training = data.t_main)
dev.off()


