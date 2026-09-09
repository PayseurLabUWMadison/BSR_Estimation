setwd('/mnt/ssd/wspurley/ABCRF_Eval/simbanks')

library(abcrf)

popname <- 'sim_bank_const_size_alpha3_Xloci100'
filename <- 'ABC_combined_dadi_terms_const_size_alpha3_Xloci100_fembias'

start_time <- Sys.time()

## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

## read in simbank using notation from abcrf documentation
data <- read.csv('noNA_sim_bank_const_size_alpha3_Xloci100_v5_no_phase.csv')
# read in observed data
obs_data <- read.csv("../test_scenarios/ABC_combined_obs_data_const_size_pf0_75_alpha3_Xloci100.csv")
print("Read in Data")

## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\



data$name = popname
data <- data[,which(!grepl("Prop_0", colnames(data)))]
simbank_name <- paste("dadi_terms",data$name[1], sep="_")


## ABC-RF Models ####
model.BSR <- readRDS(paste("../models/BSR_model",simbank_name,"no_phase.rds", sep='_'))
model.NC_NA <- readRDS(paste("../models/NC_NA_model",simbank_name,"no_phase.rds", sep='_'))
model.t <- readRDS(paste("../models/t_model",simbank_name,"no_phase.rds", sep='_'))


### Training Data ####
data.BSR <- data.frame(BSR = data[, "beta"], data[, 2:41])
data.NC_NA <- data.frame(NC_NA = data[, 'NC_NA_ratio'], data[,2:41])
data.t <- data.frame(t = data[,'t_NA_ratio'], data[,2:41])



## Output the posterior distribution summaries #####

# reorder obs_data to match colnames of simbank data
obs_data <- obs_data[,names(data)]
obs_data$name = popname
obs_sum_stats <- obs_data[, c(2:41)]
obs_true_vals <- obs_data[,c(42:48)]

print("Obs Sum Stats: ")
print(colnames(obs_sum_stats))
print("True Val Col Names:")
print(colnames(obs_true_vals))

print("Beginning Predictions ")

## Iterate over simulated observed data to predict posterior distribution, each row of observed data is a replicate simulation so we will output posterior estimates for each parameter based on the patterns of variation from each replicate ## 

for (i in 1:nrow(obs_data)) {
  tmpdf <- obs_sum_stats[i, ]
  ## Create Preditions for Each Parameter, also extract 95% CI taking 2.5th and 97.5th quantiles
  pred.obsBSR <- predict(object = model.BSR, obs = tmpdf, training = data.BSR, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
  pred.obst <- predict(object = model.t, obs = tmpdf, training = data.t, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)
  pred.obsNC_NA <- predict(object = model.NC_NA, obs = tmpdf, training = data.NC_NA, quantiles = c(0.025,0.975), paral = TRUE, ncores = 6)

  ## Create Dataframe for Each Parameter
  if (i == 1) {
    BSR_posterior_pred <- as.data.frame(pred.obsBSR)
    t_posterior_pred <- as.data.frame(pred.obst)
    NC_NA_posterior_pred <- as.data.frame(pred.obsNC_NA)
    next
  }
  ## Append predictions for the replicate to the parameter dataframe
  BSR_posterior_pred <- rbind(BSR_posterior_pred, as.data.frame(pred.obsBSR))
  t_posterior_pred <- rbind(t_posterior_pred, as.data.frame(pred.obst))
  NC_NA_posterior_pred <- rbind(NC_NA_posterior_pred, as.data.frame(pred.obsNC_NA))

}


## Save datasets of predictions for each parameter #### 
write.csv(BSR_posterior_pred, file = paste("../posterior_ests/BSR_posterior_pred", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(t_posterior_pred, file = paste("../posterior_ests/t_posterior_pred_", filename, "no_phase.csv", sep = '_'), row.names = F)
write.csv(NC_NA_posterior_pred, file = paste("../posterior_ests/NC_NA_posterior_pred", filename, "no_phase.csv", sep = '_'), row.names = F)


### Compile predicted value and simulated values into one dataframe for output
abc_perform <- data.frame(matrix(nrow = nrow(BSR_posterior_pred), ncol=6))
colnames(abc_perform) <- c("pred_beta", "pred_t", "pred_NC_NA", "sim_beta", "sim_t", "sim_NC_NA")

## popululate output dataframe with expected value of posterior (median)
abc_perform$pred_beta <- BSR_posterior_pred$median
abc_perform$pred_t <- t_posterior_pred$median
abc_perform$pred_NC_NA <- NC_NA_posterior_pred$median
abc_perform$sim_beta <- obs_true_vals$beta
abc_perform$sim_t <- obs_true_vals$t_NA_ratio
abc_perform$sim_NC_NA <- obs_true_vals$NC_NA_ratio

write.csv(abc_perform, file=paste("../posterior_ests/abcrf_", filename, "_no_phase_posterior_perfomance.csv", sep=""), row.names = F)

