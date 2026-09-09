setwd('/mnt/payseurgenomes/will/Maple_Ridge_Inference/simbank')

library(abcrf)
library(ggplot2)
library(patchwork)

## for troubleshooting
# setwd("~/ABC_Eval/simulations/")
# data <- read.csv("sim_bank_pop_red_new_alpha1_pilot.csv")


### Data with no phase known####

data_no_phase <- read.csv('noNA_simbank_Maple_Ridge_10kbgenes_20kindep_alpha2_Xloci1319_v5_no_phase.csv')
data_no_phase <- data_no_phase[,which(!grepl("Prop_0", colnames(data_no_phase)))]

simbank_name <- "Maple_Ridge_10kbgenes_20kindep"

 # Create tailored datasets for each parameter of interest (need separate simulations for each parameter in practice)
data_no_phase.BSR <- data.frame(BSR = data_no_phase[, "beta"], data_no_phase[, 2:50])
data_no_phase.NC_NA <- data.frame(NC_NA = data_no_phase[, 'NC_NA_ratio'], data_no_phase[,2:50])
data_no_phase.NM_NA <- data.frame(NM_NA = data_no_phase[,'Nmain_NA_ratio'], data_no_phase[,2:50])
data_no_phase.t_curr <- data.frame(t_curr = data_no_phase[,'t_curr_NA_ratio'], data_no_phase[,2:50])
data_no_phase.t_main <- data.frame(t_main = data_no_phase[,'t_main_NA_ratio'], data_no_phase[,2:50])

##  Read in models for each parameter (dadi terms)
model_no_phase.BSR <- readRDS(paste("../models/BSR_model_dadi_terms",simbank_name,"no_phase.rds", sep='_'))
model_no_phase.NC_NA <- readRDS(paste("../models/NC_NA_model_dadi_terms",simbank_name,"no_phase.rds", sep='_'))
model_no_phase.NM_NA <- readRDS(paste("../models/NM_NA_model_dadi_terms",simbank_name,"no_phase.rds", sep='_'))
model_no_phase.t_curr <- readRDS(paste("../models/t_curr_model_dadi_terms",simbank_name,"no_phase.rds", sep='_'))
model_no_phase.t_main <- readRDS(paste("../models/t_main_model_dadi_terms",simbank_name,"no_phase.rds", sep='_'))


## Plot out of bag simulations compared to predicted values for those sims in model
png(filename = paste("../plots/",simbank_name,"_dadi_terms_no_phase_oob_estimates_BSR.png",sep=""), width = 1400, height = 1200, res = 200)
plot(data_no_phase$beta,
      model_no_phase.BSR$model.rf$predictions,
      xlab=expression('beta'),
      ylab=expression("out-of-bag "*hat(beta)), main = "BSR Out of Bag Estimates")
abline(a=0,b=1,col="red")
dev.off()

png(filename = paste("../plots/",simbank_name,"_dadi_terms_no_phase_oob_estimates_NC_NA.png",sep=""), width = 1400, height = 1200, res = 200)
plot(data_no_phase$NC_NA_ratio,
      model_no_phase.NC_NA$model.rf$predictions,
      xlab=expression('NC:NA'),
      ylab=expression("out-of-bag "*hat("NC:NA")),
      main = "Current Size : Ancestral Size of Bag Estimates")
abline(a=0,b=1,col="red")
dev.off()

png(filename = paste("../plots/",simbank_name,"_dadi_terms_no_phase_oob_estimates_NM_NA.png",sep=""), width = 1400, height = 1200, res = 200)
plot(data_no_phase$Nmain_NA_ratio,
      model_no_phase.NM_NA$model.rf$predictions,
      xlab=expression('NM:NA'),
      ylab=expression("out-of-bag "*hat('NM:NA')), main = "Mainland Anc : Ancestral Size Out of Bag Estimates")
abline(a=0,b=1,col="red")
dev.off()

png(filename = paste("../plots/",simbank_name,"_dadi_terms_no_phase_oob_estimates_t_curr.png",sep=""), width = 1400, height = 1200, res = 200)
plot(data_no_phase$t_curr_NA_ratio,
      model_no_phase.t_curr$model.rf$predictions,
      xlab=expression('time'),
      ylab=expression("out-of-bag "*hat(time)), main = "Mainland -> Current Time Change Out of Bag Estimates")
abline(a=0,b=1,col="red")
dev.off()

png(filename = paste("../plots/",simbank_name,"_dadi_terms_no_phase_oob_estimates_t_main.png",sep=""), width = 1400, height = 1200, res = 200)
plot(data_no_phase$t_main_NA_ratio,
     model_no_phase.t_main$model.rf$predictions,
     xlab=expression('time'),
     ylab=expression("out-of-bag "*hat(time)), main = "Ancestral -> Mainland Time Change Out of Bag Estimates")
abline(a=0,b=1,col="red")
dev.off()

### Plot the out-of-bag mean sq error for number of trees for each parameter model
png(filename = paste("../plots/",simbank_name,"_no_phase_errorOOB_BSR.png", sep = ""), width = 1400, height = 1200, res = 200)
errorOOB_BSR <- err.regAbcrf(object = model_no_phase.BSR, training = data_no_phase.BSR)
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_errorOOB_NC_NA.png",sep=""), width = 1400, height = 1200, res = 200)
errorOOB_N_A <- err.regAbcrf(object=model_no_phase.NC_NA, training = data_no_phase.NC_NA)
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_errorOOB_NM_NA.png", sep=""), width = 1400, height = 1200, res = 200)
errorOOB_N_C <- err.regAbcrf(object = model_no_phase.NM_NA, training = data_no_phase.NM_NA)
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_errorOOB_t_curr.png", sep=""), width = 1400, height = 1200, res = 200)
errorOOB_t <- err.regAbcrf(object = model_no_phase.t_curr, training = data_no_phase.t_curr)
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_errorOOB_t_main.png", sep=""), width = 1400, height = 1200, res = 200)
errorOOB_t <- err.regAbcrf(object = model_no_phase.t_main, training = data_no_phase.t_main)
dev.off()


## What are the most useful summary statistics? - Stats with higher importance on X axis are most meaningful
png(filename = paste("../plots/",simbank_name,"_no_phase_BSR_sumstat_importance.png",sep=""), width = 1400, height = 1200, res = 200)
plot(x=model_no_phase.BSR, n.var = 25, main = "Variable Importance for BSR")
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_NC_NA_sumstat_importance.png", sep=""), width = 1400, height = 1200, res = 200)
plot(x=model_no_phase.NC_NA, n.var = 25, main = "Variable Importance for Epoch 2 Size Change")
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_NM_NA_sumstat_importance.png",sep=""), width = 1400, height = 1200, res = 200)
plot(x=model_no_phase.NM_NA, n.var = 25, main = "Variable Importance for Epoch 1 Size Change")
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_t_curr_sumstat_importance.png", sep=""), width = 1400, height = 1200, res = 200)
plot(x=model_no_phase.t_curr, n.var= 25, main = "Variable Importance for Epoch 2 Timing")
dev.off()

png(filename = paste("../plots/",simbank_name,"_no_phase_t_main_sumstat_importance.png", sep=""), width = 1400, height = 1200, res = 200)
plot(x=model_no_phase.t_main, n.var= 25, main = "Variable Importance for Epoch 1 Timing")
dev.off()






