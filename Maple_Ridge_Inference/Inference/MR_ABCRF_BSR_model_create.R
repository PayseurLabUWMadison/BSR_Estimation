setwd('/mnt/payseurgenomes/will/Maple_Ridge_Inference/simbank/')
library(abcrf)
start_time <- Sys.time()

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

## read in simbank using notation from abcrf documentation

data <- read.csv('fixed_simbank_Maple_Ridge_alpha2_Xloci154_v5.csv')

data <- na.omit(data)
data <- data[,which(!grepl("Prop_0", colnames(data)))]

data$name = 'Maple_Ridge'

write.csv(data, 'noNA_simbank_Maple_Ridge_alpha2_Xloci154_v5.csv', row.names = FALSE)

## Create no phase dataset 
### RUN THIS ONLY WHEN USING SIMBANK FOR THE FIRST TIME 
# data_no_phase <- data[, which(!grepl("M_", colnames(data)))]
# data_no_phase <- data_no_phase[, which(!grepl("K_", colnames(data_no_phase)))]
# write.csv(data_no_phase, 'noNA_simbank_Maple_Ridge_alpha2_Xloci154_v5_no_phase.csv', row.names = FALSE)
data_no_phase <- read.csv('noNA_simbank_Maple_Ridge_alpha2_Xloci154_v5_no_phase.csv')

print("Read in Data")

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

## for troubleshooting
# setwd("~/ABC_Eval/simulations/")
# data_no_phase <- read.csv("~/ABC_Eval/simulations/noNA_pilot_Xloci100_red_a1_no_phase.csv")
# data_no_phase <- data_no_phase[,which(!grepl("Prop_0", colnames(data_no_phase)))]
# data <- read.csv('sim_bank_pop_red_new_alpha1_v5.csv')


simbank_name <- paste("dadi_terms_", data$name[1], sep="")

### Data with Phase Known ####
## Create tailored datasets for each parameter of interest (need separate simulations for each parameter in practice)
data.BSR <- data.frame(BSR = data[, "beta"], data[, 2:58])
data.NC_NA <- data.frame(NC_NA = data[, 'NC_NA_ratio'], data[,2:58])
data.NM_NA <- data.frame(NM_NA = data[, 'Nmain_NA_ratio'], data[,2:58])
data.t_curr <- data.frame(t_curr = data[, 't_curr_NA_ratio'], data[,2:58])
data.t_main <- data.frame(t_main = data[, 't_main_NA_ratio'], data[,2:58])

print(colnames(data.BSR))

##SAVE MODELS TO ALLOW FOR EASIER EXPLORATION####
#saveRDS(model.BSR, file = "BSR_abcrf_model_phase_pop_red_a1.rds")
#####

## Create RF model for parameters
# Saves model as object in terminal so I don't need to regenerate it each time
model.BSR <- regAbcrf(formula = BSR~., data = data.BSR, ntree=500, paral=T, ncores=4, sampsize = 200000)
model.BSR
saveRDS(model.BSR, file = paste("../models/BSR_model",simbank_name,"phase.rds", sep = "_"))

model.t_curr <- regAbcrf(formula=t_curr~., data=data.t_curr, ntree= 500, paral=T, ncores=4, sampsize = 200000)
model.t_curr
saveRDS(model.t_curr, file = paste("../models/t_curr_model",simbank_name,"phase.rds", sep = "_"))

model.t_main <- regAbcrf(formula=t_main~., data=data.t_main, ntree= 500, paral=T, ncores=4, sampsize = 200000)
model.t_main
saveRDS(model.t_main, file = paste("../models/t_main_model",simbank_name,"phase.rds", sep = "_"))

model.NC_NA <- regAbcrf(formula = NC_NA~., data=data.NC_NA, ntree=500, paral=T, ncores=4, sampsize = 200000)
model.NC_NA
saveRDS(model.NC_NA, file=paste("../models/NC_NA_model",simbank_name,"phase.rds", sep = "_"))

model.NM_NA <- regAbcrf(formula = NM_NA~., data=data.NM_NA, ntree=500, paral=T, ncores=4, sampsize = 200000)
model.NM_NA
saveRDS(model.NM_NA, file=paste("../models/NM_NA_model",simbank_name,"phase.rds", sep = "_"))



## Data with no phase known####

## Create tailored datasets for each parameter of interest (need separate simulations for each parameter in practice)
data_no_phase.BSR <- data.frame(BSR = data_no_phase[, "beta"], data_no_phase[, 2:50])
data_no_phase.NC_NA <- data.frame(NC_NA = data_no_phase[, 'NC_NA_ratio'], data_no_phase[,2:50])
data_no_phase.NM_NA <- data.frame(NM_NA = data_no_phase[, 'Nmain_NA_ratio'], data_no_phase[,2:50])
data_no_phase.t_curr <- data.frame(t_curr = data_no_phase[, 't_curr_NA_ratio'], data_no_phase[,2:50])
data_no_phase.t_main <- data.frame(t_main = data_no_phase[, 't_main_NA_ratio'], data_no_phase[,2:50])

print(colnames(data_no_phase.BSR))

## Create RF model for parameters
model_no_phase.BSR <- regAbcrf(formula = BSR~., data = data_no_phase.BSR, ntree=500, sampsize=200000, paral=T, ncores=4)
model_no_phase.BSR
saveRDS(model_no_phase.BSR, file = paste("../models/BSR_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.t_curr <- regAbcrf(formula=t_curr~., data = data_no_phase.t_curr, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.t_curr
saveRDS(model_no_phase.t_curr, file = paste("../models/t_curr_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.t_main <- regAbcrf(formula=t_main~., data = data_no_phase.t_main, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.t_main
saveRDS(model_no_phase.t_main, file = paste("../models/t_main_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.NM_NA <- regAbcrf(formula=NM_NA~., data = data_no_phase.NM_NA, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.NM_NA
saveRDS(model_no_phase.NM_NA, file = paste("../models/NM_NA_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.NC_NA <- regAbcrf(formula=NC_NA~., data = data_no_phase.NC_NA, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.NC_NA
saveRDS(model_no_phase.NC_NA, file = paste("../models/NC_NA_model",simbank_name,"no_phase.rds", sep = "_"))

end_time <- Sys.time()

runtime <- end_time - start_time

print("Model Creation for Phase and No Phase Runtime = ")
print(runtime)


