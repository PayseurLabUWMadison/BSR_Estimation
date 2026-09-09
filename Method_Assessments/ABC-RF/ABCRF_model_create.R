library(abcrf)
start_time <- Sys.time()

## Goal: Create ABC-RF models for each parameter by inputting the reference table constructed using abc_sims_v5.py in parallel on CHTC. 
# All of our explorations of ABC-RF presented in the manuscript focus on models constructed without knowing the phase of the genomic data - we exclude the haplotype statsitics that would be impossible to calculate without resolving phase

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

## read in simbank using notation from abcrf documentation

data_no_phase <- read.csv('SIMBANK_file.csv') ## Read in bank of simulations here 

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

data_no_phase$name = 'UNIQUE Scenario Identifier' ### Edit to reflect simulated criteria

data_no_phase <- data_no_phase[,which(!grepl("Prop_0", colnames(data_no_phase)))] # Scikit-allel will include a bin of 0 in the SFS, which would include invariant sites but we constructed this SFS using variant data so these bins can be removed

simbank_name <- paste("dadi_terms_", data_no_phase$name[1], sep="")


## Create tailored datasets for each parameter of interest ### 

# First, parse each parameter into it's own dataframe with all of the associated simulated summary stats 
data_no_phase.BSR <- data.frame(BSR = data_no_phase[, "beta"], data_no_phase[, 2:41])
data_no_phase.NC_NA <- data.frame(NC_NA = data_no_phase[, 'NC_NA_ratio'], data_no_phase[,2:41])
data_no_phase.t <- data.frame(t = data_no_phase[,'t_NA_ratio'], data_no_phase[,2:41])

## Create RF model for parameters with each forest containing 500 trees that sample 200k replicates with replacement 
# Save model as RDS file to break the process up into seperate processes saving time and ensuring we don't need to run this process again. 
model_no_phase.BSR <- regAbcrf(formula = BSR~., data = data_no_phase.BSR, ntree=500, sampsize=200000, paral=T, ncores=4)
model_no_phase.BSR
saveRDS(model_no_phase.BSR, file = paste("BSR_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.t <- regAbcrf(formula=t~., data = data_no_phase.t, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.t
saveRDS(model_no_phase.t, file = paste("t_model",simbank_name,"no_phase.rds", sep = "_"))

model_no_phase.NC_NA <- regAbcrf(formula=NC_NA~., data = data_no_phase.NC_NA, ntree= 500, sampsize=200000, paral=T, ncores=4)
model_no_phase.NC_NA
saveRDS(model_no_phase.NC_NA, file = paste("NC_NA_model",simbank_name,"no_phase.rds", sep = "_"))

end_time <- Sys.time()

runtime <- end_time - start_time # report runtime for the whole processs 

print("Model Creation for Phase and No Phase Runtime = ")
print(runtime)


