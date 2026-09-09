setwd("~/Maple_Ridge_Inference/dadi/")


###### From windows at least 10kb away from genes  #####

## read in optimized parameters from both dadi models 
MR_10kb_fixed_alpha <- read.csv('MR_10kbgenesfix_fixed_alpha_opts.csv')
MR_10kb_fixed_beta_alpha <- read.csv('MR_10kbgenesfix_fixed_beta_opts.csv')

MR_10kb_fixed_alpha_top100 <- MR_10kb_fixed_alpha[1:100, ]
MR_10kb_fixed_beta_alpha_top100 <- MR_10kb_fixed_beta_alpha[1:100, ]

## Use highest LL values for each model to perform LRT 
# LRT Test stat = -2 * LL(nested) - LL(complex)
test_stat10k = -2 * (MR_10kb_fixed_beta_alpha$ll_opt[1] - MR_10kb_fixed_alpha$ll_opt[1])
pval10k = pchisq(test_stat10k, df=1, lower.tail = F)


###Subsampled Autosome 10kb results ##### 

MR_10kdown_fixed_alpha <- read.csv('MR_10kbgenesfix_downsampledA_fixed_alpha_opts.csv')
MR_10kdown_fixed_beta <- read.csv('MR_10kbgenesfix_downsampledA_fixed_beta_opts.csv')

MR_10kdown_fixed_alpha_top100 <- MR_10kdown_fixed_alpha[1:100, ]
MR_10kdown_fixed_beta_top100 <- MR_10kdown_fixed_beta[1:100, ]

## Use highest LL values for each model to perform LRT 
# LRT Test stat = -2 * LL(nested) - LL(complex)
test_stat10kdown = -2 * (MR_10kdown_fixed_beta$ll_opt[1] - MR_10kdown_fixed_alpha$ll_opt[1])
pval10kdown = pchisq(test_stat10k, df=1, lower.tail = F)


#### Windows at least 5kb away from genes ####
MR_5k_fixed_alpha <- read.csv('MR_5kbgenes_fixed_alpha_opts.csv')
MR_5k_fixed_beta <- read.csv('MR_5kbgenes_fixed_beta_opts.csv')

test_stat_5k <- -2 * (MR_5k_fixed_beta$ll_opt[1] - MR_5k_fixed_alpha$ll_opt[1])
pval5k <- pchisq(test_stat_5k, df = 1, lower.tail = F)

MR_5k_fixed_alpha_top100 <- MR_5k_fixed_alpha[1:100, ]
MR_5k_fixed_beta_top100 <- MR_5k_fixed_beta[1:100, ]


#### Lowered Independence Threshold to 20kb #### 

MR_10k_20kbindep_fixed_alpha <- read.csv('MR_10kbgenes_20kindep_fixed_alpha_opts.csv')
MR_10k_20kbindep_fixed_beta <- read.csv('MR_10kbgenes_20kindep_fixed_beta_opts.csv')

MR_10k_20kbindep_fixed_alpha_top100 <- MR_10k_20kbindep_fixed_alpha[1:100,]
MR_10k_20kbindep_fixed_beta_top100 <- MR_10k_20kbindep_fixed_beta[1:100, ]

test_stat_20kbindep <-  -2 * (MR_10k_20kbindep_fixed_beta$ll_opt[1] - MR_10k_20kbindep_fixed_alpha$ll_opt[1])
pval20kbindep <- pchisq(test_stat_20kbindep, df = 1, lower.tail = F)



### 20kb indep and larger grid search 
MR_10k_20kindep_more_grid_alpha <- read.csv('MR_10kbgenes_20kindep_more_grid_fixed_alpha_opts.csv')
MR_10k_20kindep_more_grid_beta <- read.csv('MR_10kbgenes_20kindep_more_grid_fixed_beta_opts.csv')

## does not increase the concordance between the best fitting parameters 


## Convert the top values for each model into demographic model #### 
mut_rate = 5e-9
length_full_A = 8264095
length_down_A = 3973968
length_5k_A = 7046053
length_20kb_A = 7117308

### Make demographic model in terms of actual sizes #####
#theta = 4*Ne*mut_rate*length_full_A
## Na = theta / (4*mut_rate*length_full_A) 

## 10 kb full 
MR_10kb_fixed_alpha_top100$N_A <- MR_10kb_fixed_alpha_top100$theta_opt / (4*mut_rate*length_full_A)
MR_10kb_fixed_alpha_top100$N_main <- MR_10kb_fixed_alpha_top100$numain * MR_10kb_fixed_alpha_top100$N_A
MR_10kb_fixed_alpha_top100$N_curr <- MR_10kb_fixed_alpha_top100$nucurr * MR_10kb_fixed_alpha_top100$N_A
MR_10kb_fixed_alpha_top100$t_main <- MR_10kb_fixed_alpha_top100$Tmain * (2*MR_10kb_fixed_alpha_top100$N_A)
MR_10kb_fixed_alpha_top100$t_curr <- MR_10kb_fixed_alpha_top100$Tcurr * (2*MR_10kb_fixed_alpha_top100$N_A)

MR_10kb_fixed_beta_alpha_top100$N_A <- MR_10kb_fixed_beta_alpha_top100$theta_opt / (4*mut_rate*length_full_A)
MR_10kb_fixed_beta_alpha_top100$N_main <- MR_10kb_fixed_beta_alpha_top100$numain * MR_10kb_fixed_beta_alpha_top100$N_A
MR_10kb_fixed_beta_alpha_top100$N_curr <- MR_10kb_fixed_beta_alpha_top100$nucurr * MR_10kb_fixed_beta_alpha_top100$N_A
MR_10kb_fixed_beta_alpha_top100$t_main <- MR_10kb_fixed_beta_alpha_top100$Tmain * (2*MR_10kb_fixed_beta_alpha_top100$N_A)
MR_10kb_fixed_beta_alpha_top100$t_curr <- MR_10kb_fixed_beta_alpha_top100$Tcurr * (2*MR_10kb_fixed_beta_alpha_top100$N_A)


## 10kb down 
MR_10kdown_fixed_alpha_top100$N_A <- MR_10kdown_fixed_alpha_top100$theta_opt / (4*mut_rate*length_down_A)
MR_10kdown_fixed_alpha_top100$N_main <- MR_10kdown_fixed_alpha_top100$numain * MR_10kdown_fixed_alpha_top100$N_A
MR_10kdown_fixed_alpha_top100$N_curr <- MR_10kdown_fixed_alpha_top100$nucurr * MR_10kdown_fixed_alpha_top100$N_A
MR_10kdown_fixed_alpha_top100$t_main <- MR_10kdown_fixed_alpha_top100$Tmain * (2*MR_10kdown_fixed_alpha_top100$N_A)
MR_10kdown_fixed_alpha_top100$t_curr <- MR_10kdown_fixed_alpha_top100$Tcurr * (2*MR_10kdown_fixed_alpha_top100$N_A)

MR_10kdown_fixed_beta_top100$N_A <- MR_10kdown_fixed_beta_top100$theta_opt / (4*mut_rate*length_down_A)
MR_10kdown_fixed_beta_top100$N_main <- MR_10kdown_fixed_beta_top100$numain * MR_10kdown_fixed_beta_top100$N_A
MR_10kdown_fixed_beta_top100$N_curr <- MR_10kdown_fixed_beta_top100$nucurr * MR_10kdown_fixed_beta_top100$N_A
MR_10kdown_fixed_beta_top100$t_main <- MR_10kdown_fixed_beta_top100$Tmain * (2*MR_10kdown_fixed_beta_top100$N_A)
MR_10kdown_fixed_beta_top100$t_curr <- MR_10kdown_fixed_beta_top100$Tcurr * (2*MR_10kdown_fixed_beta_top100$N_A)


## 5kb 
MR_5k_fixed_alpha_top100$N_A <- MR_5k_fixed_alpha_top100$theta_opt / (4*mut_rate*length_5k_A)
MR_5k_fixed_alpha_top100$N_main <- MR_5k_fixed_alpha_top100$numain * MR_5k_fixed_alpha_top100$N_A
MR_5k_fixed_alpha_top100$N_curr <- MR_5k_fixed_alpha_top100$nucurr * MR_5k_fixed_alpha_top100$N_A
MR_5k_fixed_alpha_top100$t_main <- MR_5k_fixed_alpha_top100$Tmain * (2*MR_5k_fixed_alpha_top100$N_A)
MR_5k_fixed_alpha_top100$t_curr <- MR_5k_fixed_alpha_top100$Tcurr * (2*MR_5k_fixed_alpha_top100$N_A)

MR_5k_fixed_beta_top100$N_A <- MR_5k_fixed_beta_top100$theta_opt / (4*mut_rate*length_5k_A)
MR_5k_fixed_beta_top100$N_main <- MR_5k_fixed_beta_top100$numain * MR_5k_fixed_beta_top100$N_A
MR_5k_fixed_beta_top100$N_curr <- MR_5k_fixed_beta_top100$nucurr * MR_5k_fixed_beta_top100$N_A
MR_5k_fixed_beta_top100$t_main <- MR_5k_fixed_beta_top100$Tmain * (2*MR_5k_fixed_beta_top100$N_A)
MR_5k_fixed_beta_top100$t_curr <- MR_5k_fixed_beta_top100$Tcurr * (2*MR_5k_fixed_beta_top100$N_A)

## 10kb_20kbindep
MR_10k_20kbindep_fixed_alpha_top100$N_A <- MR_10k_20kbindep_fixed_alpha_top100$theta_opt / (4*mut_rate*length_20kb_A)
MR_10k_20kbindep_fixed_alpha_top100$N_main <- MR_10k_20kbindep_fixed_alpha_top100$numain * MR_10k_20kbindep_fixed_alpha_top100$N_A
MR_10k_20kbindep_fixed_alpha_top100$N_curr <- MR_10k_20kbindep_fixed_alpha_top100$nucurr * MR_10k_20kbindep_fixed_alpha_top100$N_A
MR_10k_20kbindep_fixed_alpha_top100$t_main <- MR_10k_20kbindep_fixed_alpha_top100$Tmain * (2*MR_10k_20kbindep_fixed_alpha_top100$N_A)
MR_10k_20kbindep_fixed_alpha_top100$t_curr <- MR_10k_20kbindep_fixed_alpha_top100$Tcurr * (2*MR_10k_20kbindep_fixed_alpha_top100$N_A)

MR_10k_20kbindep_fixed_beta_top100$N_A <- MR_10k_20kbindep_fixed_beta_top100$theta_opt / (4*mut_rate*length_20kb_A)
MR_10k_20kbindep_fixed_beta_top100$N_main <- MR_10k_20kbindep_fixed_beta_top100$numain * MR_10k_20kbindep_fixed_beta_top100$N_A
MR_10k_20kbindep_fixed_beta_top100$N_curr <- MR_10k_20kbindep_fixed_beta_top100$nucurr * MR_10k_20kbindep_fixed_beta_top100$N_A
MR_10k_20kbindep_fixed_beta_top100$t_main <- MR_10k_20kbindep_fixed_beta_top100$Tmain * (2*MR_10k_20kbindep_fixed_beta_top100$N_A)
MR_10k_20kbindep_fixed_beta_top100$t_curr <- MR_10k_20kbindep_fixed_beta_top100$Tcurr * (2*MR_10k_20kbindep_fixed_beta_top100$N_A)


## write csv for top 100 lists 

write.csv(MR_10kb_fixed_alpha_top100, file = 'MR_10kb_fixed_alpha_top100_inferences_corrected.csv', row.names = F)
write.csv(MR_10kb_fixed_beta_alpha_top100, file = "MR_10kb_fixed_beta_top100_inferences_corrected.csv", row.names = F)

write.csv(MR_10kdown_fixed_alpha_top100, file = 'MR_10kbdown_fixed_alpha_top100_inferences_corrected.csv', row.names = F)
write.csv(MR_10kdown_fixed_beta_top100, file = "MR_10kbdown_fixed_beta_top100_inferences_corrected.csv", row.names = F)

write.csv(MR_5k_fixed_alpha_top100, file = 'MR_5kb_fixed_alpha_top100_inferences_corrected.csv', row.names = F)
write.csv(MR_5k_fixed_beta_top100, file = "MR_5kb_fixed_beta_top100_inferences_corrected.csv", row.names = F)

write.csv(MR_10k_20kbindep_fixed_alpha_top100, file = 'MR_10kb_20kbindep_fixed_alpha_top100_inferences_corrected.csv', row.names = F)
write.csv(MR_10k_20kbindep_fixed_beta_top100, file = "MR_10kb_20kbindep_fixed_beta_top100_inferences_corrected.csv", row.names = F)


