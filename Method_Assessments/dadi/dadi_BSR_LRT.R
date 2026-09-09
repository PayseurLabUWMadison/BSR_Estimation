## Goal: Perform LRT for scenarios where simulated beta ≠ 1 using fixed alpha model and nested model where alpha and beta are fixed 
library(ggplot2)
library(patchwork)

## Read in all simulated scenarios with a BSR that is not 1:1 ## 

## used for labelling plots
scenarios <- c("Constant Size, Beta = 1/3", "Constant Size, Beta = 3", "95% Reduction, Beta = 1/3", "95% Reduction, Beta = 3", "20x Expansion, Beta = 1/3", "20x Expansion, Beta = 3", "Constant Size, Beta = 1/3, Alpha = 3", "Constant Size, Beta = 3, Alpha = 3", "95% Reduction, Beta = 1/3, Alpha = 3", "95% Reduction, Beta = 3, Alpha = 3", "20x Expansion, Beta = 1/3, Alpha = 3", "20x Expansion, Beta = 3, Alpha = 3")

# used for labelling files 
prefixes <- c('const_a1_b03', 'const_a1_b3', 'red_a1_b03', 'red_a1_b3', 'exp_a1_b03', 'exp_a1_b3', 'const_a3_b03', 'const_a3_b3', 'red_a3_b03', 'red_a3_b3', 'exp_a3_b03', 'exp_a3_b3')

# const_size alpha 1
const_a1_b03_alpha <- read.csv('dadi_combined_const_size_pf0_25_alpha1_Xloci100_fixed_alpha_opts.csv')
const_a1_b03_beta <- read.csv('dadi_combined_const_size_pf0_25_alpha1_Xloci100_fixed_beta_opts.csv')
const_a1_b3_alpha <- read.csv('dadi_combined_const_size_pf0_75_alpha1_Xloci100_fixed_alpha_opts.csv')
const_a1_b3_beta <- read.csv('dadi_combined_const_size_pf0_75_alpha1_Xloci100_fixed_beta_opts.csv')

# pop red alpha 1
red_a1_b03_alpha <- read.csv('dadi_combined_pop_red_pf0_25_alpha1_Xloci100_fixed_alpha_opts.csv')
red_a1_b03_beta <- read.csv('dadi_combined_pop_red_pf0_25_alpha1_Xloci100_fixed_alphabeta_opts.csv')
red_a1_b3_alpha <- read.csv('dadi_combined_pop_red_pf0_75_alpha1_Xloci100_fixed_alpha_opts.csv')
red_a1_b3_beta <- read.csv('dadi_combined_pop_red_pf0_75_alpha1_Xloci100_fixed_alphabeta_opts.csv')

# pop exp alpha 1 
exp_a1_b03_alpha <- read.csv('dadi_combined_pop_exp_pf0_25_alpha1_Xloci100_fixed_alpha_opts.csv')
exp_a1_b03_beta <- read.csv('dadi_combined_pop_exp_pf0_25_alpha1_Xloci100_fixed_alphabeta_opts.csv')
exp_a1_b3_alpha <- read.csv('dadi_combined_pop_exp_pf0_75_alpha1_Xloci100_fixed_alpha_opts.csv')
exp_a1_b3_beta <- read.csv('dadi_combined_pop_exp_pf0_75_alpha1_Xloci100_fixed_alphabeta_opts.csv')

# const_size alpha 3
const_a3_b03_alpha <- read.csv('dadi_combined_const_size_pf0_25_alpha3_Xloci100_fixed_alpha_opts.csv')
const_a3_b03_beta <- read.csv('dadi_combined_const_size_pf0_25_alpha3_Xloci100_fixed_beta_opts.csv')
const_a3_b3_alpha <- read.csv('dadi_combined_const_size_pf0_75_alpha3_Xloci100_fixed_alpha_opts.csv')
const_a3_b3_beta <- read.csv('dadi_combined_const_size_pf0_75_alpha3_Xloci100_fixed_beta_opts.csv')

# pop red alpha 3
red_a3_b03_alpha <- read.csv('dadi_combined_pop_red_pf0_25_alpha3_Xloci100_fixed_alpha_opts.csv')
red_a3_b03_beta <- read.csv('dadi_combined_pop_red_pf0_25_alpha3_Xloci100_fixed_alphabeta_opts.csv')
red_a3_b3_alpha <- read.csv('dadi_combined_pop_red_pf0_75_alpha3_Xloci100_fixed_alpha_opts.csv')
red_a3_b3_beta <- read.csv('dadi_combined_pop_red_pf0_75_alpha3_Xloci100_fixed_alphabeta_opts.csv')

# pop exp alpha 3 
exp_a3_b03_alpha <- read.csv('dadi_combined_pop_exp_pf0_25_alpha3_Xloci100_fixed_alpha_opts.csv')
exp_a3_b03_beta <- read.csv('dadi_combined_pop_exp_pf0_25_alpha3_Xloci100_fixed_alphabeta_opts.csv')
exp_a3_b3_alpha <- read.csv('dadi_combined_pop_exp_pf0_75_alpha3_Xloci100_fixed_alpha_opts.csv')
exp_a3_b3_beta <- read.csv('dadi_combined_pop_exp_pf0_75_alpha3_Xloci100_fixed_alphabeta_opts.csv')

## Create a list of all scenarios 
fix_alphas <- list(const_a1_b03_alpha, const_a1_b3_alpha, red_a1_b03_alpha, red_a1_b3_alpha, exp_a1_b03_alpha, exp_a1_b3_alpha, const_a3_b03_alpha, const_a3_b3_alpha, red_a3_b03_alpha, red_a3_b3_alpha, exp_a3_b03_alpha, exp_a3_b3_alpha)
fix_betas <- list(const_a1_b03_beta, const_a1_b3_beta, red_a1_b03_beta, red_a1_b3_beta, exp_a1_b03_beta, exp_a1_b3_beta, const_a3_b03_beta, const_a3_b3_beta, red_a3_b03_beta, red_a3_b3_beta, exp_a3_b03_beta, exp_a3_b3_beta)

## Calculate the test statistic for each scenario using the equation: -2*(LL(nested) - LL(complex))
test_stats <- list(NA)
## LRT Test stat = -2 * LL(nested) - LL(complex)
for (i in c(1:length(fix_alphas))){
  test_stats[[i]] <- data.frame(matrix(nrow = nrow(const_a1_b03_alpha), ncol=0))
  test_stats[[i]]$test_stat <- -2 * (fix_betas[[i]]$ll_opt_fix - fix_alphas[[i]]$ll_opt_fix)
  test_stats[[i]]$scenario <- scenarios[i]
  ## Chi square test with one degree of freedom to determine p-value 
  test_stats[[i]]$pval <- pchisq(test_stats[[i]]$test_stat, df=1, lower.tail = F)
}


### Iterate over each sceanrio to create plots of distribution for all p-values
pval_plots <- list(NA)

for (i in c(1:length(test_stats))) {
  pval_plots[[i]] <- ggplot(test_stats[[i]], aes(x=pval))+
  geom_histogram(fill='black')+
  ggtitle(paste(test_stats[[i]]$scenario[1], "LRT Results"))+
  geom_vline(xintercept = 0.05, linetype = 'dashed')+
  theme_bw()+
  xlab("P-value")+
  theme(strip.text = element_text(size=16), axis.text = element_text(size = 15), axis.title = element_text(size = 18), legend.text = element_text(size = 16), legend.title = element_text(size=19))
}


for (i in c(1:length(pval_plots))){
  ggsave(filename = paste('../plots/',prefixes[i], '_LRT_pvals.png', sep = ""), pval_plots[[i]], width = 6, height = 6, units = 'in')
}


## Split up the plots for each demographic model 
const <- pval_plots[[1]] + pval_plots[[2]]
red <- pval_plots[[3]] + pval_plots[[4]]
exp <- pval_plots[[5]] + pval_plots[[6]]

ggsave(filename ='../plots/const_a1_LRT_pvals.png', const, width = 12, height = 4, units = 'in')
ggsave(filename ='../plots/red_a1_LRT_pvals.png', red, width = 12, height = 4, units = 'in')
ggsave(filename ='../plots/exp_a1_LRT_pvals.png', exp, width = 12, height = 4, units = 'in')


## how many replicates in each scenario have significant p-value? - make the table describing power to infer BSR 
dadi_power_df <- data.frame(matrix(ncol = length(test_stats), nrow=1))
colnames(dadi_power_df) <- prefixes

for (i in c(1:length(test_stats))){
  dadi_power_df[,i] = length(which(test_stats[[i]]$pval < 0.05)) / 100
}

write.csv(dadi_power_df, file = 'dadi_BSR_power_table.csv', row.names = F)
