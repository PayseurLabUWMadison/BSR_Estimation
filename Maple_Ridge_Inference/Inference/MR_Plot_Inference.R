library(ggplot2)
library(patchwork)
library(reshape2)

setwd("~/Maple_Ridge_Inference/")

##### \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
### ABC-RF INFERENCE #####

# edit to change scenario of focus 
prefix = '10kbgenes_20kindep_alpha2_Xloci1319'

# read in posterior estimates for each parameter 

abc_bsr <- read.csv(paste('ABC_RF/posterior_ests/BSR_posterior_pred_Maple_Ridge_',prefix,'_no_phase.csv', sep=''))
abc_bsr$parameter <- "Beta"
abc_NM_NA <- read.csv(paste('ABC_RF/posterior_ests/NM_NA_posterior_pred_Maple_Ridge_',prefix,'_no_phase.csv', sep=''))
abc_NM_NA$parameter <- "nu1"
abc_NC_NA <- read.csv(paste('ABC_RF/posterior_ests/NC_NA_posterior_pred_Maple_Ridge_',prefix,'_no_phase.csv', sep=''))
abc_NC_NA$parameter <- 'nu2'
abc_t_main <- read.csv(paste('ABC_RF/posterior_ests/t_main_posterior_pred__Maple_Ridge_',prefix,'_no_phase.csv', sep=''))
abc_t_main$parameter <- 'T1'
abc_t_curr <- read.csv(paste('ABC_RF/posterior_ests/t_curr_posterior_pred__Maple_Ridge_',prefix,'_no_phase.csv', sep=''))
abc_t_curr$parameter <- 'T2'

abc_plot_df <- rbind(abc_bsr, abc_NM_NA)
abc_plot_df <- rbind(abc_plot_df, abc_NC_NA)
abc_plot_df <- rbind(abc_plot_df, abc_t_main)
abc_plot_df <- rbind(abc_plot_df, abc_t_curr)

abc_plot_df$Method <- "ABC"


## Plot the point estimates 

abc_post_bsr <- ggplot(abc_bsr, aes(y=median, x=1))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = quantile.0.025, ymax = quantile.0.975), width = 0.25)+
  ylab('Beta (Nf:Nm)')+
  theme_bw()+
  xlab("Maple Ridge Inference")+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 0.75, alpha = 0)+
  geom_vline(xintercept = 1.25, alpha = 0)+
  ggtitle("Beta Estimate")

abc_post_NM_NA <- ggplot(abc_NM_NA, aes(y=median, x=1))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = quantile.0.025, ymax = quantile.0.975), width = 0.25)+
  ylab('Ratio of Mainland Ancestor and Ancestral Size')+
  theme_bw()+
  xlab("Maple Ridge Inference")+
  geom_vline(xintercept = 0.75, alpha = 0)+
  geom_vline(xintercept = 1.25, alpha = 0)+
  ggtitle("NM:NA Estimate")

abc_post_NC_NA <- ggplot(abc_NC_NA, aes(y=median, x=1))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = quantile.0.025, ymax = quantile.0.975), width = 0.25)+
  ylab('Ratio of Current and Ancestral Size')+
  theme_bw()+
  xlab("Maple Ridge Inference")+
  geom_vline(xintercept = 0.75, alpha = 0)+
  geom_vline(xintercept = 1.25, alpha = 0)+
  ggtitle("NC:NA Estimate")

abc_post_t_curr <- ggplot(abc_t_curr, aes(y=median, x=1))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = quantile.0.025, ymax = quantile.0.975), width = 0.25)+
  ylab('Mainland Anc -> Current Timing (2Na Generations)')+
  theme_bw()+
  xlab("Maple Ridge Inference")+
  geom_vline(xintercept = 0.75, alpha = 0)+
  geom_vline(xintercept = 1.25, alpha = 0)+
  ggtitle("Timing of NM -> NC Est.")


abc_post_t_main <- ggplot(abc_t_main, aes(y=median, x=1))+
  geom_point(size=2)+
  geom_errorbar(aes(ymin = quantile.0.025, ymax = quantile.0.975), width = 0.25)+
  ylab('Ancestral -> Mainland Anc Timing (2Na Generations)')+
  theme_bw()+
  xlab("Maple Ridge Inference")+
  geom_vline(xintercept = 0.75, alpha = 0)+
  geom_vline(xintercept = 1.25, alpha = 0)+
  ggtitle("Timing of NA -> NM Est.")

layout <- '
#AA#
BBCC
DDEE'
point_ests <- abc_post_bsr+ abc_post_NM_NA+ abc_post_NC_NA+ abc_post_t_curr+ abc_post_t_main+ plot_layout(design = layout)

ggsave(filename = paste('../plots/MR_',prefix,'_ABC_point_ests.png', sep=''), point_ests, width = 15, height = 6, units = 'in')

## \/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
### Dadi Inference #### 

dadi_inf <- read.csv('dadi/MR_10kb_20kindep_dadi_inferences.csv')
dadi_inf$Method <- "dadi"



## make the ABC and dadi plots agree 
ABC_comb_df <- abc_plot_df[c('median', 'quantile.0.025', 'quantile.0.975', 'parameter', 'Method') ]
colnames(ABC_comb_df) <- c('Estimate', 'Lower_bound', 'Upper_bound', 'Parameter', 'Method')

dadi_comb_df <- dadi_inf[2:6,c('Estimate', 'Lower.Limit_2SD', 'Upper.Limit_2SD', 'Parameter', 'Method')]
colnames(dadi_comb_df) <- c("Estimate", 'Lower_bound', 'Upper_bound', "Parameter", 'Method')

MR_comb_df <- rbind(ABC_comb_df, dadi_comb_df)

## Facet plot of each parameter inference
ggplot(MR_comb_df, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2)+
  facet_wrap(Parameter~., scales = 'free_y', nrow=3, ncol = 2, strip.position = 'top')+
  geom_errorbar(aes(ymin = Lower_bound, ymax = Upper_bound), width = 0.25)+
  theme_bw()

## Plot each parameter seperately
beta_plot <- MR_comb_df[which(grepl("Beta", MR_comb_df$Parameter)), ]
b <- ggplot(beta_plot, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=Lower_bound, ymax=Upper_bound), width=0.25)+
  ylab("Beta (Nf:Nm)")+
  theme_bw()+
  geom_hline(yintercept = 1, linetype='dashed')+
  scale_color_manual(values = c('#F8766D', '#00BA38'), 'Method')

nu1_plot <- MR_comb_df[which(grepl("nu1", MR_comb_df$Parameter)), ]
n1 <- ggplot(nu1_plot, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=Lower_bound, ymax=Upper_bound), width=0.25)+
  ylab("Size Change Epoch 1")+
  theme_bw()+
  scale_color_manual(values = c('#F8766D', '#00BA38'), 'Method')

nu2_plot <- MR_comb_df[which(grepl("nu2", MR_comb_df$Parameter)), ]
n2 <- ggplot(nu2_plot, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=Lower_bound, ymax=Upper_bound), width=0.25)+
  ylab("Size Change Epoch 2")+
  theme_bw()+
  scale_color_manual(values = c('#F8766D', '#00BA38'), 'Method')

T1_plot <- MR_comb_df[which(grepl("T1", MR_comb_df$Parameter)), ]
t1 <- ggplot(T1_plot, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=Lower_bound, ymax=Upper_bound), width=0.25)+
  ylab("Timing Epoch 1")+
  theme_bw()+
  scale_color_manual(values = c('#F8766D', '#00BA38'), 'Method')

T2_plot <- MR_comb_df[which(grepl("T2", MR_comb_df$Parameter)), ]
t2 <- ggplot(T2_plot, aes(x=Method, y=Estimate, color=Method))+
  geom_point(size=2.5)+
  geom_errorbar(aes(ymin=Lower_bound, ymax=Upper_bound), width=0.25)+
  ylab("Timing Epoch 2")+
  theme_bw()+
  scale_color_manual(values = c('#F8766D', '#00BA38'), 'Method')



layout <- "
#AA#
BBCC
DDEE
"

MR_out <- b + t1 + n1 + t2 + n2 + plot_layout(guides = 'collect', design = layout)+ plot_annotation(tag_levels = "A", title = 'Maple Ridge Inferred Parameters', theme=theme(plot.title = element_text(hjust = 0.5, size = 16),))
MR_out


ggsave('Maple_Ridge_Inferred_Params.png', MR_out, width = 8, height=7, units = 'in')
