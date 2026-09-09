import pandas as pd



num_index = 10000
prefix = "Maple_Ridge"


#######
### Create list of all files for each model
#######

all_beta_fix = []
all_alpha_fix = []
for i in range(num_index):
	all_beta_fix.append(prefix+"_beta_fix_opt_"+str(i)+'.txt')
	all_alpha_fix.append(prefix+"_alpha_fix_opt_"+str(i)+'.txt')


Fixed_Beta_dict = {"popt_fix":[], "ll_opt_fix":[], "theta_opt_fix":[]}
Fixed_Alpha_dict = {"popt_fix":[], "ll_opt_fix":[], "theta_opt_fix":[]}

#######
### Populate Dictionaries for fixed beta and fixed alpha ####
#######

### Iterate over all fixed beta files to populate dictionary

for file in all_beta_fix:

	with open(file, 'r') as f:
		for line in f:
			list_params = line.split("\t")
			popt_fix = list_params[0]
			ll_opt_fix = float(list_params[1])
			theta_opt_fix = float(list_params[2])

	popt_fix_list = []
	for param in popt_fix.split(", "):
		popt_fix_list.append(float(param))


	Fixed_Beta_dict['popt_fix'].append(popt_fix_list)
	Fixed_Beta_dict['ll_opt_fix'].append(ll_opt_fix)
	Fixed_Beta_dict['theta_opt_fix'].append(theta_opt_fix)


### Iterate over all fixed alpha files to populate dictionary

for file in all_alpha_fix:

	with open(file, 'r') as f:
		for line in f: 
			list_params = line.split("\t")
			popt_fix_a = list_params[0]
			ll_opt_fix_a = float(list_params[1])
			theta_opt_fix_a = float(list_params[2])

	popt_fix_list_a = []
	for param in popt_fix_a.split(", "):
		popt_fix_list_a.append(float(param))

	Fixed_Alpha_dict['popt_fix'].append(popt_fix_list_a)
	Fixed_Alpha_dict['ll_opt_fix'].append(ll_opt_fix_a)
	Fixed_Alpha_dict['theta_opt_fix'].append(theta_opt_fix_a)



#### Create a pandas dataframe for each model 

beta_df_dict = {'alpha':[], 'beta':[], 'Tmain':[], 'numain':[], 'Tcurr':[], 'nucurr':[], 'XseqRatio':[], 'll_opt':[], 'theta_opt':[]}
alpha_df_dict = {'alpha':[], 'beta':[], 'Tmain':[], 'numain':[], 'Tcurr':[], 'nucurr':[], 'XseqRatio':[], 'll_opt':[], 'theta_opt':[]}

## populate the dictionaries 
for i in range(len(Fixed_Alpha_dict['ll_opt_fix'])):
	alpha_df_dict['alpha'].append(Fixed_Alpha_dict['popt_fix'][i][0])
	alpha_df_dict['beta'].append(Fixed_Alpha_dict['popt_fix'][i][1])
	alpha_df_dict['Tmain'].append(Fixed_Alpha_dict['popt_fix'][i][2])
	alpha_df_dict['numain'].append(Fixed_Alpha_dict['popt_fix'][i][3])
	alpha_df_dict['Tcurr'].append(Fixed_Alpha_dict['popt_fix'][i][4])
	alpha_df_dict['nucurr'].append(Fixed_Alpha_dict['popt_fix'][i][5])
	alpha_df_dict['XseqRatio'].append(Fixed_Alpha_dict['popt_fix'][i][6])
	alpha_df_dict['ll_opt'].append(Fixed_Alpha_dict['ll_opt_fix'][i])
	alpha_df_dict['theta_opt'].append(Fixed_Alpha_dict['theta_opt_fix'][i])

	beta_df_dict['alpha'].append(Fixed_Beta_dict['popt_fix'][i][0])
	beta_df_dict['beta'].append(Fixed_Beta_dict['popt_fix'][i][1])
	beta_df_dict['Tmain'].append(Fixed_Beta_dict['popt_fix'][i][2])
	beta_df_dict['numain'].append(Fixed_Beta_dict['popt_fix'][i][3])
	beta_df_dict['Tcurr'].append(Fixed_Beta_dict['popt_fix'][i][4])
	beta_df_dict['nucurr'].append(Fixed_Beta_dict['popt_fix'][i][5])
	beta_df_dict['XseqRatio'].append(Fixed_Beta_dict['popt_fix'][i][6])
	beta_df_dict['ll_opt'].append(Fixed_Beta_dict['ll_opt_fix'][i])
	beta_df_dict['theta_opt'].append(Fixed_Beta_dict['theta_opt_fix'][i])


fixed_beta_df = pd.DataFrame(beta_df_dict)
fixed_alpha_df = pd.DataFrame(alpha_df_dict)

sorted_fixed_beta_df = fixed_beta_df.sort_values(by='ll_opt', ascending=False)
sorted_fixed_alpha_df = fixed_alpha_df.sort_values(by='ll_opt', ascending=False)

sorted_fixed_beta_df.to_csv(prefix+"_fixed_beta_opts.csv", index=True)
sorted_fixed_alpha_df.to_csv(prefix+"_fixed_alpha_opts.csv", index=True)











