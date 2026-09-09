import pandas as pd 

n_SFS = 100

prefix = 'dadi_combined_pop_red_pf0_5_alpha1_Xloci100_folded'


free_param_out_dict = {"alpha":[], "beta":[], "T":[], 'nu':[], "ll_opt":[], "theta_opt":[]}
fixed_beta_out_dict = {"alpha_fix":[], "beta_fix":[], "T_fix":[], 'nu_fix':[], "ll_opt_fix":[], "theta_opt_fix":[]}
fixed_alpha_out_dict = {"alpha_fix":[], "beta_fix":[], "T_fix":[], 'nu_fix':[], "ll_opt_fix":[], "theta_opt_fix":[]}



def get_best_fit(params):
    """
    Extract the list of [parameters, log-likelihood, theta] that best fit the data after optimization
    """
    best_ll = float('-inf')
    best_fit_index = 0
    for i, rep in enumerate(params):
        if rep > best_ll:
            best_fit = rep
            best_ll = rep
            best_fit_index = i

    return best_fit_index, best_fit



for rep in range(n_SFS):
	########
	### SPECIFY PARAMETERS
	######## 

	num_index = 100

	free_file_head = str(rep)+"_"+prefix+"_free_param_opt"
	fixed_beta_head = str(rep)+"_"+prefix+"_beta_fix_opt"
	fixed_alpha_head = str(rep)+"_"+prefix+"_alpha_fix_opt"

	########
	### Create list of file names 
	########

	all_free_files = []
	all_beta_fix = []
	all_alpha_fix = []

	for i in range(num_index):
		all_free_files.append(free_file_head+"_"+str(i)+".txt")
		all_beta_fix.append(fixed_beta_head+"_"+str(i)+".txt")
		all_alpha_fix.append(fixed_alpha_head+"_"+str(i)+".txt")


	# Dictionary of all iteration outputs 
	Free_Param_dict = {"popt":[], "ll_opt":[], "theta_opt":[]}
	Fixed_Beta_dict = {"popt_fix":[], "ll_opt_fix":[], "theta_opt_fix":[]}
	Fixed_Alpha_dict = {"popt_fix":[], "ll_opt_fix":[], "theta_opt_fix":[]}

	############ 
	### Iterate over all free param files to populate dictionary
	############

	for file in all_free_files:

		with open(file, 'r') as f:
			for line in f:
				list_params = line.split("\t")
				popt = list_params[0]
				ll_opt = float(list_params[1])
				theta_opt = float(list_params[2])

		popt_list = []
		for param in popt.split(", "):
			popt_list.append(float(param))

		Free_Param_dict['popt'].append(popt_list)
		Free_Param_dict['ll_opt'].append(ll_opt)
		Free_Param_dict['theta_opt'].append(theta_opt)


	############ 
	### Iterate over all fixed beta files to populate dictionary
	############

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

	############ 
	### Iterate over all fixed alpha files to populate dictionary
	############

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


	

	###########
	### Extract the best-fitting parameters 
	###########

	best_i_free, best_ll_free = get_best_fit(Free_Param_dict["ll_opt"])

	best_i_fix, best_ll_fix = get_best_fit(Fixed_Beta_dict['ll_opt_fix'])

	best_i_fix_a, best_ll_fix_a = get_best_fit(Fixed_Alpha_dict['ll_opt_fix'])


	# Extract the best fit parameters, ll, and theta for free param model
	free_best_alpha = Free_Param_dict['popt'][best_i_free][0]
	free_best_beta = Free_Param_dict['popt'][best_i_free][1]
	free_best_T = Free_Param_dict['popt'][best_i_free][2]
	free_best_nu = Free_Param_dict['popt'][best_i_free][3]
	free_best_ll = Free_Param_dict['ll_opt'][best_i_free]
	free_best_theta = Free_Param_dict['theta_opt'][best_i_free]

	# Extract the best fit parameters, ll, and theta for fixed beta model
	fix_best_alpha = Fixed_Beta_dict['popt_fix'][best_i_fix][0]
	fix_best_beta = Fixed_Beta_dict['popt_fix'][best_i_fix][1]
	fix_best_T = Fixed_Beta_dict['popt_fix'][best_i_fix][2]
	fix_best_nu = Fixed_Beta_dict['popt_fix'][best_i_fix][3]
	fix_best_ll = Fixed_Beta_dict['ll_opt_fix'][best_i_fix]
	fix_best_theta = Fixed_Beta_dict['theta_opt_fix'][best_i_fix]

	# Extract best fit params, ll, and theta for fixed alpha model 
	fix_a_best_alpha = Fixed_Alpha_dict['popt_fix'][best_i_fix][0]
	fix_a_best_beta = Fixed_Alpha_dict['popt_fix'][best_i_fix][1]
	fix_a_best_T = Fixed_Alpha_dict['popt_fix'][best_i_fix][2]
	fix_a_best_nu = Fixed_Alpha_dict['popt_fix'][best_i_fix][3]
	fix_a_best_ll = Fixed_Alpha_dict['ll_opt_fix'][best_i_fix]
	fix_a_best_theta = Fixed_Alpha_dict['theta_opt_fix'][best_i_fix]

	# Populate output dictionary 
	free_param_out_dict['alpha'].append(free_best_alpha)
	free_param_out_dict['beta'].append(free_best_beta)
	free_param_out_dict['T'].append(free_best_T)
	free_param_out_dict['nu'].append(free_best_nu)
	free_param_out_dict['ll_opt'].append(free_best_ll)
	free_param_out_dict['theta_opt'].append(free_best_theta)

	fixed_beta_out_dict['alpha_fix'].append(fix_best_alpha)
	fixed_beta_out_dict['beta_fix'].append(fix_best_beta)
	fixed_beta_out_dict['T_fix'].append(fix_best_T)
	fixed_beta_out_dict['nu_fix'].append(fix_best_nu)
	fixed_beta_out_dict['ll_opt_fix'].append(fix_best_ll)
	fixed_beta_out_dict['theta_opt_fix'].append(fix_best_theta)

	fixed_alpha_out_dict['alpha_fix'].append(fix_a_best_alpha)
	fixed_alpha_out_dict['beta_fix'].append(fix_a_best_beta)
	fixed_alpha_out_dict['T_fix'].append(fix_a_best_T)
	fixed_alpha_out_dict['nu_fix'].append(fix_a_best_nu)
	fixed_alpha_out_dict['ll_opt_fix'].append(fix_a_best_ll)
	fixed_alpha_out_dict['theta_opt_fix'].append(fix_a_best_theta)

	


free_df = pd.DataFrame(free_param_out_dict)
fixed_df = pd.DataFrame(fixed_beta_out_dict)
fix_a_df = pd.DataFrame(fixed_alpha_out_dict)

free_df.to_csv(prefix+"_free_param_opts.csv", index=True)
fixed_df.to_csv(prefix+"_fixed_beta_opts.csv", index=True)
fix_a_df.to_csv(prefix+"_fixed_alpha_opts.csv", index=True)


