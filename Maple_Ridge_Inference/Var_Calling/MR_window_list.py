import pandas as pd

working_dir = '/Users/willspurley/X_Chr_Var_Calling/MR_windows'

A_data = pd.read_csv(working_dir+'/Maple_Ridge_Autosome_ind_windows.csv')
X_data = pd.read_csv(working_dir+'/Maple_Ridge_X_ind_windows.csv')


A_list = [str(a)+":"+str(b)+"-"+str(c) for a,b,c in zip(A_data['Chromosome'], A_data['window_start'], A_data['window_end'])]
X_list = [str(a)+":"+str(b)+"-"+str(c) for a,b,c in zip(X_data['Chromosome'], X_data['window_start'], X_data['window_end'])]



with open(working_dir+'/Maple_Ridge_Autosome_ind_windows.list', 'w') as file: 
	for locus in A_list:
		file.write(locus)
		file.write('\n')

with open(working_dir+'/Maple_Ridge_X_ind_windows.list', 'w') as file: 
	for locus in X_list:
		file.write(locus)
		file.write('\n')