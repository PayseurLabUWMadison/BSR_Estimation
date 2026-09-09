
## ABC-RF Simbank Construction
# Parallelized simulation framework that can draw from a prior distriubiton of paramters

# Step zero: input priors 
# Step one: Specify the demographic model
# Step two: draw from prior distributions of all parameters
# Step three: simulate using these parameters
# Step four: summary statsitics of genomci loci 
# Step five: compile output file of mean and variance of sum. stats and parameters used for inference


import msprime
import allel
import numpy as np 
import pandas as pd 
import os
import argparse
import time
from fractions import Fraction
from decimal import Decimal, getcontext

start_time = time.time()

parser = argparse.ArgumentParser()
parser.add_argument("--rec_rate", type=float, help="Recombination rate input into simulations. Default = 5e-9")
parser.add_argument("--mut_rate", type=float, help="Mutation rate input into simulations. Defualt = 5e-9")
parser.add_argument("--Aloci", type=int, help='Number of Autosomal loci to simulate, default 1000')
parser.add_argument("--name", type=str, help="Name to include as the output file. Make it iterable if running in parallel") # some sort of unique identifier to add to row name of output file 
parser.add_argument("--start", type=int, help="Integer to start the count for loop. Useful for parallelization")
parser.add_argument("--file", type=int, help="Number to add to the end of file to help with merge step")
parser.add_argument("--alpha", type=float, help="ratio of male:female mutation rates. 1 indicates an equal mutation rate for males and females")
parser.add_argument("--Xloci", type=int, help="Number of X loci to simulate, if different from autosome. If not provided, assume X and autsome have same number of loci")
args = parser.parse_args()


if args.mut_rate == None:
	mut_rate = 5e-9
else:
	mut_rate = args.mut_rate

if args.rec_rate == None:
	rec_rate = 5e-9
else:
	rec_rate = args.rec_rate

if args.Aloci == None:
	Aloci = 1000
else: 
	Aloci = args.Aloci

if args.Xloci == None:
	Xloci = Aloci 
else:
	Xloci = args.Xloci


if args.start == None:
	start = 0
else:
	start = int(args.start)

if args.alpha == None:
	alpha = 1
else:
	alpha = args.alpha


if args.name == None:
	raise ValueError("Please enter a name for the output file")


# fixed parameters 
ind = 10
length = 10000
#assume alpha of 1 for now (will need some sort of estimate of it before simulating)


# # FUNCTIONS FOR ADJUSTING VALUES ## 

## Scale A parameters based on BSR 
def scale_A_N(N, pf):
	scalar = 4*pf*(1-pf)
	return N * scalar

## Scale X parameters based on BSR 
def scale_X_N(N, pf):
	scalar = (9 * pf * (1-pf))/(2 * (2-pf))
	return N * scalar

#Scale per-site, per gen recombination rate on X 
def scale_X_c(c, pf): 
	scalar = (2*pf)/(1+pf)
	return c * scalar
def scale_X_mu(mut_rate, alpha):
	alp_num = 2 * (2+alpha)
	alp_den = 3 * (1+alpha)
	alp_eq = alp_num/alp_den # allows for differences in mutation rate between sexes 
	return mut_rate * alp_eq

def decimal_to_fraction(decimal_num):
    """
    Converts a decimal number (float, string, or Decimal) to a simplified fraction.

    For floats, it is recommended to pass a string representation to avoid 
    floating-point precision issues.
    """
    
    # Convert to string if the input is a float to maintain precision
    if isinstance(decimal_num, float):
        decimal_str = str(decimal_num)
        # Use the string constructor for accuracy
        return Fraction(decimal_str)

    # Directly use the constructor for string or Decimal inputs
    return Fraction(decimal_num)

def pf_to_beta(pf):
    """
    Convert the input proportion of females (float) into the ratio of females to males
    """
    pf_frac = decimal_to_fraction(pf)
    Nf = pf_frac.numerator
    Nm = pf_frac.denominator - pf_frac.numerator
    
    beta = Nf / Nm
    
    if beta < 1:
        return round(beta, 9) 
    else:
        return(round(beta, 3))

out_dict = {}


## Simulate 100 replicates at a time, allowing for faster scaling to get a reference table of 1 million replicates 
for rep in range(start, start+100):

	name = args.name + str(rep)


	## Create prior distributions from uniform distributions and sample from them for each parameter ## 

	# Ancestral Ne
	NA_val = int(np.random.uniform(5000, 25000, 1)[0])

	# Breeding sex ratio
	pf_val = round(np.random.uniform(0.1, 0.9, 1)[0], 3)

	beta_val = pf_to_beta(pf_val)
	# Current Ne
	NC_val = int(np.random.uniform(50, 5000, 1)[0])

	# ratio of current size to ancestral size 
	NC_NA_ratio = NC_val / NA_val ## Formatted to match the dadi notation

	# Time of demographic shift 
	t_val = int(np.random.uniform(10, 20000, 1)[0])
	t_NA_ratio = t_val / (2*NA_val) ## Formatted to match the dadi notation
	## Adjust Ns for the X and Autosome ## 

	# print("Ancestral Ne: "+ str(NA_val))
	# print("Current Ne: "+ str(NC_val))
	# print("pf: "+ str(pf_val))
	# print("time of change: "+str(t_val)+"\n")

	#Ancestral Ne adjustment
	NA_A = scale_A_N(NA_val, pf_val)
	NA_X = scale_X_N(NA_val, pf_val)

	#Current Ne adjustment
	NC_A = scale_A_N(NC_val, pf_val)
	NC_X = scale_X_N(NC_val, pf_val)


	#scale recombination for the X 
	rec_rate_X = scale_X_c(rec_rate, pf_val)

	#scale mutation rate for X according to provided alpha 
	mut_rate_X = scale_X_mu(mut_rate, alpha)

	## Create the demography objects using adjusted parameters ## 
	Demo_A = msprime.Demography()
	Demo_A.add_population(name="A", initial_size = NC_A)
	Demo_A.add_population_parameters_change(time=t_val, initial_size=NA_A)

	Demo_X = msprime.Demography()
	Demo_X.add_population(name = "X", initial_size = NC_X)
	Demo_X.add_population_parameters_change(time=t_val, initial_size=NA_X)



	## Run Simulations following same appraoch as abc_dadi_combined_obs_sims.py ## 

	## Autosome first 
	anc_A = msprime.sim_ancestry(samples = ind, sequence_length = length, recombination_rate = rec_rate, ploidy = 2, num_replicates=Aloci, demography = Demo_A)

	# make a dictionary to store the summary statistics and site frequency spectra
	loci_dict_A = {"pi":[], "theta":[], "r2":[], "r2_var":[], "K":[], "M":[], "TD":[]} 

	SFS_dict_A = {}
	FSFS_dict_A = {}

	num_ind = ind
	nchr = num_ind*2

	# Set the keys for SFS_dictA
	for i in range(0, nchr+1):
		SFS_dict_A[i] = []
	#Set the keys for FSFS_dictA
	for i in range(0, num_ind+1):
		FSFS_dict_A[i] = []

	## dictionary for the SFS and FSFS 
	Prop_SFS_dict_A = {}
	Prop_FSFS_dict_A = {}

	## iterate over each bin of the SFS to create a key for each bin 
	for key in SFS_dict_A.keys():
		Prop_SFS_dict_A["A_Prop_"+str(key)] = []

	for key in FSFS_dict_A.keys():
		Prop_FSFS_dict_A["A_Prop_"+str(key)] = []

	FSFS_list = []
	SFS_list = []
	for ts in anc_A: ## Main simulation loop 
		#add mutations to the genealogy of loci i 
		mts = msprime.sim_mutations(ts, mut_rate)

		# convert to vcf to allow for summary statsitic calculation
		vcf = mts.as_vcf()

		with open("temp_vcf.vcf", "w") as f:
			f.write(vcf)
		data = allel.read_vcf("temp_vcf.vcf")

		os.remove("temp_vcf.vcf") # above chunk just allows allel to read the file into memory (as_vcf command not sufficient for some unknown reason)


		#if no SNPs, reflect that with stats that can be populated and include NAs for the others
		if data == None:
			loci_dict_A['pi'].append(0)
			loci_dict_A['theta'].append(0)
			loci_dict_A['r2'].append(np.nan)
			loci_dict_A['r2_var'].append(np.nan)
			loci_dict_A['K'].append(1)
			loci_dict_A['M'].append(1)
			loci_dict_A['TD'].append(np.nan)
			continue

		## Create data structures for scikit allel calculation of patterns of variation
		POS = data['variants/POS']
		gt = allel.GenotypeArray(data['calldata/GT'])
		ac = gt.count_alleles()

		if len(POS) == 0: 
			loci_dict_A['pi'].append(0)
			loci_dict_A['theta'].append(0)
			loci_dict_A['r2'].append(np.nan)
			loci_dict_A['r2_var'].append(np.nan)
			loci_dict_A['K'].append(1)
			loci_dict_A['M'].append(1)
			loci_dict_A['TD'].append(np.nan)
			continue

		## remove any abberant variants with more than two alleles 
		biallelic_check = ac.is_biallelic_01()[:]
		ac_biallelic = ac.compress(biallelic_check,  axis=0)[:, :2]
		gt_biallelic = gt.compress(biallelic_check, axis=0)
		POS_bial = POS[biallelic_check]


		#Remove non-segregating sites from biallelic SNPs 
		seg_sites = ac.is_segregating()
		bi_seg_sites=seg_sites[biallelic_check]
		ac_seg_bi = ac_biallelic.compress(bi_seg_sites, axis=0)
		gt_seg_bi = gt_biallelic.compress(bi_seg_sites, axis=0)
		POS_seg_bi = POS_bial[bi_seg_sites]
		gn_seg_bi = gt_seg_bi.to_n_alt()

		# if after filtering out non-seg sites there are no sites, stop here
		if len(POS_seg_bi) == 0: 
			loci_dict_A['pi'].append(0)
			loci_dict_A['theta'].append(0)
			loci_dict_A['r2'].append(np.nan)
			loci_dict_A['r2_var'].append(np.nan)
			loci_dict_A['K'].append(1)
			loci_dict_A['M'].append(1)
			loci_dict_A['TD'].append(np.nan)
			continue

		## needed for haplotype statistic calculation
		h = gt_seg_bi.to_haplotypes()
		dac_seg_bi = ac_seg_bi[:,1]

		#Calculate pi 
		pi, windows, n_bases, counts = allel.windowed_diversity(pos = POS_seg_bi, ac = ac_seg_bi, size = length, start=1, stop=length)
		loci_dict_A['pi'].append(round(float(pi[0]), 8))

		#Calculate theta 
		theta, windows, n_bases, counts = allel.windowed_watterson_theta(pos=POS_seg_bi, ac=ac_seg_bi, size=length, start=1, stop=length)
		loci_dict_A["theta"].append(round(float(theta[0]), 8))

		#Calculate Tajima's D
		D, windows, counts = allel.windowed_tajima_d(pos = POS_seg_bi, ac = ac_seg_bi, size= length, start=1, stop=length)
		loci_dict_A["TD"].append(round(float(D[0]), 8))

		#r squared calc
		#Remove singletons from gn and POS to get better estimation of r2
		f_alt = gn_seg_bi.sum(axis=1)
		gn_nosingle = np.array(f_alt != 1)
		gn_r2 = gn_seg_bi.compress(gn_nosingle, axis=0)
		POS_r2 = POS_seg_bi[gn_nosingle]
		#Further filter for r2 to remove any sites where all ind are heterozygous 
		to_remove = np.all(gn_r2==1, axis=1)
		gn_r2_final = gn_r2[~to_remove]
		POS_r2_final = POS_r2[~to_remove]

		#The r2 output requires a comparison between two SNPs, so if only one SNP remains after filtering singletons, correlations can't be found
		if len(POS_r2_final) <= 1:
			loci_dict_A['r2'].append(np.nan)
			loci_dict_A['r2_var'].append(np.nan)

		else:
			#Create r2 array for whole locus
			r_squared = allel.rogers_huff_r(gn_r2_final) ** 2
			#Calc mean r2
			meanr2 = float(np.mean(r_squared))
			loci_dict_A['r2'].append(meanr2)
			#Calc variance in r2
			varr2 = float(np.var(r_squared))
			loci_dict_A['r2_var'].append(varr2)

		#Calculate K
		unique_haps = h.distinct_counts() 
		K=len(unique_haps)
		loci_dict_A['K'].append(K)

		#Calculate M
		sum_haps = sum(unique_haps)
		freq_hap=max(unique_haps)
		M = freq_hap/sum_haps
		loci_dict_A['M'].append(M)

		#Calculate Folded SFS ((OUTPUT AS ALLELE COUNTS NOT PROPORTION))
		FSFS = allel.sfs_folded(ac_seg_bi, n=nchr)			
		FSFS_list.append(FSFS)
						
		#Calculate Unfolded SFS (OUTPUT AS ALLELE COUNTS NOT PROPORTION)
		# SFS = allel.sfs(dac_seg_bi, n = nchr)	
		# SFS_list.append(SFS)

	# Create intermediate dictionaries to collect each bin in its own list instead of each SFS in a list
	# ex {"1":[ALL SINGLETONS], "2":[ALL DOUBLETONS]} etc 
	SFS_bin_dict = {}
	for key in SFS_dict_A.keys():
		SFS_bin_dict[key] = []
	FSFS_bin_dict = {}
	for key in FSFS_dict_A.keys():
		FSFS_bin_dict[key] = []

	#Populate above dictionary with 
	for SFS in SFS_list:
		for i in range(len(SFS)):
			SFS_bin_dict[i].append(float(SFS[i]))
	for FSFS in FSFS_list:
		for i in range(len(FSFS)):
			FSFS_bin_dict[i].append(float(FSFS[i]))
			
	## Take the mean of each bin and output it to the correct bin in SFS/FSFS dict_A
	for key, bin_list in SFS_bin_dict.items():
		SFS_dict_A[key].append(sum(bin_list))
	for key, bin_list in FSFS_bin_dict.items():
		FSFS_dict_A[key].append(sum(bin_list))



	## Now for the X chromosome
	

	anc_X = msprime.sim_ancestry(samples = ind, sequence_length = length, recombination_rate = rec_rate_X, ploidy = 2, num_replicates=Xloci, demography = Demo_X)


	# dictionary for collecting summary statsitics of each locus 
	loci_dict_X = {"pi":[], "theta":[], "r2":[], "r2_var":[], "K":[], "M":[], "TD":[]}
	SFS_dict_X = {}
	FSFS_dict_X = {}

	# Set the keys for SFS_dictX
	for i in range(0, nchr+1):
		SFS_dict_X[i] = []

	#Set the keys for FSFS_dictX
	for i in range(0, num_ind+1):
		FSFS_dict_X[i] = []

	Prop_SFS_dict_X = {}
	Prop_FSFS_dict_X = {}

	for key in SFS_dict_X.keys():
		Prop_SFS_dict_X["X_Prop_"+str(key)] = []

	for key in FSFS_dict_X.keys():
		Prop_FSFS_dict_X["X_Prop_"+str(key)] = []


	FSFS_list = []
	SFS_list = []

	for ts in anc_X: 
		#add mutations to the genealogy of loci i 
		mts = msprime.sim_mutations(ts, mut_rate_X)

		# convert to vcf to allow for summary statsitic calculation
		vcf = mts.as_vcf()

		with open("temp_vcf.vcf", "w") as f:
			f.write(vcf)
		data = allel.read_vcf("temp_vcf.vcf")

		os.remove("temp_vcf.vcf") # above chunk just allows allel to read the file into memory (as_vcf command not sufficient for some unknown reason)

		if data == None:
			loci_dict_X['pi'].append(0)
			loci_dict_X['theta'].append(0)
			loci_dict_X['r2'].append(np.nan)
			loci_dict_X['r2_var'].append(np.nan)
			loci_dict_X['K'].append(1)
			loci_dict_X['M'].append(1)
			loci_dict_X['TD'].append(np.nan)
			continue

		POS = data['variants/POS']
		gt = allel.GenotypeArray(data['calldata/GT'])
		ac = gt.count_alleles()

		if len(POS) == 0: 
			loci_dict_X['pi'].append(0)
			loci_dict_X['theta'].append(0)
			loci_dict_X['r2'].append(np.nan)
			loci_dict_X['r2_var'].append(np.nan)
			loci_dict_X['K'].append(1)
			loci_dict_X['M'].append(1)
			loci_dict_X['TD'].append(np.nan)
			continue

		biallelic_check = ac.is_biallelic_01()[:]
		ac_biallelic = ac.compress(biallelic_check,  axis=0)[:, :2]
		gt_biallelic = gt.compress(biallelic_check, axis=0)
		POS_bial = POS[biallelic_check]


		#Remove non-segregating sites from biallelic SNPs 
		seg_sites = ac.is_segregating()
		bi_seg_sites=seg_sites[biallelic_check]
		ac_seg_bi = ac_biallelic.compress(bi_seg_sites, axis=0)
		gt_seg_bi = gt_biallelic.compress(bi_seg_sites, axis=0)
		POS_seg_bi = POS_bial[bi_seg_sites]
		gn_seg_bi = gt_seg_bi.to_n_alt()

		# if after filtering out non-seg sites there are no sites, stop here
		if len(POS_seg_bi) == 0: 
			loci_dict_X['pi'].append(0)
			loci_dict_X['theta'].append(0)
			loci_dict_X['r2'].append(np.nan)
			loci_dict_X['r2_var'].append(np.nan)
			loci_dict_X['K'].append(1)
			loci_dict_X['M'].append(1)
			loci_dict_X['TD'].append(np.nan)
			continue


		h = gt_seg_bi.to_haplotypes()
		dac_seg_bi = ac_seg_bi[:,1]

		#Calculate pi 
		pi, windows, n_bases, counts = allel.windowed_diversity(pos = POS_seg_bi, ac = ac_seg_bi, size = length, start=1, stop=length)
		loci_dict_X['pi'].append(round(float(pi[0]), 8))

		#Calculate theta 
		theta, windows, n_bases, counts = allel.windowed_watterson_theta(pos=POS_seg_bi, ac=ac_seg_bi, size=length, start=1, stop=length)
		loci_dict_X["theta"].append(round(float(theta[0]), 8))

		#Calculate Tajima's D
		D, windows, counts = allel.windowed_tajima_d(pos = POS_seg_bi, ac = ac_seg_bi, size= length, start=1, stop=length)
		loci_dict_X["TD"].append(round(float(D[0]), 8))

		#r squared calc
		#Remove singletons from gn and POS to get better estimation of r2
		f_alt = gn_seg_bi.sum(axis=1)
		gn_nosingle = np.array(f_alt != 1)
		gn_r2 = gn_seg_bi.compress(gn_nosingle, axis=0)
		POS_r2 = POS_seg_bi[gn_nosingle]
		#Further filter for r2 to remove any sites where all ind are heterozygous 
		to_remove = np.all(gn_r2==1, axis=1)
		gn_r2_final = gn_r2[~to_remove]
		POS_r2_final = POS_r2[~to_remove]

		#The r2 output requires a comparison between two SNPs, so if only one SNP remains after filtering singletons, correlations can't be found
		if len(POS_r2_final) <= 1:
			loci_dict_X['r2'].append(np.nan)
			loci_dict_X['r2_var'].append(np.nan)

		else:
			#Create r2 array for whole locus
			r_squared = allel.rogers_huff_r(gn_r2_final) ** 2
			#Calc mean r2
			meanr2 = float(np.mean(r_squared))
			loci_dict_X['r2'].append(meanr2)
			#Calc variance in r2
			varr2 = float(np.var(r_squared))
			loci_dict_X['r2_var'].append(varr2)

		#Calculate K
		unique_haps = h.distinct_counts() 
		K=len(unique_haps)
		loci_dict_X['K'].append(K)

		#Calculate M
		sum_haps = sum(unique_haps)
		freq_hap=max(unique_haps)
		M = freq_hap/sum_haps
		loci_dict_X['M'].append(M)

		#Calculate Folded SFS (OUTPUT AS ALLELE COUNTS NOT PROPORTION)
		FSFS = allel.sfs_folded(ac_seg_bi, n=num_ind*2)
		FSFS_list.append(FSFS)  

		# #Calculate Unfolded SFS (OUTPUT AS ALLELE COUNTS NOT PROPORTION) 
		# SFS = allel.sfs(dac_seg_bi, n = num_ind*2)
		# SFS_list.append(SFS)



	# SFS_bin_dict = {}
	for key in SFS_dict_X.keys():
		SFS_bin_dict[key] = []
	FSFS_bin_dict = {}
	for key in FSFS_dict_X.keys():
		FSFS_bin_dict[key] = []

	# #Populate above dictionary
	for SFS in SFS_list:
		for i in range(len(SFS)):
			SFS_bin_dict[i].append(float(SFS[i]))
	for FSFS in FSFS_list:
		for i in range(len(FSFS)):
			FSFS_bin_dict[i].append(float(FSFS[i]))
  
	## Take the sum of each bin and output it to the correct bin in SFS/FSFS dict_X
	for key, bin_list in SFS_bin_dict.items():
		SFS_dict_X[key].append(sum(bin_list))

	for key, bin_list in FSFS_bin_dict.items():
		FSFS_dict_X[key].append(sum(bin_list))



	## Calculate the proportion of variants in each bin of folded SFS for the autosome and X chromosome ##
	
	# Calculate the total number of variants in this replicate for the X and autosome 
	A_total_vars = 0
	for val in FSFS_dict_A.values():
		A_total_vars += val[0]
	X_total_vars = 0
	for val in FSFS_dict_X.values():
		X_total_vars += val[0]

	# Determine the proporiton of variants in each bin, outputting into Prop dictionary 
	for key, val in FSFS_dict_A.items():
		prop_bin = val[0] / A_total_vars
		Prop_FSFS_dict_A['A_Prop_'+str(key)] = prop_bin

	for key, val in FSFS_dict_X.items():
		prop_bin = val[0] / X_total_vars
		Prop_FSFS_dict_X['X_Prop_'+str(key)] = prop_bin



	#Create dataframes from dictionaries to make it easier to work with data 
	loci_A = pd.DataFrame(loci_dict_A)
	loci_X = pd.DataFrame(loci_dict_X)


	# calculate interlocus means and variances, excluding any NA's 
	mean_stats_A = loci_A.mean()
	mean_stats_X = loci_X.mean()
	var_stats_A = loci_A.var()
	var_stats_X = loci_X.var()

	# construct the output dictionary keys 
	if rep == start:
		for i in range(len(mean_stats_A)):
			out_dict["A"+mean_stats_A.index[i]+"_m"] = [] 
			out_dict["X"+mean_stats_X.index[i]+"_m"] = [] 
			out_dict["A"+var_stats_A.index[i]+"_v"] =  []
			out_dict['X'+var_stats_X.index[i]+"_v"] = []
		for key in Prop_FSFS_dict_A.keys():
			out_dict[key] = []
		for key in Prop_FSFS_dict_X.keys():
			out_dict[key] = []
		out_dict['NA'] = []
		out_dict['NC'] = []
		out_dict['pf'] = []
		out_dict['t'] = []
		out_dict['NC_NA_ratio'] = []
		out_dict['t_NA_ratio'] = []
		out_dict['beta'] = []
		out_dict['name'] = []
		print(out_dict)



	# Fill in reported values for the given replciate 
	for i in range(len(mean_stats_A)):
		out_dict["A"+mean_stats_A.index[i]+"_m"].append(mean_stats_A.iloc[i])
		out_dict["X"+mean_stats_X.index[i]+"_m"].append(mean_stats_X.iloc[i])
		out_dict["A"+var_stats_A.index[i]+"_v"].append(var_stats_A.iloc[i])
		out_dict['X'+var_stats_X.index[i]+"_v"].append(var_stats_X.iloc[i])

	#Append the proportion of variants in each frequency bin in the folded sfs
	for key, val in Prop_FSFS_dict_A.items():
		out_dict[key].append(val)

	for key, val in Prop_FSFS_dict_X.items():
		out_dict[key].append(val)

	# add the parameters used to simualte the data in order to keep track later 
	out_dict['NA'].append(NA_val)
	out_dict['NC'].append(NC_val)
	out_dict['pf'].append(pf_val)
	out_dict['t'].append(t_val)
	out_dict['NC_NA_ratio'].append(NC_NA_ratio)
	out_dict['t_NA_ratio'].append(t_NA_ratio)
	out_dict['beta'].append(beta_val)
	out_dict['name'].append(name)


out_df = pd.DataFrame(out_dict)

out_df = out_df.set_index('name')

#  Write a csv file of out_df

end_time = time.time()
runtime = end_time - start_time

print("Run time (s):", runtime)
out_df.to_csv(args.name+"_"+str(args.file)+"v5.csv")




