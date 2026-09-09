import msprime
import allel
import numpy as np 
import pandas as pd 
import os
import argparse
import copy 
from fractions import Fraction
from decimal import Decimal, getcontext
import statistics 


## Goal: Using msprime coalescent simulation framework, generate input data for ABC-RF, dadi, and Musharoff et al. approaches. 
## We simulate three demographic models: a constant size, a 95% population reduction, and a 20x expansion. Each model is simulated alongside a BSR of 1/3, 1, and 3, but is designed to allow for any BSR using the proportion of females in the population (pf). 
## We have tailored this script to be run in parallel, with each replicate of observed data output individually before combining csvs together  
## Define the main funtion of the script - simulate X and Autosome to output both observed datasets for ABCRF and dadi inference 

def main():
	auto = Autosome(ind=Aind, rec_rate=rec_rate, mut_rate_A=mut_rate, length=length, pf=pf, demography=demo, alpha=alpha)
	auto_sims = auto.simulate(Aloci=Aloci, NA=NA, NC=NC, t=t, pf=pf, prefix=args.prefix)

	X_chrom = X_chr(ind=Xind, rec_rate_f=rec_rate, mut_rate_A=mut_rate, length=length, pf=pf, demography=demo, alpha=alpha)
	X_sims = X_chrom.simulate(Xloci=Xloci, NA=NA, NC=NC, t=t, pf=pf, prefix=args.prefix)

	final_df = pd.concat([auto_sims, X_sims], axis=1)

	final_df.to_csv("ABC_obs_data_"+args.prefix+"_"+str(args.rep)+".csv")



# -----------------------------
### Parse Command Line Arguments for population parameters ### 
# -----------------------------

# Allow for command line arguments for all parameters (use these to call arguments in CHTC jobs)
parser = argparse.ArgumentParser()
parser.add_argument("--rec_rate", type=float, help="Recombination rate input into simulations. Default = 5e-9")
parser.add_argument("--mut_rate", type=float, help='Mutation rate input into simulations. Default = 5e-9')
parser.add_argument("--pop_size", type=int, help = "Effective population size of the simulated population. Ignored if demography is specified. Default = 10,000")
parser.add_argument("--alpha", type = float, help="Male mutation rate / female mutation rate. Default = 1")
parser.add_argument("--demo", type=int, help="*REQUIRED* Specify the demography object to use for simulations [1-12] (both inclusive), enter 0 for None")
parser.add_argument("--prefix", type=str, help="*REQUIRED* Prefix to be included in output file names.")
parser.add_argument("--rec_rate_X", type=float, help="Recombination rate for the X chromosome in females. If not specified, input rec_rate_X = rec_rate_A")
parser.add_argument("--pf", type=float, help="BSR, Specified as the proportion of females Nf/(Nm+Nf). Default 0.5")
parser.add_argument("--rep", type = int, help = '*REQUIRED* Replicate this run is for parallel run on CHTC')
parser.add_argument('--name', type=str, help="population name")
parser.add_argument('--Aloci', type=int, help='Number of autosomal independent 10kb loci to simulate. Default 1000')
parser.add_argument('--Xloci', type=int, help = 'Number of X independent 10kb loci to simulate. If not specified, use same number of loci as Autosome')


args = parser.parse_args()

### Demographic Models to Simulate (This framework suits a two-epoch model, but can be adapted for three-epoch by altering the output demographic parameters)

if args.demo == 1:
	## create number 1: 95% reduction, 1000 gens ago
	pop_red = msprime.Demography()
	pop_red.add_population(name="A", initial_size=500)
	pop_red.add_population_parameters_change(time=1000, initial_size=10000)
	demo = pop_red
	NA = 10000
	NC = 500
	t = 1000

elif args.demo == 2:
	## create number 2: 95% expansion, 1000 generations ago
	pop_exp = msprime.Demography()
	pop_exp.add_population(name='A', initial_size = 200000)
	pop_exp.add_population_parameters_change(time=1000, initial_size = 10000)
	demo = pop_exp
	NA = 10000
	NC = 200000
	t = 1000

elif args.demo == 3:
	const = msprime.Demography()
	const.add_population(name='A', initial_size = 10000)
	demo = const
	NA = 10000
	NC = 10000
	t = 1000
	

# Define the prefix used for output 
prefix = args.prefix
if prefix == None:
	raise ValueError('Please specify an output file prefix using --prefix')

## Constraint for parallelization
rep = str(args.rep) 
if rep == None:
	raise ValueError("Please specify rep using --rep")

# Fixed Variables [3] - 
Aind = 10 # number of autosomal diploid 
Xind = 10
length = 10000

## Number of autosomal loci to explore 
Aloci = args.Aloci
if Aloci == None:
	Aloci = 1000

## Number of X-linked loci to explore. If not specified, assume equal numbers of X and autosomal loci 
Xloci = args.Xloci
if Xloci == None:
	Xloci = Aloci 


# Recombination rate to use, if not provided in command line set to 5e-8
if args.rec_rate != None:
	rec_rate = args.rec_rate
else:
	rec_rate = 5e-8

# If a unique value of recombination for the X is provided, use it. If not, set it equal to Autosome (Before adjusting based on BSR, see below)
if args.rec_rate_X == None:
	rec_rate_X = rec_rate
else:
	rec_rate_X = args.rec_rate_X

# Mutation rate to use, if not provided in command line set to 5e-8
if args.mut_rate != None:
	mut_rate = args.mut_rate
else:
	mut_rate = 5e-8

# Male:Female mutation rate
if args.alpha == None:
	alpha = 1
else:
	alpha = args.alpha

if args.name == None:
	name = prefix
else:
	name = args.name

# If demographic hisotry not specified, use constant size value. If that not specified, default to 10000
if args.demo == 0 or args.demo == None:
	if args.pop_size != None:
		size = args.pop_size
	else:
		size = 10000
else: size = None

# BSR to use for simulations, default is 0.5
if args.pf == None:
	pf = 0.5
else:
	pf = args.pf

# -----------------------------
### Functions to confer agreement between ABCRF and dadi parameters
# -----------------------------




# convert the proportion of females (Nf/(Nm+Nf)) to the BSR we interpret as (beta = Nf/Nm)
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



## Function to Create FS object from pandas dataframe - used to input into dadi optimization
def create_SFS_file(sfs, spec_type, filename):
    
    # create first line of sfs file 
    sfs_len = len(sfs)
    sfs_header = (str(sfs_len)+" "+str(spec_type))
    
    # create second lines of sfs file 
    sfs_data = ""
    for i in range(sfs_len):
        sfs_data = sfs_data+str(sfs[i])
        if i == sfs_len:
            continue
        sfs_data = sfs_data+" "
        
    # create thrid line of sfs file 
    mask_data = []
    for i in range(len(sfs)):
        if i == 0 or i == len(sfs)-1:
            mask_data.append(1)
        else:
            mask_data.append(0)
    mask_line = ""
    for i in range(sfs_len): 
        mask_line = mask_line+str(mask_data[i])
        if i == sfs_len:
            continue
        mask_line = mask_line+" "
    
    # Write the file using provided file name 
    with open(filename, 'w') as f:
        f.write(sfs_header+"\n"+sfs_data+"\n"+mask_line)
    print("Wrote frequency spectrum file: ", filename)


# -----------------------------
### Class for simulating Autosomes based on demographic history and BSR 
# -----------------------------


class Autosome:
	"""Parameters:
	 ___________________
	 ind: number of individuals to sample
	 rec_rate: sex-averaged per-site recombination rate 
	 mut_rate_A: autosomal per-site recombination rate 
	 ploidy: ploidy of sample (defaults to 2)
	 length: sequence length simulated segment
	 pf: proportion of breeding individuals that are female. 
	 pop_size: effective population size of simulated population (default 1, will be modified based on pf)
	 alpha: ratio of male per-site mutation rate to female per-site mutation rate (defaults to 1)
	 pf_out: proportion of females expressed as a string, replacing the period with an underscore (eg. "0_5" instead of 0.5)
	 ___________________
	 If running simualtions for populations with demograhpy object, the pop_size command will be ignored

	"""
	def __init__(self, ind=10, rec_rate=None, mut_rate_A=None, length=None, pf=0.5, pop_size=1, alpha = 1, demography=None):
		self.ind = ind
		self.ploidy = 2
		self.length = length
		self.pf = pf
		self.demography = copy.deepcopy(demography)
		aut_var = 4*pf*(1-pf) #adjust parameters according to the proportion of females in the population 


		# Modify effective population size based on breeding sex ratio
			## For demography object, adjust all epochs using the same BSR 
		if self.demography == None:
			self.pop_size = aut_var * pop_size
		else:
			for pop in self.demography.values():
				pop.initial_size *= aut_var
			for event in self.demography.events:
				if event.initial_size == None:
					continue
				event.initial_size *= aut_var

		self.rec_rate_A = rec_rate
		self.mut_rate_A = mut_rate_A
		self.alpha = alpha


		self.parent_path = os.getcwd()

 		# Report parameters for autosome after modifying for BSR 
		print("Autosome parameters (pf = "+str(self.pf)+"):")
		print("Mutation rate: ", self.mut_rate_A)
		print("Recombination rate: ", self.rec_rate_A)
		if self.demography == None:
			print("Effective Population Size: ", self.pop_size)
		else:
			print(self.demography)


	def simulate(self, Aloci, NA, NC, t, pf, prefix="Autosome", name="pop1"):

		
		self.prefix = prefix
		self.Aloci = Aloci
		self.pf = pf
		self.NA = NA
		self.NC = NC
		self.t = t 
		self.name = name
		num_ind = self.ind
		nchr = num_ind*self.ploidy
		seq_length = self.length

		## Include dadi terms in observed data 
		NC_NA = NC / NA
		t_NA = t / (2*NA)
		beta = pf_to_beta(pf)


		abc_out_dict = {"name":[],"NA":[], "NC":[], "t":[], "pf":[], "beta":[], "NC_NA_ratio":[], "t_NA_ratio":[]}
		abc_out_dict['name'].append(self.name)
		abc_out_dict["NA"].append(self.NA)
		abc_out_dict["NC"].append(self.NC)
		abc_out_dict["t"].append(self.t)
		abc_out_dict["pf"].append(self.pf)
		abc_out_dict['beta'].append(beta)
		abc_out_dict['NC_NA_ratio'].append(NC_NA)
		abc_out_dict['t_NA_ratio'].append(t_NA)




		abc_loci_dict_A = {"pi":[], "theta":[], "r2":[], "r2_var":[], "K":[], "M":[], "TD":[]} 

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

		Prop_SFS_dict_A = {}
		Prop_FSFS_dict_A = {}

		for key in SFS_dict_A.keys():
			Prop_SFS_dict_A["A_Prop_"+str(key)] = []

		for key in FSFS_dict_A.keys():
			Prop_FSFS_dict_A["A_Prop_"+str(key)] = []

		FSFS_list = []
		SFS_list = []


		if self.demography == None:
			#Simulate genealogy using the specified parameters in msprime. Replicates represent independent loci
			self.anc = msprime.sim_ancestry(samples = self.ind, sequence_length = self.length, recombination_rate = self.rec_rate_A, 
			ploidy = self.ploidy,
			num_replicates = self.Aloci, population_size = self.pop_size)
		else:
			self.anc = msprime.sim_ancestry(samples = self.ind, sequence_length = self.length, recombination_rate = self.rec_rate_A, 
			ploidy = self.ploidy, 
			num_replicates = self.Aloci, demography = self.demography)


		# _______________________________________________________________________________________________________________________________ #

		counter = 0

		for ts in self.anc: # loop for every locus in the chromosome
			# Add mutations to the genalogy of the ith locus
			mut_ts = msprime.sim_mutations(ts, rate = self.mut_rate_A) 

			vcf=mut_ts.as_vcf() #export as vcf to make readable file for scikit allel (mimic real data)

			with open(self.parent_path+"/temp_vcf"+str(counter)+".vcf", "w") as f:
				f.write(vcf)

			data = allel.read_vcf(self.parent_path+"/temp_vcf"+str(counter)+".vcf")

			os.remove(self.parent_path+"/temp_vcf"+str(counter)+".vcf")

			counter += 1

			#if no SNPs, reflect that with stats that can be populated and include NAs for the others
			if data == None:
				abc_loci_dict_A['pi'].append(0)
				abc_loci_dict_A['theta'].append(0)
				abc_loci_dict_A['r2'].append(np.nan)
				abc_loci_dict_A['r2_var'].append(np.nan)
				abc_loci_dict_A['K'].append(1)
				abc_loci_dict_A['M'].append(1)
				abc_loci_dict_A['TD'].append(np.nan)
				continue

			POS = data['variants/POS']
			gt = allel.GenotypeArray(data['calldata/GT'])
			ac = gt.count_alleles()

			if len(POS) == 0: 
				abc_loci_dict_A['pi'].append(0)
				abc_loci_dict_A['theta'].append(0)
				abc_loci_dict_A['r2'].append(np.nan)
				abc_loci_dict_A['r2_var'].append(np.nan)
				abc_loci_dict_A['K'].append(1)
				abc_loci_dict_A['M'].append(1)
				abc_loci_dict_A['TD'].append(np.nan)
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
				abc_loci_dict_A['pi'].append(0)
				abc_loci_dict_A['theta'].append(0)
				abc_loci_dict_A['r2'].append(np.nan)
				abc_loci_dict_A['r2_var'].append(np.nan)
				abc_loci_dict_A['K'].append(1)
				abc_loci_dict_A['M'].append(1)
				abc_loci_dict_A['TD'].append(np.nan)
				continue


			h = gt_seg_bi.to_haplotypes()
			dac_seg_bi = ac_seg_bi[:,1]

			#Calculate pi 
			pi, windows, n_bases, counts = allel.windowed_diversity(pos = POS_seg_bi, ac = ac_seg_bi, size = length, start=1, stop=length)
			abc_loci_dict_A['pi'].append(round(float(pi[0]), 8))

			#Calculate theta 
			theta, windows, n_bases, counts = allel.windowed_watterson_theta(pos=POS_seg_bi, ac=ac_seg_bi, size=length, start=1, stop=length)
			abc_loci_dict_A["theta"].append(round(float(theta[0]), 8))

			#Calculate Tajima's D
			D, windows, counts = allel.windowed_tajima_d(pos = POS_seg_bi, ac = ac_seg_bi, size= length, start=1, stop=length)
			abc_loci_dict_A["TD"].append(round(float(D[0]), 8))

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
				abc_loci_dict_A['r2'].append(np.nan)
				abc_loci_dict_A['r2_var'].append(np.nan)

			else:
				#Create r2 array for whole locus
				r_squared = allel.rogers_huff_r(gn_r2_final) ** 2
				#Calc mean r2
				meanr2 = float(np.mean(r_squared))
				abc_loci_dict_A['r2'].append(meanr2)
				#Calc variance in r2
				varr2 = float(np.var(r_squared))
				abc_loci_dict_A['r2_var'].append(varr2)

			#Calculate K
			unique_haps = h.distinct_counts() 
			K=len(unique_haps)
			abc_loci_dict_A['K'].append(K)

			#Calculate M
			sum_haps = sum(unique_haps)
			freq_hap=max(unique_haps)
			M = freq_hap/sum_haps
			abc_loci_dict_A['M'].append(M)

			#Calculate Folded SFS ((OUTPUT AS ALLELE COUNTS NOT PROPORTION))
			FSFS = allel.sfs_folded(ac_seg_bi, n=nchr)
			FSFS_list.append(FSFS)
			#Calculate Unfolded SFS (OUTPUT AS ALLELE COUNTS NOT PROPORTION)
			SFS = allel.sfs(dac_seg_bi, n = nchr)
						
			SFS_list.append(SFS)
                        
 
		# -----------------
		### Compile FSFS and SFS for dadi output ### 
		# -----------------

		dadi_SFS_bin_dict = {}
		for key in SFS_dict_A.keys():
			dadi_SFS_bin_dict[key] = []
		dadi_FSFS_bin_dict = {}
		for key in FSFS_dict_A.keys():
			dadi_FSFS_bin_dict[key] = []

		#Populate above dictionary with 
		for SFS in SFS_list:
			for i in range(len(SFS)):
				dadi_SFS_bin_dict[i].append(float(SFS[i]))
		for FSFS in FSFS_list:
			for i in range(len(FSFS)):
				dadi_FSFS_bin_dict[i].append(float(FSFS[i]))
				
		## Take the sum of each bin and output it to the correct bin in SFS/FSFS dict_A
		for key, bin_list in dadi_SFS_bin_dict.items():
			SFS_dict_A[key].append(sum(bin_list))

		for key, bin_list in dadi_FSFS_bin_dict.items():
			FSFS_dict_A[key].append(sum(bin_list))

		# Compile SFS and FSFS 
		dadi_sfs_df = pd.DataFrame(SFS_dict_A)
		dadi_fsfs_df = pd.DataFrame(FSFS_dict_A)

		fsfs = dadi_fsfs_df.values.tolist()
		fsfs_list = fsfs[0]

		create_SFS_file(sfs_list, spec_type = 'folded', filename="folded_dadi_auto_"+rep+"_"+prefix+".fs")


		# ------------------
		### Compile ABC Dataframe ###
		# ------------------

		# Convert FSFS to proportion of variants
		A_total_vars = 0
		for val in FSFS_dict_A.values():
			A_total_vars += val[0]

		for key, val in FSFS_dict_A.items():
			prop_bin = val[0] / A_total_vars
			Prop_FSFS_dict_A['A_Prop_'+str(key)] = prop_bin	

		abc_loci_A = pd.DataFrame(abc_loci_dict_A)
		mean_stats_A = abc_loci_A.mean()
		var_stats_A = abc_loci_A.var()

		for i in range(len(mean_stats_A)):
			abc_out_dict["A"+mean_stats_A.index[i]+"_m"] = [] 
			abc_out_dict["A"+var_stats_A.index[i]+"_v"] =  []
		for key in Prop_FSFS_dict_A.keys():
			abc_out_dict[key] = []

		for i in range(len(mean_stats_A)):
			abc_out_dict["A"+mean_stats_A.index[i]+"_m"].append(mean_stats_A.iloc[i])
			abc_out_dict["A"+var_stats_A.index[i]+"_v"].append(var_stats_A.iloc[i])

		#Append the proportion of variants in each frequency bin in the folded sfs
		for key, val in Prop_FSFS_dict_A.items():
			abc_out_dict[key].append(val)

		abc_out_df = pd.DataFrame(abc_out_dict)
		abc_out_df = abc_out_df.set_index('name')


		return abc_out_df



## /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/ ##


class X_chr:
	"""Parameters:
	________________
	ind: number of individuals to sample
	rec_rate_f: female-specific per-site recombination rate (will be scaled with changing pf)
	mut_rate_A: autosomal per-site mutation rate (will differ from autosomal value if alpha does not equal 1)
	length: sequence length of simulated segment
	pf: proportion of breeding individuals that are female
	pop_size: effective population size of simulated population – must equal autosomal entry for proper comparison
	alpha: ratio of male per-site mutation rate to female per-site mutation rate (defaults to 1)
	demography: msprime demography object for the desired population. Must include initial_sizes representing the Ne of that population size
	pf_out: proportion of females expressed as a string, replacing the period with an underscore (eg. "0_5" instead of 0.5)
	_________________
	If running simualtions for populations with demograhpy object, the pop_size command will be ignored

	"""
	def __init__(self, ind=10, rec_rate_f=None, mut_rate_A=None, length=None, pf=0.5, pop_size=1, alpha = 1, demography=None):
		self.ind = ind
		self.ploidy = 2
		self.length = length
		self.pf = pf
		self.alpha = alpha
		self.demography = copy.deepcopy(demography)
		self.parent_path = os.getcwd()


		x_var = (9 * pf * (1-pf))/(2 * (2-pf)) #adjust parameters based on the proportion of females in the population 
		alp_num = 2 * (2+alpha)
		alp_den = 3 * (1+alpha)
		alp_eq = alp_num/alp_den # allows for differences in mutation rate between sexes 

		if self.demography == None:	
			self.pop_size = x_var * pop_size
		else:
			for pop in self.demography.values():
				pop.initial_size *= x_var
			for event in self.demography.events:
				if event.initial_size == None:
					continue
				event.initial_size *= x_var

		self.rec_rate_X = (2*pf)/(1+pf) * rec_rate_f #scaling recombination rate down for all chr based on pf 
		self.mut_rate_X = mut_rate_A * alp_eq #scales mutation rate accounting for differing rates among males and females


		print("X chromosome parameters:")
		print("Mutation rate: ", self.mut_rate_X)
		print("Recombination rate: ", self.rec_rate_X)
		if self.demography == None:
			print("Effective Population Size: ", self.pop_size, "\n")
		else:
			print(self.demography, "\n")


	def simulate(self, Xloci, NA, NC, t, pf, prefix="X_chr", name="pop1"):

		
		self.prefix = prefix
		self.Xloci = Xloci
		self.pf = pf
		self.NA = NA
		self.NC = NC
		self.t = t 
		self.name = name
		num_ind = self.ind
		nchr = num_ind*self.ploidy
		seq_length = self.length

		## Include dadi terms in observed data 
		NC_NA = NC / NA
		t_NA = t / (2*NA)
		beta = pf_to_beta(pf)


		abc_out_dict = {"name":[]}
		abc_out_dict['name'].append(self.name)



		abc_loci_dict_X = {"pi":[], "theta":[], "r2":[], "r2_var":[], "K":[], "M":[], "TD":[]} 

		SFS_dict_X = {}
		FSFS_dict_X = {}

		num_ind = ind
		nchr = num_ind*2

		# Set the keys for SFS_dictA
		for i in range(0, nchr+1):
			SFS_dict_X[i] = []
		#Set the keys for FSFS_dictA
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


		if self.demography == None:
			#Simulate genealogy using the specified parameters in msprime. Replicates represent independent loci
			self.anc = msprime.sim_ancestry(samples = self.ind, sequence_length = self.length, recombination_rate = self.rec_rate_X, 
			ploidy = self.ploidy,
			num_replicates = self.Xloci, population_size = self.pop_size)
		else:
			self.anc = msprime.sim_ancestry(samples = self.ind, sequence_length = self.length, recombination_rate = self.rec_rate_X, 
			ploidy = self.ploidy, 
			num_replicates = self.Xloci, demography = self.demography)


		# _______________________________________________________________________________________________________________________________ #

		counter = 0

		for ts in self.anc: # loop for every locus in the chromosome
			# Add mutations to the genalogy of the ith locus
			mut_ts = msprime.sim_mutations(ts, rate = self.mut_rate_X) 

			vcf=mut_ts.as_vcf() #export as vcf to make readable file for scikit allel (mimic real data)

			with open(self.parent_path+"/temp_vcf"+str(counter)+".vcf", "w") as f:
				f.write(vcf)

			data = allel.read_vcf(self.parent_path+"/temp_vcf"+str(counter)+".vcf")

			os.remove(self.parent_path+"/temp_vcf"+str(counter)+".vcf")

			counter += 1

			#if no SNPs, reflect that with stats that can be populated and include NAs for the others
			if data == None:
				abc_loci_dict_X['pi'].append(0)
				abc_loci_dict_X['theta'].append(0)
				abc_loci_dict_X['r2'].append(np.nan)
				abc_loci_dict_X['r2_var'].append(np.nan)
				abc_loci_dict_X['K'].append(1)
				abc_loci_dict_X['M'].append(1)
				abc_loci_dict_X['TD'].append(np.nan)
				continue

			POS = data['variants/POS']
			gt = allel.GenotypeArray(data['calldata/GT'])
			ac = gt.count_alleles()

			if len(POS) == 0: 
				abc_loci_dict_X['pi'].append(0)
				abc_loci_dict_X['theta'].append(0)
				abc_loci_dict_X['r2'].append(np.nan)
				abc_loci_dict_X['r2_var'].append(np.nan)
				abc_loci_dict_X['K'].append(1)
				abc_loci_dict_X['M'].append(1)
				abc_loci_dict_X['TD'].append(np.nan)
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
				abc_loci_dict_X['pi'].append(0)
				abc_loci_dict_X['theta'].append(0)
				abc_loci_dict_X['r2'].append(np.nan)
				abc_loci_dict_X['r2_var'].append(np.nan)
				abc_loci_dict_X['K'].append(1)
				abc_loci_dict_X['M'].append(1)
				abc_loci_dict_X['TD'].append(np.nan)
				continue


			h = gt_seg_bi.to_haplotypes()
			dac_seg_bi = ac_seg_bi[:,1]

			#Calculate pi 
			pi, windows, n_bases, counts = allel.windowed_diversity(pos = POS_seg_bi, ac = ac_seg_bi, size = length, start=1, stop=length)
			abc_loci_dict_X['pi'].append(round(float(pi[0]), 8))

			#Calculate theta 
			theta, windows, n_bases, counts = allel.windowed_watterson_theta(pos=POS_seg_bi, ac=ac_seg_bi, size=length, start=1, stop=length)
			abc_loci_dict_X["theta"].append(round(float(theta[0]), 8))

			#Calculate Tajima's D
			D, windows, counts = allel.windowed_tajima_d(pos = POS_seg_bi, ac = ac_seg_bi, size= length, start=1, stop=length)
			abc_loci_dict_X["TD"].append(round(float(D[0]), 8))

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
				abc_loci_dict_X['r2'].append(np.nan)
				abc_loci_dict_X['r2_var'].append(np.nan)

			else:
				#Create r2 array for whole locus
				r_squared = allel.rogers_huff_r(gn_r2_final) ** 2
				#Calc mean r2
				meanr2 = float(np.mean(r_squared))
				abc_loci_dict_X['r2'].append(meanr2)
				#Calc variance in r2
				varr2 = float(np.var(r_squared))
				abc_loci_dict_X['r2_var'].append(varr2)

			#Calculate K
			unique_haps = h.distinct_counts() 
			K=len(unique_haps)
			abc_loci_dict_X['K'].append(K)

			#Calculate M
			sum_haps = sum(unique_haps)
			freq_hap=max(unique_haps)
			M = freq_hap/sum_haps
			abc_loci_dict_X['M'].append(M)

			#Calculate Folded SFS ((OUTPUT AS ALLELE COUNTS NOT PROPORTION))
			FSFS = allel.sfs_folded(ac_seg_bi, n=nchr)
			FSFS_list.append(FSFS)

			#Calculate Unfolded SFS (OUTPUT AS ALLELE COUNTS NOT PROPORTION)
			SFS = allel.sfs(dac_seg_bi, n = nchr)
						
			SFS_list.append(SFS)
                        

		# -----------------
		### Compile FSFS and SFS for dadi output ### 
		# -----------------

		dadi_SFS_bin_dict = {}
		for key in SFS_dict_X.keys():
			dadi_SFS_bin_dict[key] = []
		dadi_FSFS_bin_dict = {}
		for key in FSFS_dict_X.keys():
			dadi_FSFS_bin_dict[key] = []

		#Populate above dictionary with 
		for SFS in SFS_list:
			for i in range(len(SFS)):
				dadi_SFS_bin_dict[i].append(float(SFS[i]))
		for FSFS in FSFS_list:
			for i in range(len(FSFS)):
				dadi_FSFS_bin_dict[i].append(float(FSFS[i]))
				
		## Take the sum of each bin and output it to the correct bin in SFS/FSFS dict_X
		for key, bin_list in dadi_SFS_bin_dict.items():
			SFS_dict_X[key].append(sum(bin_list))

		for key, bin_list in dadi_FSFS_bin_dict.items():
			FSFS_dict_X[key].append(sum(bin_list))

		# Compile SFS and FSFS 
		dadi_sfs_df = pd.DataFrame(SFS_dict_X)
		dadi_fsfs_df = pd.DataFrame(FSFS_dict_X)

		sfs = dadi_sfs_df.values.tolist()
		sfs_list = sfs[0]

		create_SFS_file(sfs_list, spec_type = 'unfolded', filename="dadi_X_"+rep+"_"+prefix+".fs")


		# ------------------
		### Compile ABC Dataframe ###
		# ------------------

		# Convert FSFS to proportion of variants
		X_total_vars = 0
		for val in FSFS_dict_X.values():
			X_total_vars += val[0]

		for key, val in FSFS_dict_X.items():
			prop_bin = val[0] / X_total_vars
			Prop_FSFS_dict_X['X_Prop_'+str(key)] = prop_bin	

		abc_loci_X = pd.DataFrame(abc_loci_dict_X)
		mean_stats_X = abc_loci_X.mean()
		var_stats_X = abc_loci_X.var()

		for i in range(len(mean_stats_X)):
			abc_out_dict["X"+mean_stats_X.index[i]+"_m"] = [] 
			abc_out_dict["X"+var_stats_X.index[i]+"_v"] =  []
		for key in Prop_FSFS_dict_X.keys():
			abc_out_dict[key] = []

		for i in range(len(mean_stats_X)):
			abc_out_dict["X"+mean_stats_X.index[i]+"_m"].append(mean_stats_X.iloc[i])
			abc_out_dict["X"+var_stats_X.index[i]+"_v"].append(var_stats_X.iloc[i])

		#Append the proportion of variants in each frequency bin in the folded sfs
		for key, val in Prop_FSFS_dict_X.items():
			abc_out_dict[key].append(val)

		abc_out_df = pd.DataFrame(abc_out_dict)
		abc_out_df = abc_out_df.set_index('name')


		return abc_out_df





main()

