## Goal: Extract 10kb windows that are independent for the autosomes and the X chromosome. Using VCF filtered for X chromosome quality, bed files masking genic regions, 
## and scikit-allel, read in variants. Once variants are read, use bed file to mark all sites in genome falling in masked regions as inaccessible. Then, compute sumamry stasitics 
## for 10kb windows along the genome. Set threshold for minimum inclusion of inaccessible sites for windows, remove windows close enough to still be LD. Output csv of retained windows with their summary statsitcs and metadata 

import allel
import pandas as pd 
import numpy as np 
import statistics

##steps used to create masked fasta 
#cat merged_beds/P_man_merged_inv_repeats_10k_genes.bed Maple_Ridge_X_regions.bed | bedtools sort | bedtools merge > Maple_Ridge_10kgenes_final_merged.bed 
# bedtools maskfasta -fi GCF_003704035.1_HU_Pman_2.1.3_genomic.fna -bed ../bed_files/Maple_Ridge_10kgenes_final_merged.bed -fo MR_10kgenes_final_masked_genomic.fna



working_dir = '/mnt/payseurgenomes/will/SR_Processing/'

#troubleshoot: ['NC_056031.1'] #
chroms  = ['NC_056008.1', 'NC_056009.1', 'NC_056010.1', 'NC_056011.1', 'NC_056012.1', 'NC_056013.1', 'NC_056014.1', 'NC_056015.1', 'NC_056016.1', 'NC_056017.1', 'NC_056018.1', 'NC_056019.1', 'NC_056020.1', 'NC_056021.1', 'NC_056022.1', 'NC_056023.1', 'NC_056024.1', 'NC_056025.1', 'NC_056026.1', 'NC_056027.1', 'NC_056028.1', 'NC_056029.1', 'NC_056030.1', 'NC_056031.1']

window_size = 10000
masked_fasta = working_dir+"reference_files/MR_5kgenes_final_masked_genomic.fna"
num_ind = 17

out_file_prefix = "MR_5kgenes"

# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

######
## Define Functions
######

def parse_mask(masked_fasta):

	contig = []
	start_idx = []
	stop_idx = []

	with open(masked_fasta) as m:  # open the fasta file 
		for idx, line in enumerate(m.readlines()):
			if ">" in line: 
				line = line.strip(">")
				line = line.strip("\n")

				contig.append(line)
				start_idx.append(idx)
				if idx > 0:
					stop_idx.append(idx - 1)
		stop_idx.append(idx)


	# make a dictionary with keys being each contig and values being the line numbers that contig falls on 
	contig_dict = dict(zip(contig, [list(tup) for tup in zip(start_idx, stop_idx)]))
	return contig_dict



def make_accessible(contig_dict, chr, masked_fasta):
	# start/stop indicies for chromosome in masked fasta 
	coords = contig_dict[chr] 

	# boolean array that will inform whether or not the location on the chromosome is accessible or not
	accessible = []


	with open(masked_fasta) as m: 
		lines = m.readlines()[(coords[0]+1):(coords[1]+1)] # lines of the nucleotides (without the header) that will be turned into T/F array 
		for line in lines: 
			# convert line into list of nucleotides
			line = line.strip("\n")
			elements = list(line)

			# convert list to boolean array, where any base that is not an 'N' is True 
			bool_elements = [(e != 'N') for e in elements]

			#append each base to the accessible array 
			[accessible.append(b) for b in bool_elements]

	return accessible


# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\



for chrom in chroms:

	if chrom == 'NC_056031.1':
		nchr = 24
	else:
		nchr = 34


	print("Starting processing for Chromosome ", chrom)

	VCF = working_dir+'VCFs/unrelated_5kbgenes_maple_ridge_'+chrom+'.vcf.gz' ##EDIT AS NEEDED 

	if chrom == 'NC_056031.1':
		indep_thresh = 500000
	else:
		indep_thresh = 250000
	
	#Read in VCF 
	data = allel.read_vcf(VCF)

	#create the accessiblity matrix for chromosome 
	contig_dict_chr = parse_mask(masked_fasta)
	accessible_mat = make_accessible(contig_dict_chr, chrom, masked_fasta)

	windows_start = 1
	windows_stop = len(accessible_mat)


	POS = data['variants/POS']
	# Genotype of each 
	gt = allel.GenotypeArray(data['calldata/GT'])
	
	## Mask the genotypes that fall within inacessible regions of the chromosome ## 
	masked_pos = []
	for position in POS:
		i = position-1 #POS is 1 based but want the 0-based index 
		masked_pos.append(accessible_mat[i])

	mask_pos = np.array(masked_pos)

	accessibility_mask = np.broadcast_to(mask_pos[:, np.newaxis], (len(POS),num_ind))

	gt.mask = ~accessibility_mask
	masked_gt = gt.fill_masked()

	ac = masked_gt.count_alleles()


	# filter out non-biallelic variants from the chromosome 
	biallelic = ac.is_biallelic_01()[:]
	ac_biallelic = ac.compress(biallelic,  axis=0)[:, :2]
	gt_biallelic = masked_gt.compress(biallelic, axis=0)
	POS_bial = POS[biallelic]

	# Ensure there are segregating variants in the vcf (should be the case, but just to be careful)
	seg_sites = ac.is_segregating()
	bi_seg_sites=seg_sites[biallelic]

	#remove the variants from the allele count data that are not biallelic or segregating  
	ac_seg_bi = ac_biallelic.compress(bi_seg_sites, axis=0)
	POS_seg_bi = POS_bial[bi_seg_sites]
	gt_seg_bi = gt_biallelic.compress(bi_seg_sites, axis=0)
	gn_seg_bi = gt_seg_bi.to_n_alt(fill=-1)

	## Diversity Caclulations ###

	# calculate nucleotide diversity for 10kb windows
	pi, pi_windows, pi_n_bases, pi_counts = allel.windowed_diversity(pos = POS_seg_bi, ac = ac_seg_bi, size = window_size, start=windows_start, stop=windows_stop, is_accessible=accessible_mat)

	## calculate other summary stastiscs using the same windows as pi calculation
	# Tajima's D
	TD, windows, counts = allel.windowed_tajima_d(pos = POS_seg_bi, ac = ac_seg_bi, windows=pi_windows)

	# Watterson's theta
	theta_W, windows, n_bases, counts = allel.windowed_watterson_theta(pos=POS_seg_bi, ac=ac_seg_bi, windows=pi_windows, is_accessible=accessible_mat)


	## R squared Calculations ###  

	#Remove singletons from gn and POS to get better estimation of r2
	f_alt = gn_seg_bi.sum(axis=1)
	gn_nosingle = np.array(f_alt != 1)
	gn_r2 = gn_seg_bi.compress(gn_nosingle, axis=0)
	POS_r2 = POS_seg_bi[gn_nosingle]


	#Further filter for r2 to remove any sites where all ind are heterozygous 
	to_remove = np.all(gn_r2==1, axis=1)
	gn_r2_final = gn_r2[~to_remove]
	POS_r2_final = POS_r2[~to_remove]


	## variance in r2 for a window
	def var_r2(gnw):
		r_squared = allel.rogers_huff_r(gnw) ** 2
		if len(gnw) == 0:
			return np.nan
		else:
			return np.nanvar(r_squared)

	## mean r2 in a window using same approach
	def mean_r2(gnw):
		r_squared = allel.rogers_huff_r(gnw) ** 2
		return np.nanmean(r_squared)

	r2_m, windows, counts = allel.windowed_statistic(POS_r2_final, gn_r2_final, mean_r2, windows=pi_windows, fill=np.nan)
	r2_v, windows, counts = allel.windowed_statistic(POS_r2_final, gn_r2_final, var_r2, windows=pi_windows, fill=np.nan)

	### Folded SFS Calculations ### 
    
	def folded_SFS(macw):
		num = nchr+1
		fsfs = allel.sfs_folded(macw, num)
		return fsfs

	fsfs, windows, counts = allel.windowed_statistic(POS_seg_bi, ac_seg_bi, statistic = folded_SFS, windows=pi_windows, fill=np.full(nchr//2+1, np.nan))

	# break windows list into start and end components 
	window_start = []
	window_end = []
	for window in pi_windows:
		window_start.append(window[0])
		window_end.append(window[-1])

	# convert window lists into pandas dataframe TODO: ADD IN 'mean_r2':r2_vals
	div_dict = {'pi':pi, 'Tajimas_D':TD, 'Watterson_theta':theta_W, 'mean_r2':r2_m, 'var_r2':r2_v, 'window_start':window_start, 'window_end':window_end, 'n_bases':pi_n_bases, 'var_counts':pi_counts}
	
	# Add FSFS counts to div_dict ## 

	# define fsfs keys in div dict 
	for i in range(fsfs.shape[1]):
		div_dict["Number_vars_bin_"+str(i)] = []

	# populate fsfs keys with the respective bins

	for spec in fsfs:
		for i in range(len(spec)):
			div_dict["Number_vars_bin_"+str(i)].append(spec[i])
      
    # Convert the dictionary to pandas dataframe 
	div_df = pd.DataFrame(div_dict)

	# Exclude windows where there are more than 90% inacessible regions
	div_df_no_mask = div_df[div_df['n_bases']>= 1000]

	print("Number of windows for "+chrom+ " with at least 1000 accessible bases: ", str(len(div_df_no_mask)))

	div_df_no_mask.to_csv(working_dir+'tables/'+out_file_prefix+"_"+chrom+'_div_windows_1000accessible.csv')


	# Extract windows that are independent of each other
	ind_window_start = []

	for i in range(len(div_df_no_mask)):
		# check distance from start of current window to end of previous 
		if i > 0 and abs(div_df_no_mask['window_start'].iloc[i] - div_df_no_mask['window_end'].iloc[i-1]) >= indep_thresh:
			ind_window_start.append(div_df_no_mask['window_start'].iloc[i])
		elif i < len(div_df_no_mask) - 1 and abs(div_df_no_mask['window_end'].iloc[i] - div_df_no_mask['window_start'].iloc[i+1]) >= indep_thresh:
			ind_window_start.append(div_df_no_mask['window_start'].iloc[i])


	div_df_no_mask_indep = div_df_no_mask[div_df_no_mask['window_start'].isin(ind_window_start)]

	print("Number of independent windows with at least 1000 accessible bases: ", len(div_df_no_mask_indep))

	div_df_no_mask_indep.to_csv(working_dir+'tables/'+out_file_prefix+"_"+chrom+'_div_windows_1000accessible_'+str(int(indep_thresh/1000))+'k_indep.csv')

	print("************************************", '\n')





