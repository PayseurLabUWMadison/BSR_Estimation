import pandas as pd 
import dadi 

obs_dir = '/Users/willspurley/Maple_Ridge_Inference/Obs_Data/'
prefix = "MR_10kbgenes_20kindep"


## Function to Create FS object from list 
def create_SFS_file(sfs, spec_type, filename):

	# create first line of sfs file 
	if spec_type == 'folded':
		sfs_len = (len(sfs)*2) - 1
		fold = 1
	else:
		sfs_len = len(sfs)
		fold = 0
	sfs_header = (str(sfs_len)+" "+str(spec_type))
    

    # create second lines of sfs file 
	sfs_data = ""
	for i in range(sfs_len):		
		if (i > len(sfs)-1) & (fold == 1):
			sfs_data = sfs_data+str(0)+" "
		else:
			sfs_data = sfs_data+str(sfs[i])+" "
        
    # create thrid line of sfs file 
	mask_data = []
	for i in range(sfs_len):
		if i == 0:
			mask_data.append(1)
		elif i == sfs_len:
			continue
		elif (i > len(sfs)-1) & (fold == 1):
			mask_data.append(1)
		else:
			mask_data.append(0)
	mask_line = ""
	for i in range(sfs_len-1): 
		mask_line = mask_line+str(mask_data[i])
		mask_line = mask_line+" "


	# Write the file using provided file name 
	with open(filename, 'w') as f:
		f.write(sfs_header+"\n"+sfs_data+"\n"+mask_line)
	print("Wrote frequency spectrum file: ", filename)



# Run example: create_SFS_file(sfs_list, spec_type = 'folded', filename="folded_SFS.fs")


## read in autosome and X FSFS for Maple Ridge
MR_A_FSFS = pd.read_csv(obs_dir+prefix+"_Autosome_FSFS.csv")

MR_X_FSFS = pd.read_csv(obs_dir+prefix+"_X_FSFS.csv")


# convert to lists 
A_FSFS_list = MR_A_FSFS.iloc[0].to_list()


X_FSFS_list = MR_X_FSFS.iloc[0].to_list()



create_SFS_file(A_FSFS_list, spec_type="folded", filename = obs_dir+prefix+"_Autosome_FSFS.fs")
create_SFS_file(X_FSFS_list, spec_type='folded', filename=obs_dir+prefix+"_X_FSFS.fs")







