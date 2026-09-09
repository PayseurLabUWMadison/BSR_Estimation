import dadi 
import scipy
import sys
import numpy as np
import nlopt
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument("--auto_data", type=str, help="Simulated spectra for the autosome")
parser.add_argument("--x_data", type=str, help='Simulated spectra for the X')
parser.add_argument("--index", type=int, help='Index for replication on CHTC')
parser.add_argument('--prefix', type=str, help='outfile prefix. Include information about simulations and model')
parser.add_argument('--fold', type = int, help='folded (1) or unfolded (0) SFS')
parser.add_argument('--Aseq', type = int, help='Length of sequence used for calculation of Autosomal SFS')
parser.add_argument('--Xseq', type = int, help='Length of sequence used for calculation of X SFS')

args = parser.parse_args()

alpha = 2 #From Wilson-Sayres et al. 2011 
a_fix = 2

if args.Xseq == None:
    XseqRatio = 1
else:
    XseqRatio = args.Xseq / args.Aseq # number of accessible bases for X and A 


# Test run: python MR_HT_dadi_optimize.py --auto_data ../Obs_Data/Autosome_FSFS_MR.fs --x_data ../Obs_Data/X_FSFS_MR.fs --prefix Maple_Ridge --fold 1 --index 0 

def concatenate_spectra(fsA, fsX, folded=False):
    """
    Concatenate two spectra, fsA and fsX, into a single spectrum.
    """
    if fsA.ndim != fsX.ndim:
        raise ValueError("Spectra must have the same number of dimensions.")
    if fsA.shape[0] != fsX.shape[0]:
        raise ValueError("Spectra must have the same shape.")
    
    fsCombined = dadi.Spectrum(np.ma.concatenate([fsA, fsX], axis=0), pop_ids=fsA.pop_ids)
    if hasattr(fsA, 'extrap_x'):
        fsCombined.extrap_x = fsA.extrap_x
   # Make sure the output spectra has the correct metadata for being folded and unfolded 
    if folded: 
        fsCombined.folded = True
    else:
        fsCombined.folded = False
    
    return fsCombined



@ dadi.Numerics.make_extrap_func
def MR_three_epoch_AX(params, ns, pts): ### TODO: Make this a three epoch model, with 5 parameters 
    """
    Two-epoch model for autosomes and X chromosome.

    params: [alpha, beta, Tmain, numain, Tcurr, nucurr, XseqRatio]
    alpha: male to female mutation rate ratio (mu_m/mu_f).
    beta: breeding ratio (Nf/Nm).
    T: time of the epoch in units of 2 * NaT, where NaT is the total ancestral population size.
    nu: relative *total* population size.
    XseqRatio: Relative sequence length of X chromosome to autosome
    """
    alpha, beta, Tmain, numain, Tcurr, nucurr = params[:6]
    if len(params) == 7:
        XseqRatio = params[6]
    else:
        XseqRatio = 1.0

    ns = [_//2 for _ in ns]

    xx = dadi.Numerics.default_grid(pts)

    phiA = dadi.PhiManip.phi_1D(xx, beta=beta)
    phiX = dadi.PhiManip.phi_1D_X(xx, alpha=alpha, beta=beta)

    phiA = dadi.Integration.one_pop(phiA, xx, T=Tmain, nu=numain, beta=beta)
    phiX = dadi.Integration.one_pop_X(phiX, xx, T=Tmain, nu=numain, alpha=alpha, beta=beta)

    phiA = dadi.Integration.one_pop(phiA, xx, T=Tcurr, nu=nucurr, beta=beta)
    phiX = dadi.Integration.one_pop_X(phiX, xx, T=Tcurr, nu=nucurr, alpha=alpha, beta=beta)

    fsA = dadi.Spectrum.from_phi(phiA, ns, (xx,))

    fsX = XseqRatio * dadi.Spectrum.from_phi(phiX, ns, (xx,))

    if args.fold == 1:
        fsA_fold = fsA.fold()
        fsX_fold = fsX.fold()
        fsCombined = concatenate_spectra(fsA_fold, fsX_fold, folded=True)
    else:
        fsCombined = concatenate_spectra(fsA, fsX, folded=False)

    return fsCombined



##########
### Read in Data
##########

Autosome_data = dadi.Spectrum.from_file(args.auto_data)

X_data = dadi.Spectrum.from_file(args.x_data)

# project autosome data down to size of X chromosome
Auto_data_down = Autosome_data.project(X_data.sample_sizes)

data_fs = concatenate_spectra(Auto_data_down, X_data, folded=True)


ns = Auto_data_down.sample_sizes * 2

pts_l = [30,40,50]

demo_model_AX = dadi.Numerics.make_extrap_func(MR_three_epoch_AX)


###########
### Initial Parameters
###########
# define starting parameters 
# [alpha, beta, Tmain, numain, Tcurr, nucurr, XseqRatio] 
start_params = [2, 1, 10, 5, 0.02, 0.5, XseqRatio]


# define boundaries of optimization (advice from doc: if optimization is close to bounds, expand)
lower_bounds = [1, 0.02, 2, 1, 0.001, 0.001, XseqRatio]
upper_bounds = [3, 8, 20, 25, 1, 2, XseqRatio]


## Fixed Alpha Model
popt_fix_a = [0,0,0,0,0,0,0]
#ensure beta does not stray from set bounds 
while popt_fix_a[1] < lower_bounds[1] or popt_fix_a[1] > upper_bounds[1] or popt_fix_a[2] < lower_bounds[2] or popt_fix_a[2] > upper_bounds[2] or popt_fix_a[3] < lower_bounds[3] or popt_fix_a[3] > upper_bounds[3] or popt_fix_a[4] < lower_bounds[4] or popt_fix_a[4] > upper_bounds[4] or popt_fix_a[5] < lower_bounds[5] or popt_fix_a[5] > upper_bounds[5]:
    p0_fix_a = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)

    print("Running fixed alpha optimization using initial parameters: ", p0_fix_a)

    popt_fix_a = dadi.Inference.optimize(p0_fix_a, data_fs, demo_model_AX, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, maxiter=600, fixed_params = [a_fix, None, None, None, None, None, XseqRatio])

print("Optimized fixed alpha parameters: ", popt_fix_a)

model_fix_a = demo_model_AX(popt_fix_a, ns, pts_l)

ll_opt_fix_a = dadi.Inference.ll_multinom(model_fix_a, data_fs)

theta_opt_fix_a = dadi.Inference.optimal_sfs_scaling(model_fix_a, data_fs)

print("Finished optimizing fixed alpha model")




## Nested Fixed Beta model (fix beta at 1)


#############
### Fixed Beta Model
#############

popt_fix_b = [0,0,0,0,0]
while popt_fix_b[2] < lower_bounds[2] or popt_fix_b[2] > upper_bounds[2] or popt_fix_b[3] < lower_bounds[3] or popt_fix_b[3] > upper_bounds[3] or popt_fix_b[4] < lower_bounds[4] or popt_fix_b[4] > upper_bounds[4] or popt_fix_b[5] < lower_bounds[5] or popt_fix_b[5] > upper_bounds[5]:
    p0_fix_b = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)

    print("Running fixed beta optimization using initial parameters: ", p0_fix_b)

    popt_fix_b = dadi.Inference.optimize(p0_fix_b, data_fs, demo_model_AX, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, maxiter=600, fixed_params = [a_fix, 1, None, None, None, None, XseqRatio])

print("Optimized fixed beta parameters: ", popt_fix_b)

model_fix_b = demo_model_AX(popt_fix_b, ns, pts_l)

ll_opt_fix_b = dadi.Inference.ll_multinom(model_fix_b, data_fs)

theta_opt_fix_b = dadi.Inference.optimal_sfs_scaling(model_fix_b, data_fs)

print("Finished optimizing fixed beta model")


###########
### Output Optimal Parameters
###########

popt_fix_out_a = popt_fix_a.tolist()
popt_fix_out_b = popt_fix_b.tolist()

output_alpha_fixed = [popt_fix_out_a, str(ll_opt_fix_a), str(theta_opt_fix_a)]
output_beta_fixed = [popt_fix_out_b, str(ll_opt_fix_b), str(theta_opt_fix_b)]

alpha_fix_name = args.prefix+"_alpha_fix_opt_"+str(args.index)+".txt"
beta_fix_name = args.prefix+"_beta_fix_opt_"+str(args.index)+".txt"

with open(alpha_fix_name, 'w') as f: 
    for item in output_alpha_fixed:
        if type(item) == list:
            f.write(', '.join(map(str, item))+"\t")
        else:
            f.write(item+'\t')

with open(beta_fix_name, 'w') as f:
	for item in output_beta_fixed:
		if type(item) == list:
			f.write(', '.join(map(str, item))+"\t")
		else:
			f.write(item+'\t')




