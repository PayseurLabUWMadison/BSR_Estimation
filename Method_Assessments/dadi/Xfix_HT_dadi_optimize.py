
import dadi 
import scipy
import sys
import numpy as np
import nlopt
import argparse 

## Goal: Optimize demographic parameters and the BSR simulatneously using concatenated autosomal and X-linked SFS 
# We compute three different models for this optimization process: 
#   Model 1 = Free Parameter model (alpha, BSR, and demography all freely inferred)
#   Model 2 = Fixed Alpha model (alpha fixed at correct simulated value, BSR and demography inferred)
#   Model 3 = Nested model (alpha fixed at correct value, BSR fixed at 1, demography is inferred)

# using FS files for A and X constructed by abc_dadi_combined_obs_sims.py, concatentate spectra and run the optimiation for each of these three models 

#python Xfix_HT_dadi_optimize.py --auto_data dadi_auto_53_pop_red_pf0_75_alpha1_Xloci100.fs --x_data dadi_X_53_pop_red_pf0_75_alpha1_Xloci100.fs --index 1 --prefix TESTF --XseqRatio 0.1 --fold 1

parser = argparse.ArgumentParser()
parser.add_argument("--auto_data", type=str, help="Simulated spectra for the autosome")
parser.add_argument("--x_data", type=str, help='Simulated spectra for the X')
parser.add_argument("--index", type=int, help='Index for replication on CHTC')
parser.add_argument('--prefix', type=str, help='outfile prefix. Include information about simulations and model')
parser.add_argument("--XseqRatio", type=float, help='Ratio of X loci to Autosomal loci')
parser.add_argument('--fold', type = int, help='folded (1) or unfolded (0) SFS')

args = parser.parse_args()

print("Fold arg: ", args.fold)

## Setting parameter thresholds based on simulated alpha and demographic scenario
alpha = 1
nu = 'red'

if alpha == 3:
    a_start = 3
    a_lower = 1
    a_upper = 5
    a_fix = 3
elif alpha == 1:
    a_start = 1
    a_lower = 0.5
    a_upper = 3
    a_fix = 1

if nu == 'red':
    nu_start = 0.05
    nu_lower = 0.001
    nu_upper = 1
elif nu == 'exp':
    nu_start = 20
    nu_lower = 1
    nu_upper = 100
elif nu == 'const':
    nu_start = 1
    nu_lower = 0.25
    nu_upper = 4

if args.XseqRatio == None:
    XseqRatio = 1
else:
    XseqRatio = args.XseqRatio


#########
### Define Useful Functions
#########

## From Ryan Gutenkunst, combined A and X fs objects into one 
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

## From Ryan Gutenkunst, function for inference of BSR and demogarphic params using concantenated spectrum 
@ dadi.Numerics.make_extrap_func
def two_epoch_AX(params, ns, pts):
    """
    Two-epoch model for autosomes and X chromosome.

    params: [alpha, beta, T, nu]
    alpha: male to female mutation rate ratio (mu_m/mu_f).
    beta: breeding ratio (Nf/Nm).
    T: time of the epoch in units of 2 * NaT, where NaT is the total ancestral population size.
    nu: relative *total* population size.
    """
    alpha, beta, T, nu = params[:4]
    if len(params) == 5:
        XseqRatio = params[4]
    else:
        XseqRatio = 1.0

    ns = [_//2 for _ in ns]

    xx = dadi.Numerics.default_grid(pts)

    phiA = dadi.PhiManip.phi_1D(xx, beta=beta)
    phiX = dadi.PhiManip.phi_1D_X(xx, alpha=alpha, beta=beta)

    phiA = dadi.Integration.one_pop(phiA, xx, T=T, nu=nu, beta=beta)
    phiX = dadi.Integration.one_pop_X(phiX, xx, T=T, nu=nu, alpha=alpha, beta=beta)

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

## Read in A and X fs seperately 
Autosome_data = dadi.Spectrum.from_file(args.auto_data)
X_data = dadi.Spectrum.from_file(args.x_data)


## Simulated SFS was unfolded, but can be folded using built in dadi functions
if args.fold == 1:
    folded_A = Autosome_data.fold() 
    folded_X = X_data.fold()
    data_fs = concatenate_spectra(folded_A, folded_X, folded=True)
else:
    data_fs = concatenate_spectra(Autosome_data, X_data)

print(data_fs)

###########
### Specify Model
###########

## Sample size needs to be 2x actual value for inference of concatenated SFS to work 
ns = [40]
## grid search can still operate using minimum value of original sample size + 10 
pts_l = [30, 40, 50]

demo_model_AX = dadi.Numerics.make_extrap_func(two_epoch_AX) # create model function


###########
### Initial Parameters
###########
# define starting parameters 
# [alpha, beta, T, nu, Xseqratio] 
start_params = [a_start, 1, 0.05, nu_start, XseqRatio]

# define boundaries of optimization (advice from dadi documentation: if optimization is close to bounds, expand)
lower_bounds = [a_lower, 0.02, 0.001, nu_lower, XseqRatio]
upper_bounds = [a_upper, 7, 2, nu_upper, XseqRatio]

###########
### Free Parameter Model
###########
# on rare occasions, the optimization process will return negative values for some parameters. To prevent this from crashing the process, rerun optimization until value falls within bounds 

popt = [0,0,0,0,0]
while popt[0] < lower_bounds[0] or popt[0] > upper_bounds[0] or popt[1] < lower_bounds[1] or popt[1] > upper_bounds[1]:
    # randomize starting point of optimization    
    p0 = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)

    print("Running free parameter optimization using initial parameters: ", p0)
    # optimize parameters based on simulated data_fs 
    popt = dadi.Inference.optimize(p0, data_fs, demo_model_AX, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, maxiter=600, verbose=0, fixed_params = [None, None, None, None, XseqRatio])

print("Optimized free model parameters: ", popt, '\n')

# fit model SFS based on best-fitting parameters from optimization 
model_fs = demo_model_AX(popt, ns, pts_l)

#calculate log-likelihood of model fit to data 
ll_opt = dadi.Inference.ll_multinom(model_fs, data_fs)

# return estimate of theta for interpetation of parameters
theta_opt = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)

print("Finished optimizing free paramter model", '\n')


############
### Fixed Alpha Model
############

popt_fix_a = [0,0,0,0,0]
#ensure beta does not stray from set bounds 
while popt_fix_a[1] < lower_bounds[1] or popt_fix_a[1] > upper_bounds[1] or popt_fix_a[2] < lower_bounds[2] or popt_fix_a[2] > upper_bounds[2]:
    p0_fix_a = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)

    print("Running fixed alpha optimization using initial parameters: ", p0_fix_a)

    popt_fix_a = dadi.Inference.optimize(p0_fix_a, data_fs, demo_model_AX, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, maxiter=600, fixed_params = [a_fix, None, None, None, XseqRatio])

print("Optimized fixed alpha parameters: ", popt_fix_a)

model_fix_a = demo_model_AX(popt_fix_a, ns, pts_l)

ll_opt_fix_a = dadi.Inference.ll_multinom(model_fix_a, data_fs)

theta_opt_fix_a = dadi.Inference.optimal_sfs_scaling(model_fix_a, data_fs)

print("Finished optimizing fixed alpha model")

#############
### Nested Model
#############

popt_fix_b = [0,0,0,0,0]
while popt_fix_b[0] < lower_bounds[0] or popt_fix_b[0] > upper_bounds[0] or popt_fix_b[2] < lower_bounds[2] or popt_fix_b[2] > upper_bounds[2] or popt_fix_b[3] < lower_bounds[3] or popt_fix_b[3] > upper_bounds[3]:
    p0_fix_b = dadi.Misc.perturb_params(start_params, fold=1, upper_bound=upper_bounds, lower_bound=lower_bounds)

    print("Running nested model optimization using initial parameters: ", p0_fix_b)

    popt_fix_b = dadi.Inference.optimize(p0_fix_b, data_fs, demo_model_AX, pts_l, lower_bound=lower_bounds, upper_bound=upper_bounds, maxiter=600, fixed_params = [a_fix, 1, None, None, XseqRatio])

print("Optimized nested model parameters: ", popt_fix_b)

model_fix_b = demo_model_AX(popt_fix_b, ns, pts_l)

ll_opt_fix_b = dadi.Inference.ll_multinom(model_fix_b, data_fs)

theta_opt_fix_b = dadi.Inference.optimal_sfs_scaling(model_fix_b, data_fs)

print("Finished optimizing nested model")


###########
### Output Optimal Parameters
###########

popt_out = popt.tolist()
popt_fix_out_a = popt_fix_a.tolist()
popt_fix_out_b = popt_fix_b.tolist()


output_all_free = [popt_out, str(ll_opt), str(theta_opt)]
output_alpha_fixed = [popt_fix_out_a, str(ll_opt_fix_a), str(theta_opt_fix_a)]
output_beta_fixed = [popt_fix_out_b, str(ll_opt_fix_b), str(theta_opt_fix_b)]


free_file_name = args.prefix+"_free_param_opt_"+str(args.index)+".txt"
alpha_fix_name = args.prefix+"_alpha_fix_opt_"+str(args.index)+".txt"
nested_name = args.prefix+"_nested_opt_"+str(args.index)+".txt"


with open(free_file_name, 'w') as f:
	for item in output_all_free:
		if type(item) == list:
			f.write(', '.join(map(str, item))+"\t")
		else:
			f.write(item+'\t')

with open(alpha_fix_name, 'w') as f:
    for item in output_alpha_fixed:
        if type(item) == list:
            f.write(', '.join(map(str, item))+"\t")
        else:
            f.write(item+'\t')


with open(nested_name, 'w') as f:
    for item in output_beta_fixed:
        if type(item) == list:
            f.write(', '.join(map(str, item))+"\t")
        else:
            f.write(item+'\t')











