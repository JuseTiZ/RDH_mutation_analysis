import nlopt
import dadi
import numpy
from dadi import Numerics, PhiManip, Integration, Demographics1D
from dadi.Spectrum_mod import Spectrum
import argparse
import pickle
import os
import fcntl
dadi.Integration.timescale_factor = 1e-6

def get_args():

    parser = argparse.ArgumentParser(
        description="Demographic inference with dadi")
    
    input_parser = parser.add_argument_group("Input")
    input_parser.add_argument("-v", "--vcf", type=str,
                        help="Input vcf for spectrum generation (synonymous vcf).")
    input_parser.add_argument("--popinfo", type=str,
                        help="Input population info file for spectrum generation.")
    input_parser.add_argument("--data_dict", type=str, default=None,
                        help="Input datadict for spectrum generation (synonymous).")
    input_parser.add_argument("--fs", type=str,
                        help="Input projected synonymous 1D spectrum (.fs).")
    
    param_parser = parser.add_argument_group("Parameters")
    param_parser.add_argument("-s", "--seed", type=int, default=42,
                        help="First seed for perturbing parameters. Default: 42")
    param_parser.add_argument("--starts", type=int, default=1,
                        help="Number of optimization starts with consecutive seeds. Use 100 for the full analysis. Default: 1")
    param_parser.add_argument("--model", choices=['bpg', 'two_epoch', 'constant'], default='bpg',
                        help="Demographic model to fit. Default: bpg")
    param_parser.add_argument("--popid", type=str,
                        help="Population ID for infer 1D demography.")
    param_parser.add_argument("--ns", type=int,
                        help="Sample size to project down.")
    param_parser.add_argument("--fold_fs", action='store_true',
                        help="Fold the frequency spectrum generated from variant data.")
    
    output_parser = parser.add_argument_group("Output")
    output_parser.add_argument("-o", "--out", type=str, required=True,
                        help="Output file name for demographic inference results.")
    
    args = parser.parse_args()
    if args.starts < 1:
        parser.error('--starts must be at least 1')
    if args.fs is None and args.ns is None:
        parser.error('--ns is required when --fs is not provided')
    return args


def main():

    args = get_args()

    fold_state = 'unfold' if not args.fold_fs else 'fold'
    polarized = False if args.fold_fs else True
    fs_file = args.out + f'.samplesize_{args.ns}.{fold_state}.fs'
    if args.fs is not None:
        print("Loading frequency spectrum from {}".format(args.fs))
        data_fs = dadi.Spectrum.from_file(args.fs)

    elif not os.path.exists(fs_file):
        if args.popid is None:
            raise ValueError('--popid is required to generate a frequency spectrum.')
        if args.data_dict is not None:
            print("Using data dictionary from {}".format(args.data_dict))
            syn_dd = pickle.load(open(args.data_dict, 'rb'))
        else:
            print("Generating data dictionary from VCF: {}".format(args.vcf))
            syn_dd = dadi.Misc.make_data_dict_vcf(args.vcf, args.popinfo)
        
        pop_ids, ns = [args.popid], [args.ns]
        data_fs = dadi.Spectrum.from_data_dict(syn_dd, pop_ids, ns, polarized=polarized)
        # Save the frequency spectrum
        tmp_file = fs_file + f'.{os.getpid()}.tmp'
        try:
            data_fs.to_file(tmp_file)
            os.replace(tmp_file, fs_file)
        finally:
            if os.path.exists(tmp_file):
                os.remove(tmp_file)
        print("Frequency spectrum saved to {}".format(fs_file))

    else:
        print("Loading frequency spectrum from {}".format(fs_file))
        data_fs = dadi.Spectrum.from_file(fs_file)

    ns = data_fs.sample_sizes
    fold_state = 'fold' if data_fs.folded else 'unfold'
    print("Sample sizes: {}".format(ns))
    pts_l = [max(ns)+120, max(ns)+130, max(ns)+140]

    def bottleneck_plateau_growth(params, ns, pts):

        nu1, nu2, nuc, T1, T2, Tc = params
        xx  = Numerics.default_grid(pts)
        phi = PhiManip.phi_1D(xx)      

        phi = Integration.one_pop(phi, xx, T1, nu1)
        phi = Integration.one_pop(phi, xx, T2, nu2)
        nu_func = lambda t: nu2 * numpy.exp(numpy.log(nuc/nu2) * (t / Tc))
        phi = Integration.one_pop(phi, xx, Tc, nu_func)

        fs = Spectrum.from_phi(phi, ns, (xx,))
        return fs

    # Define the demographic model and parameter bounds based on the selected model.
    if args.model == 'bpg':
        demo_model = bottleneck_plateau_growth
        param_names = ['nu1', 'nu2', 'nuc', 'T1', 'T2', 'Tc']
        params = [0.1, 1, 10, 0.01, 0.05, 0.05]
        lower_bounds = [1e-2, 1e-2, 1, 1e-3, 1e-3, 1e-3]
        upper_bounds = [5, 5, 200, 1, 1, 1]
        fixed_params = [None, None, None, 0.01, None, None]
    elif args.model == 'two_epoch':
        demo_model = Demographics1D.two_epoch
        param_names = ['nu', 'T']
        params = [2, 0.1]
        lower_bounds = [1.00001, 0.0001]
        upper_bounds = [200, 2]
        fixed_params = None
    else:
        demo_model = Demographics1D.snm
        param_names = []

    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    model_suffix = '' if args.model == 'bpg' else '.' + args.model
    result_file = args.out + f'.samplesize_{ns[0]}.{fold_state}{model_suffix}.1d_demo_fits.result'
    header = '\t'.join(['#ll_model'] + param_names + ['theta0']) + '\n'

    # Run multiple optimization starts with different seeds to find the best fit.
    n_starts = 1 if args.model == 'constant' else args.starts
    valid_starts = 0
    for seed in range(args.seed, args.seed + n_starts):
        try:
            if args.model == 'constant':
                popt = []
            else:
                numpy.random.seed(seed)
                p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                                            lower_bound=lower_bounds)
                popt, _ = dadi.Inference.opt(p0, data_fs, demo_model_ex, pts_l,
                                            lower_bound=lower_bounds,
                                            upper_bound=upper_bounds,
                                            fixed_params=fixed_params,
                                            algorithm=nlopt.LN_BOBYQA,
                                            maxeval=600, verbose=1)
            if not numpy.isfinite(popt).all():
                raise ValueError('Nonfinite demographic parameters.')

            model_fs = demo_model_ex(popt, ns, pts_l)
            ll_model = dadi.Inference.ll_multinom(model_fs, data_fs)
            theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)
            # Check for non-finite likelihood or theta0 values, and ensure theta0 is positive.
            if not numpy.isfinite([ll_model, theta0]).all() or theta0 <= 0:
                raise ValueError('Invalid likelihood or theta.')
            
        except Exception as e:
            print(f'Skipping seed {seed}: {e}')
            continue

        res = [ll_model] + list(popt) + [theta0]
        with open(result_file, 'a+') as f:
            # Write one header when SLURM jobs append to a common result file.
            fcntl.flock(f, fcntl.LOCK_EX)
            f.seek(0)
            existing_header = f.readline()
            if not existing_header:
                f.write(header)
            elif existing_header != header:
                raise ValueError(f'Unexpected header in {result_file}')
            f.write('\t'.join(map(str, res)) + '\n')
        valid_starts += 1

    # If no valid demographic fits were found, raise an error.
    if not valid_starts:
        raise RuntimeError('No valid demographic fits.')


if __name__ == "__main__":
    main()
