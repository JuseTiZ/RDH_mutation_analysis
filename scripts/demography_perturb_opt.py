import nlopt
import dadi
import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
import argparse
import pickle
import os
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
    
    param_parser = parser.add_argument_group("Parameters")
    param_parser.add_argument("-s", "--seed", type=int, default=42,
                        help="Seed for perturb parameters.")
    param_parser.add_argument("--popid", type=str,
                        help="Population ID for infer 1D demography.")
    param_parser.add_argument("--ns", type=int,
                        help="Sample size to project down.")
    param_parser.add_argument("--fold_fs", action='store_true',
                        help="Fold the frequency spectrum.")
    
    output_parser = parser.add_argument_group("Output")
    output_parser.add_argument("-o", "--out", type=str,
                        help="Output file name for demographic inference results.")
    
    return parser.parse_args()


def main():

    args = get_args()
    numpy.random.seed(args.seed)

    fold_state = 'unfold' if not args.fold_fs else 'fold'
    polarized = False if args.fold_fs else True
    fs_file = args.out + f'.samplesize_{args.ns}.{fold_state}.fs'
    if not os.path.exists(fs_file):
        if args.data_dict is not None:
            print("Using data dictionary from {}".format(args.data_dict))
            syn_dd = pickle.load(open(args.data_dict, 'rb'))
        else:
            print("Generating data dictionary from VCF: {}".format(args.vcf))
            syn_dd = dadi.Misc.make_data_dict_vcf(args.vcf, args.popinfo)
        
        pop_ids, ns = [args.popid], [args.ns]
        data_fs = dadi.Spectrum.from_data_dict(syn_dd, pop_ids, ns, polarized=polarized)
        # Save the frequency spectrum
        data_fs.to_file(fs_file)
        print("Frequency spectrum saved to {}".format(fs_file))

    else:
        print("Loading frequency spectrum from {}".format(fs_file))
        data_fs = dadi.Spectrum.from_file(fs_file)
        ns = data_fs.sample_sizes

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

    demo_model = bottleneck_plateau_growth
    # Wrap the demographic model in a function that utilizes grid points which increases dadi's ability to more accurately generate a model frequency spectrum.
    demo_model_ex = dadi.Numerics.make_extrap_func(demo_model)

    # nu1, nu2, nuc, T1, T2, Tc
    params = [0.1, 1, 10, 0.01, 0.05, 0.05]
    lower_bounds = [1e-2, 1e-2, 1, 1e-3, 1e-3, 1e-3]
    upper_bounds = [5, 5, 200, 1, 1, 1]

    p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bounds,
                                lower_bound=lower_bounds)
    popt, ll_model = dadi.Inference.opt(p0, data_fs, demo_model_ex, pts_l,
                                lower_bound=lower_bounds,
                                upper_bound=upper_bounds,
                                fixed_params=[None, None, None, 0.01, None, None],
                                algorithm=nlopt.LN_BOBYQA,
                                maxeval=600, verbose=1)

    model_fs = demo_model_ex(popt, ns, pts_l)
    theta0 = dadi.Inference.optimal_sfs_scaling(model_fs, data_fs)
    res = [ll_model] + list(popt) + [theta0]

    result_file = args.out + f'.samplesize_{args.ns}.{fold_state}.1d_demo_fits.result'
    if not os.path.exists(result_file):
        with open(result_file, 'w') as f:
            f.write('#ll_model\tnu1\tnu2\tnuc\tT1\tT2\tTc\ttheta0\n')
    else:
        print(f"Appending to existing result file: {result_file}")

    with open(result_file, 'a') as f:
        f.write('\t'.join(map(str, res)) + '\n')


if __name__ == "__main__":
    main()