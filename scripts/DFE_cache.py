import dadi.DFE as DFE
import pickle
import glob
import json
import os
import dadi
import numpy
import pandas as pd
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
dadi.Integration.timescale_factor = 1e-6

def get_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate a 1D DFE cache with dadi")
    
    input_parser = parser.add_argument_group("Input")
    fs_parser = input_parser.add_mutually_exclusive_group(required=True)
    fs_parser.add_argument("--name", type=str,
                        help="Input file prefix for 1D DFE cache generation.")
    fs_parser.add_argument("--fs", type=str,
                        help="Input projected 1D spectrum (.fs).")
    input_parser.add_argument("--opt_result", type=str, required=True,
                        help="Demographic optimization result (.result or demography.json).")
    
    param_parser = parser.add_argument_group("Parameters")
    param_parser.add_argument("--model", choices=['bpg', 'two_epoch', 'constant'], default='bpg',
                        help="Demographic model used for the optimization. Default: bpg")
    param_parser.add_argument("--cpu", type=int, default=32,
                        help="Number of CPUs to use. Default: 32")

    output_parser = parser.add_argument_group("Output")
    output_parser.add_argument("-o", "--out", type=str,
                        help="Output cache file. By default, write next to the input spectrum.")

    return parser.parse_args()


def main():

    args = get_args()

    fs_file = args.fs
    if fs_file is None:
        fs_files = glob.glob(f'{args.name}*.fs')
        fs_file = fs_files[0]

    data_fs = dadi.Spectrum.from_file(fs_file)
    ns = data_fs.sample_sizes

    def bottleneck_plateau_growth_sel(params, ns, pts):

        nu1, nu2, nuc, T1, T2, Tc, gamma = params
        xx  = Numerics.default_grid(pts)
        phi = PhiManip.phi_1D(xx, gamma=gamma)   
          
        phi = Integration.one_pop(phi, xx, T1, nu1, gamma=gamma)
        phi = Integration.one_pop(phi, xx, T2, nu2, gamma=gamma)
        nu_func = lambda t: nu2 * numpy.exp(numpy.log(nuc/nu2) * (t / Tc))
        phi = Integration.one_pop(phi, xx, Tc, nu_func, gamma=gamma)

        fs = Spectrum.from_phi(phi, ns, (xx,))
        return fs

    def constant_sel(params, ns, pts):

        gamma = params[-1]
        xx  = Numerics.default_grid(pts)
        phi = PhiManip.phi_1D(xx, gamma=gamma)

        fs = Spectrum.from_phi(phi, ns, (xx,))
        return fs

    pts_l = [max(ns)+140, max(ns)+150, max(ns)+160]
    if args.model == 'bpg':
        demo_sel_model = bottleneck_plateau_growth_sel
        param_names = ['nu1', 'nu2', 'nuc', 'T1', 'T2', 'Tc']
    elif args.model == 'two_epoch':
        demo_sel_model = DFE.DemogSelModels.two_epoch_sel
        param_names = ['nu', 'T']
    else:
        demo_sel_model = constant_sel
        param_names = []

    if args.opt_result.endswith('.json'):
        with open(args.opt_result, 'r') as f:
            opti_result = json.load(f)
        if not opti_result['valid']:
            raise ValueError('The demographic fit in the JSON file is invalid.')
        popt = opti_result['params']
    else:
        opti_result_df = pd.read_csv(args.opt_result, sep='\t', float_precision='round_trip')
        fit_columns = ['#ll_model'] + param_names + ['theta0']
        opti_result_df = opti_result_df[numpy.isfinite(opti_result_df[fit_columns]).all(axis=1)]

        # Get maximum likelihood parameters
        opti_result_df = opti_result_df.sort_values(by='#ll_model', ascending=False, kind='stable')
        popt = opti_result_df.iloc[0][param_names].to_list()

    print('running...')
    print('Best params:')
    print(popt)

    cache1d = DFE.Cache1D(popt, ns, demo_sel_model, pts=pts_l, gamma_bounds=[1e-4, 7500], gamma_pts=300, cpus=args.cpu, verbose=True)
    if (cache1d.spectra<0).sum() > 0:
        print(
            '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(cache1d.spectra.min())+
            '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
            )
        if cache1d.spectra.min() < -0.001:
            raise RuntimeError('Large negative values in the DFE cache.')
        
    # Save cache
    model_suffix = '' if args.model == 'bpg' else '_' + args.model
    cache_file = args.out or fs_file + model_suffix + '_1d_cache.bpkl'
    tmp_file = cache_file + f'.{os.getpid()}.tmp'
    try:
        with open(tmp_file, 'wb') as fid:
            pickle.dump(cache1d, fid, protocol=2)
        os.replace(tmp_file, cache_file)
    finally:
        if os.path.exists(tmp_file):
            os.remove(tmp_file)

if __name__ == "__main__":
    main()
