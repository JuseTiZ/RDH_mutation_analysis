import dadi.DFE as DFE
import pickle
import glob
import dadi
import numpy
import pandas as pd
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
dadi.Integration.timescale_factor = 1e-6

def get_args():
    import argparse

    parser = argparse.ArgumentParser(
        description="Demographic inference with dadi")
    
    input_parser = parser.add_argument_group("Input")
    input_parser.add_argument("--name", type=str,
                        help="Input file prefix for 1D DFE cache generation.")
    input_parser.add_argument("--opt_result", type=str,
                        help="Optimization result file for 1D demography.")
    
    param_parser = parser.add_argument_group("Parameters")
    param_parser.add_argument("--cpu", type=int, default=32,
                        help="Number of CPUs to use. Default: 32")

    return parser.parse_args()


def main():

    args = get_args()

    fs_file = glob.glob(f'{args.name}*.fs')
    if len(fs_file) != 1:
        print(f'Error: Found {len(fs_file)} fs files matching {args.name}*.fs')
        return
    fs_file = fs_file[0]

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

    pts_l = [max(ns)+140, max(ns)+150, max(ns)+160]
    demo_sel_model = bottleneck_plateau_growth_sel

    opti_result_df = pd.read_csv(args.opt_result, sep='\t')
    # Get maximum likelihood parameters
    opti_result_df = opti_result_df.sort_values(by='#ll_model', ascending=False)
    popt = opti_result_df.iloc[0, 1:7].to_list()

    print('running...')
    print('Best params:')
    print(popt)

    cache1d = DFE.Cache1D(popt, ns, demo_sel_model, pts=pts_l, gamma_bounds=[1e-5, 500], gamma_pts=300, cpus=args.cpu, verbose=True)
    if (cache1d.spectra<0).sum() > 0:
        print(
            '!!!WARNING!!!\nPotentially large negative values!\nMost negative value is: '+str(cache1d.spectra.min())+
            '\nIf negative values are very negative (<-0.001), rerun with larger values for pts_l'
            )
    # Save cache
    fid = open(fs_file+'_1d_cache.bpkl', 'wb')
    pickle.dump(cache1d, fid, protocol=2)
    fid.close()

if __name__ == "__main__":
    main()