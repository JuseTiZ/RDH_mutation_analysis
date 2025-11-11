import argparse
import os
from Bio import SeqIO
from tqdm import tqdm
import json
import copy
from cyvcf2 import VCF

def get_args():

    parser = argparse.ArgumentParser(description="Calculate the mutation spectrum from given BED files.")
    
    input_parser = parser.add_argument_group("Input parameters")
    input_parser.add_argument("--genome", required=True,
                        help="The path to genome sequence file.")
    input_parser.add_argument("-b", "--bed", nargs='+', required=True,
                        help="The bed files for retaining mutations.")
    input_parser.add_argument("-n", "--name", nargs='+',
                        help="The output name of files.")
    input_parser.add_argument("-v", "--vcf",
                        help="The variant file to count variants.")
    
    params_parser = parser.add_argument_group("Parameters")
    params_parser.add_argument("--max_af", type=float, required=True)
    params_parser.add_argument("--no_filter", action='store_true')
    params_parser.add_argument("--no_singleton", action='store_true')
    params_parser.add_argument("--contain_strand", action='store_true',
                        help="Use to combine strand information.")
    
    output_parser = parser.add_argument_group("Output parameters")
    output_parser.add_argument("--output_path", default='.',
                        help="The output path of files.")

    args = parser.parse_args()
    return args


def main():

    args = get_args()
    # If no name provided, use file basename
    if not args.name:
        args.name = [os.path.basename(i) for i in args.bed]
    
    if len(args.name) != len(args.bed):
        print('Error: The number of names should be equal to the number of bed files.')
        print(args.name, args.bed)
        return
    
    print('Loading genome...')
    seqdic = SeqIO.to_dict(SeqIO.parse(args.genome, 'fasta'))

    # Load variant file
    vcf_file = VCF(args.vcf)

    nuc_list = ['A', 'T', 'C', 'G']
    nuc_set = set(nuc_list)
    # Get 3-mer list
    three_mer_list = [a+b+c for a in nuc_list for b in nuc_list for c in nuc_list]
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    # Initialize k-mer dict
    codestrand_3mer_info_init = {k: 0 for k in three_mer_list}
    codestrand_1mer_info_init = {
            'A': 0,
            'T': 0,
            'C': 0,
            'G': 0,
            'CpG': 0
        }
    # Initialize k-mer mutation dict
    three_mer_variant_info_init = {f'{a}>{a[0]}{b}{a[-1]}': {
        "mask": 0,
        "count": 0
    } for a in three_mer_list for b in nuc_list if a[1] != b}
    one_mer_variant_info_init = {f'{a}>{b}': {
        "mask": 0,
        "count": 0
    } for a in nuc_list for b in nuc_list if a != b}
    # For 1-mer, add CpG context
    for k in list(one_mer_variant_info_init.keys()):
        if k.startswith('C') or k.startswith('G'):
            one_mer_variant_info_init[f'{k[0]}(CpG)>{k[-1]}'] = {
                "mask": 0,
                "count": 0
            }

    for bed, name in zip(args.bed, args.name):

        codestrand_1mer_info = copy.deepcopy(codestrand_1mer_info_init)
        one_mer_variant_info = copy.deepcopy(one_mer_variant_info_init)

        codestrand_3mer_info = copy.deepcopy(codestrand_3mer_info_init)
        three_mer_variant_info = copy.deepcopy(three_mer_variant_info_init)

        with open(bed, 'r') as f:
            lines = f.readlines()

            for line in tqdm(lines, desc=f'Processing {name}', unit='line'):
                
                info = line.strip().split('\t')
                # Get strand information
                if args.contain_strand:
                    chr, ss, es, sd = info[0:4]
                else:
                    chr, ss, es = info[0:3]
                    sd = '+'

                # Add flanking nucleotide for CpG judgement
                region_seq = seqdic[chr].seq[(int(ss) - 1): (int(es) + 1)]
                if sd == '+':
                    region_seq = region_seq.upper()
                else:
                    region_seq = region_seq.reverse_complement().upper()

                # Get 3-mer component
                for i in range(len(region_seq) - 2):
                    three_mer = str(region_seq[i:i+3])
                    if three_mer in three_mer_list:
                        codestrand_3mer_info[three_mer] += 1

                # Get CpG for 1-mer component
                codestrand_1mer_info['CpG'] += region_seq.count('CG') * 2

                # Remove flanking CpG nucleotide
                if region_seq.startswith('CG'):
                    codestrand_1mer_info['CpG'] -= 1
                    codestrand_1mer_info[region_seq[1]] -= 1
                if region_seq.endswith('CG'):
                    codestrand_1mer_info['CpG'] -= 1
                    codestrand_1mer_info[region_seq[-2]] -= 1
                region_seq = str(region_seq[1:-1]).replace('CG', '')

                # Get 1-mer component
                for nuc in nuc_list:
                    codestrand_1mer_info[nuc] += region_seq.count(nuc)

                region = f"{chr}:{int(ss)+1}-{int(es)}"
                # Search variants in this region
                for record in vcf_file(region):
                    
                    # Get variant information
                    pos = record.POS
                    ref, alt = record.REF, record.ALT[0]
                    try:
                        ac, an = record.INFO['AC'], record.INFO['AN']
                        af = ac / an
                        if args.no_singleton and ac == 1:
                            continue
                    except KeyError:
                        af = 0
                
                    if args.no_filter:
                        filter = None
                    else:
                        filter = record.FILTER

                    mask = (af > args.max_af) or (filter is not None)

                    # If mask then add to mask count
                    add_item = 'mask' if mask else 'count'

                    # Skip if not A, T, C, G
                    if ref not in nuc_set or alt not in nuc_set:
                        continue

                    pos = pos - 1
                    mt_three_mer = seqdic[chr].seq[(int(pos) - 1): (int(pos) + 2)].upper()
                    # Get strand mutation
                    if sd == '-':
                        ref = complement[ref]
                        alt = complement[alt]
                        mt_three_mer = mt_three_mer.reverse_complement()

                    if 'N' not in mt_three_mer:
                        three_mer_variant_info[f'{mt_three_mer}>{mt_three_mer[0]}{alt}{mt_three_mer[-1]}'][add_item] += 1
                    else:
                        print(f'N found in {chr}:{pos+1} 3-mer {mt_three_mer}, skip.')

                    # Check if is CpG context
                    if 'CG' in mt_three_mer:
                        variant_key = f'{ref}(CpG)>{alt}'
                        one_mer_variant_info[variant_key][add_item] += 1
                        continue

                    variant_key = f'{ref}>{alt}'
                    one_mer_variant_info[variant_key][add_item] += 1

        # Save result
        result_file = os.path.join(args.output_path, name + '.result')
        result_dict = {
            '1mer': {
                'context': codestrand_1mer_info,
                'variant': one_mer_variant_info,
            }, 
            '3mer': {
                'context': codestrand_3mer_info,
                'variant': three_mer_variant_info,
            }
        }
        with open(result_file, 'w') as o:
            json.dump(result_dict, o)


if __name__ == '__main__':
    main()