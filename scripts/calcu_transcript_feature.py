import pandas as pd
import numpy as np
import re
import argparse
import gzip
from tqdm import tqdm
import pyBigWig
import warnings

def get_args():

    parser = argparse.ArgumentParser()

    input_group = parser.add_argument_group(title="Input")
    input_group.add_argument("-g", "--gtffile", type=str, required=True)
    input_group.add_argument("-bw", "--bigwig", nargs="+", required=True,
                             help="The path of bigwig file(s) for each sample.")
    input_group.add_argument("--bigwigname", nargs="+", required=True,
                             help="The name of each bigwig file, used for legend.")
    
    option_group = parser.add_argument_group(title="Optional")
    option_group.add_argument("--include_genetype", nargs="*",
                             help="The gene types to include in the analysis, e.g., protein_coding. Default is all genes.")
    option_group.add_argument("--include_region", nargs="+", default=['transcript', 'exon', 'CDS'],
                             help="The regions to include in the analysis. Default is ['transcript', 'exon', 'CDS'].")
    option_group.add_argument("--ensembl_gtf", action='store_true',
                              help="Whether the GTF file is in Ensembl format. Default is False.")


    output_group = parser.add_argument_group(title="Output")
    output_group.add_argument("-o", "--output", type=str, required=True,
                              help="The path of the output file.")
                            
    args = parser.parse_args()
    return args


def parse_gtf(
        gtffile, 
        gene_type_list=None,
        region_list=None,
        ensembl_gtf=False
        ):

    transcript_info_dic = {}
    print("Parsing GTF......")

    gene_type_str = 'gene_type' if not ensembl_gtf else 'gene_biotype'

    gtffile_reader = gzip.open(gtffile, 'rt') if gtffile.endswith('.gz') else open(gtffile, 'r')
    lines = gtffile_reader.readlines()
    for line in tqdm(lines):

        if line.startswith('#'):
            continue

        fields = line.strip().split('\t')
        region = fields[2]

        try:
            gene_type = re.search(rf'{gene_type_str} "(.*?)";', line).group(1)
            if gene_type_list and gene_type not in gene_type_list:
                continue
        except AttributeError:
            gene_type = None

        if region_list and region not in region_list:
            continue

        transcript_id = re.search(r'transcript_id "(.*?)";', line)
        if not transcript_id:
            continue
        transcript_id = transcript_id.group(1)

        if transcript_id not in transcript_info_dic:

            gene_name = re.search(r'gene_name "(.*?)";', line)
            gene_name = None if not gene_name else gene_name.group(1)
            
            gene_id = re.search(r'gene_id "(.*?)";', line).group(1)
            transcript_info_dic[transcript_id] = {
                'gene_name': gene_name,
                'gene_id': gene_id,
                'gene_type': gene_type,
                'strand': fields[6],
                'region': []
            }

        if region not in transcript_info_dic[transcript_id]:
            transcript_info_dic[transcript_id][region] = []
            transcript_info_dic[transcript_id]['region'].append(region)
        
        chrom = fields[0]
        start = int(fields[3]) - 1
        end = int(fields[4])
        transcript_info_dic[transcript_id][region].append((chrom, start, end))
    
    gtffile_reader.close()

    if gene_type is None:
        warnings.warn("Gene type is not found in the GTF file. Use all genes.")

    return transcript_info_dic


def main():
    args = get_args()

    bigwig_files = args.bigwig
    bigwig_names = args.bigwigname
    if len(bigwig_files) != len(bigwig_names):
        raise ValueError("The number of bigwig files and names must match.")

    print(f"GTF file: {args.gtffile}")
    print(f"Included gene types: {args.include_genetype if args.include_genetype else 'all'}")
    transcript_info_dic = parse_gtf(
        args.gtffile, args.include_genetype, args.include_region, args.ensembl_gtf
    )
    total_transcript_num = len(transcript_info_dic)

    final_df_lst = []

    bar = tqdm(total=len(bigwig_files), desc="Processing BigWig files")
    for bigwig_file, bigwig_name in zip(bigwig_files, bigwig_names):

        bar.set_description(f"Processing {bigwig_name}")
        bigwig = pyBigWig.open(bigwig_file)
        if not list(bigwig.chroms().keys())[0].startswith('chr'):
            strip_chr = True
        else:
            strip_chr = False
        bar.set_postfix({
            'progress': f"0 / {total_transcript_num}"
        })
        error_chrom = set()

        for idx, (transcript_id, info) in enumerate(transcript_info_dic.items()):

            for region in info['region']:
                region_value = []
                try:
                    for chrom, start, end in info[region]:
                        chrom = chrom[3:] if strip_chr else chrom
                        values = bigwig.values(chrom, start, end)
                        region_value.append(values)
                except:
                    error_chrom.add(info[region][0][0])
                    continue

                if len(region_value) == 0:
                    continue

                region_value = np.concatenate(region_value)
                # Calculate the nan percentage and mean value
                nan_value_num = np.sum(np.isnan(region_value))
                non_nan_length = len(region_value) - nan_value_num
                nan_percentage = nan_value_num / len(region_value)
                if nan_percentage == 1.0:
                    mean_value = np.nan
                else:
                    mean_value = np.nanmean(region_value)

                final_df_lst.append({
                    'gene id': info['gene_id'],
                    'gene name': info['gene_name'],
                    'transcript id': transcript_id,
                    'gene type': info['gene_type'],
                    'region': region,
                    'bigwig name': bigwig_name,
                    'mean value': mean_value,
                    'non-nan length': non_nan_length,
                    'nan percentage': nan_percentage,
                })
            
            if idx % 100 == 0:
                bar.set_postfix({
                    'progress': f"{idx} / {total_transcript_num}"
                })
        
        print(f"Error chromosomes in {bigwig_name}: {list(error_chrom)}")
        bigwig.close()
        bar.update(1)

    final_df = pd.DataFrame(final_df_lst)
    final_df.to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
