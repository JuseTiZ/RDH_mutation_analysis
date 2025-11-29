import argparse
from tqdm import tqdm
from calcu_transcript_feature import parse_gtf
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import gzip
import warnings

def get_args():

    parser = argparse.ArgumentParser()

    input_group = parser.add_argument_group(title="Input")
    input_group.add_argument("-g", "--gtffile", type=str, required=True)
    input_group.add_argument("--genome", type=str, required=True,)

    params_group = parser.add_argument_group(title="Parameters")
    params_group.add_argument("--ensembl_gtf", action='store_true',
                              help="Whether the GTF file is in Ensembl format. Default is False.")
    params_group.add_argument("--chr_prefix", action='store_true',
                              help="Whether to add 'chr' prefix to chromosome names. Default is False.")
    params_group.add_argument("--codon_table", default="Standard", type=str,
                              help="The codon table to use. Default is 'standard'.")
    params_group.add_argument("--valid_pc", default=None, 
                              help="Protein sequence to valid. Default is None.")
    params_group.add_argument("--delimiter", default='|', 
                              help="Delimiter to decide transcript ID. Default is '|'.")
    params_group.add_argument("--transcript_idx", default=1, type=int)
    params_group.add_argument("--output_pep", action='store_true')

    output_group = parser.add_argument_group(title="Output")
    output_group.add_argument("-o", "--output", type=str, required=True,
                              help="The path of the output 4dsite file.")
                            
    args = parser.parse_args()
    return args


def parse_codon_4dsite(
        codon_table: dict
):

    all_nuc_set = set(['A', 'T', 'C', 'G'])
    codon_4dsite_dic = {}

    for codon, aa in codon_table.items():
        for idx, nuc in enumerate(codon):
            is_4d_site = True
            for mut_nuc in all_nuc_set - {nuc}:
                new_codon = codon[:idx] + mut_nuc + codon[idx + 1:]
                if codon_table.get(new_codon) != aa:
                    is_4d_site = False
                    break
            if is_4d_site:
                if codon not in codon_4dsite_dic:
                    codon_4dsite_dic[codon] = []
                codon_4dsite_dic[codon].append(idx)
    
    return codon_4dsite_dic


def check_cds_order(cds_parts, strand):

    order = True
    ex_start = None
    ex_end = None
    for _, ss, es in cds_parts:
        if ex_start is None:
            ex_start = ss
            ex_end = es
            continue
        if strand == '+':
            if ss < ex_end:
                order = False
                break
        elif strand == '-':
            if es > ex_start:
                order = False
                break
    
    if order == False:
        warnings.warn(f"CDS parts are not in order for strand {strand}.\nMay due to sorted gtf or overlapped CDS")


def main():

    args = get_args()

    # Get transcript information from GTF file
    transcript_info_dic = parse_gtf(args.gtffile, region_list=['CDS'], ensembl_gtf=args.ensembl_gtf)
    # Get the codon table and parse 4D sites
    print(f"Parsing codon table... (using codon table: {args.codon_table})")
    codon_table = CodonTable.unambiguous_dna_by_name[args.codon_table]
    codon_4dsite_dic = parse_codon_4dsite(codon_table.forward_table)
    # Load genome
    print("Loading genome...")
    genome = SeqIO.to_dict(SeqIO.parse(args.genome, "fasta"))
    # Prepare non-4d site and 4d site set
    chrom_site_dict = {}

    if args.valid_pc is not None:
        print("Loading valid protein coding sequence...")
        handle = gzip.open(args.valid_pc, 'rt') if args.valid_pc.endswith('.gz') else open(args.valid_pc, 'r')
        pc_fasta = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
        pc_fasta = {
            k.split(args.delimiter)[args.transcript_idx]: v for k, v in pc_fasta.items()
        }
        unmatch_df_lst = []

    if args.output_pep:
        transcript_id_pep_dict = {}
        transcript_id_cds_dict = {}

    for transcript_id, info in tqdm(transcript_info_dic.items(), unit=" transcripts"):

        if 'CDS' not in info['region']:
            continue

        cds_parts = info['CDS']
        cds_strand = info['strand']
        check_cds_order(cds_parts, cds_strand)
        cds_parts = sorted(cds_parts, key=lambda x: x[1])

        chrom = cds_parts[0][0]
        # Add chromosome prefix if needed
        if args.chr_prefix:
            chrom = 'chr' + chrom

        # Skip if chromosome not in reference
        if chrom not in genome:
            continue

        if chrom not in chrom_site_dict:
            chrom_site_dict[chrom] = {
                "4d_site": set(),
                "non_4d_site": set()
            }

        # Get the sequence of the CDS
        genome_site = []
        cds_sequence = ''
        for _, start, end in cds_parts:
            cds_sequence += genome[chrom].seq[start:end]
            genome_site += range(start, end)
        cds_sequence = cds_sequence.upper()
        if cds_strand == '-':
            # Reverse complement the CDS sequence if on the negative strand
            genome_site = list(reversed(genome_site))
            cds_sequence = Seq(cds_sequence).reverse_complement()

        # Check if the CDS sequence length is a multiple of 3
        if len(cds_sequence) % 3 != 0:
            continue

        # Translate the CDS sequence to protein sequence
        try:
            protein_sequence = cds_sequence.translate(table=codon_table)
        except Exception as e:
            # If translation fails, skip this transcript
            print(f"Error translating {transcript_id}: {e}")
            continue

        # Save protein sequence if needed
        if args.output_pep:
            transcript_id_pep_dict[transcript_id] = str(protein_sequence)
            transcript_id_cds_dict[transcript_id] = cds_sequence

        # Valid if translated sequence match with annotation
        if args.valid_pc is not None:
            protein_sequence = protein_sequence.replace('*', '') # For Sec
            if transcript_id in pc_fasta:
                ref_pc_seq = pc_fasta[transcript_id].seq
                ref_pc_seq = ref_pc_seq.replace('X', '').replace('U', '') # For Sec and unknown
                # exclude non_ATG_start
                if protein_sequence[1:] != ref_pc_seq[1:]:
                    unmatch_df_lst.append({
                        'transcript_id': transcript_id,
                        'reference': str(ref_pc_seq),
                        'translate': str(protein_sequence),
                    })
                    continue
                    # raise ValueError(f'Sequence mismatch found. {transcript_id}\nREFERENCE: {str(ref_pc_seq)}\nTRANSLATE: {str(protein_sequence)}')
            else:
                unmatch_df_lst.append({
                    'transcript_id': transcript_id,
                    'reference': 'not found',
                    'translate': str(protein_sequence),
                })
                print(f'{transcript_id} not found in valid sequence.')
                
        # Save 4dsite position
        for pos in range(0, len(cds_sequence), 3):
            codon = cds_sequence[pos:pos + 3]
            genome_site_pos = genome_site[pos:pos + 3]
            codon_4dsite = codon_4dsite_dic.get(codon, [])
            for idx, genome_pos in enumerate(genome_site_pos):
                if idx not in codon_4dsite:
                    chrom_site_dict[chrom]["non_4d_site"].add(genome_pos)
                else:
                    chrom_site_dict[chrom]["4d_site"].add(genome_pos)

    # Write absolute 4D site to the output file
    print("Writing absolute 4D sites to output file...")
    with open(args.output, 'w') as output_file:
        for chrom, sites in chrom_site_dict.items():
            absolute_4d_sites = sites["4d_site"] - sites["non_4d_site"]
            for site in sorted(list(absolute_4d_sites)):
                output_file.write(f"{chrom}\t{site}\t{site+1}\n")
    
    # Save information of unmatched transcripts
    if args.valid_pc is not None:
        unmatch_df = pd.DataFrame(unmatch_df_lst)
        unmatch_df.to_csv(f'{args.output}.unmatch_transcript.tsv', index=False, sep='\t')

    # Save protein sequences if needed
    if args.output_pep:
        with open(f'{args.output}.pep.fasta', 'w') as pep_file:
            for transcript_id, protein_sequence in transcript_id_pep_dict.items():
                pep_file.write(f">{transcript_id}\n{protein_sequence}\n")
        with open(f'{args.output}.cds.fasta', 'w') as cds_file:
            for transcript_id, cds_sequence in transcript_id_cds_dict.items():
                cds_file.write(f">{transcript_id}\n{cds_sequence}\n")

if __name__ == "__main__":
    main()