# Overview

If you have any questions, please feel free to contact me via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

- [Folder `COSMIC_sig`](#folder-cosmic_sig)
- [Folder `variant_count`](#folder-variant_count)
- [`hg38.3-mer-count.txt`](#hg383-mer-counttxt)
- [`human.strand_specific.mut_info.tsv`](#humanstrand_specificmut_infotsv)

## Folder `COSMIC_sig`

This folder contains mutation profiles derived from [COSMIC SBS signatures](https://cancer.sanger.ac.uk/signatures/sbs). SBS2 is replaced with the strand-specific SBS2ss* signature reported by [Liu et al.](https://doi.org/10.1038/s41586-024-07532-8). The three files correspond to different SBS combinations used for MuSiCal refitting.

## Folder `variant_count`

This folder contains mutation spectra calculated using `script/calculate_1-3-mer_mt.py` for different gene groups at the transcript level. Variants were filtered with allele frequency (AF) ≤ 5 × 10⁻³ in the downsampled gnomAD HGDP + 1KG cohort.

For each `.result` file, site- and variant-level information at both the 1-mer and 3-mer levels is provided, including counts for each variant type and the corresponding numbers of masked sites. An example command to generate a result file is shown below:

```bash
python script/calculate_1-3-mer_mt.py \
    --genome hg38.fa \
    -b hg38.gencodeV44.repli-dep-histone_transcript.bed \
    --max_af 0.005 \
    -n hg38.gencodeV44.repli-dep-histone_transcript_gnomAD_5e-3 \
    --output_path . \
    -v gnomad.genomes.v3.1.2.hgdp_tgp.sample_3k.snv.vcf.gz
```

## `hg38.3-mer-count.txt`

This file contains the 3-mer composition of the human genome (hg38) after applying coverage and mappability filters. It is used to scale mutation spectra for different gene groups.

### How to generate

Create a script named `get_kmer_component.py` with the following content:

```python
import argparse
from Bio import SeqIO
from collections import defaultdict
from tqdm import tqdm
from multiprocessing import Pool, cpu_count

def count_kmers(sequences, k):
    kmer_counts = defaultdict(int)
    for sequence in sequences:
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if 'N' in kmer:  # Skip k-mers containing ambiguous bases
                continue
            kmer_counts[kmer] += 1
    return kmer_counts

def process_chromosome(args):
    record, k, region_list = args
    if region_list is None:
        sequences = [str(record.seq).upper()]
    else:
        sequences = []
        for start, end in region_list:
            sequences.append(str(record.seq[start:end]).upper())
    return count_kmers(sequences, k)

def count_genome_kmers(genome_file, k, num_processes, region_file):
    genome_kmer_counts = defaultdict(int)
    records = list(SeqIO.parse(genome_file, "fasta"))

    target_region = {}
    if region_file is not None:
        with open(region_file, 'r') as f:
            for line in f:
                chrom, start, end = line.strip().split('\t')[:3]
                target_region.setdefault(chrom, []).append((int(start), int(end)))

    if target_region:
        records = [record for record in records if record.id in target_region]

    if num_processes is None:
        num_processes = cpu_count()

    with Pool(num_processes) as pool:
        results = pool.imap_unordered(
            process_chromosome,
            [(record, k, target_region.get(record.id)) for record in records]
        )
        for kmer_counts in tqdm(results, total=len(records), desc="Processing chromosomes"):
            for kmer, count in kmer_counts.items():
                genome_kmer_counts[kmer] += count

    return genome_kmer_counts

def main():
    parser = argparse.ArgumentParser(description="Count k-mers in a genome")
    parser.add_argument("-g", "--genome", required=True, help="Path to genome FASTA file")
    parser.add_argument("-k", "--kmer", type=int, required=True, help="Length of k-mer")
    parser.add_argument("-o", "--output", required=True, help="Path to output file")
    parser.add_argument("-p", "--processes", type=int, help="Number of processes to use (default: all available)")
    parser.add_argument("-r", "--region", help="Target regions in BED format (default: entire genome)")
    args = parser.parse_args()

    kmer_counts = count_genome_kmers(args.genome, args.kmer, args.processes, args.region)
    kmer_counts = dict(sorted(kmer_counts.items()))
    with open(args.output, 'w') as f:
        for kmer, count in kmer_counts.items():
            f.write(f"{kmer}\t{count}\n")

if __name__ == "__main__":
    main()
```

Run the script as follows:

```bash
python get_kmer_component.py \
    -g /path/to/hg38.fa \
    -k 3 -p 16 \
    -o hg38.3-mer-count.txt \
    --region /path/to/hg38.genmap.k150e1_unique.cov_16-48.bed
```

## `human.strand_specific.mut_info.tsv`

This file contains strand-specific variant and site information for different gene groups, generated using the `--contain_strand` option in `script/calculate_1-3-mer_mt.py`. The column structure is similar to that of `figures/Fig1/data/multispecies_4dsite.1_mer_record.tsv`.