# Overview

If you have any questions, please feel free to contact us via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

- [Folder `nonsyn_SFSs_DFEs`](#folder-nonsyn_sfss_dfes)
- [`1KG_sanders.trio_TRD_analysis.tsv`](#1kg_sanderstrio_trd_analysistsv)
- [`Neville2025_*`](#neville2025_)
  
## Folder `nonsyn_SFSs_DFEs`

This folder contains **site frequency spectra (SFSs)** and the **estimated parameters of the distribution of fitness effects (DFEs)** for different gene groups.

### `gnomadv4.nfe.nonsynonymous.folded.NS2000.*.fs`

These files contain **folded SFSs of nonsynonymous variants** after downsampling using `dadi` (`dadi.Spectrum.from_data_dict`).

The file format follows the specification described in the [`dadi`](https://github.com/RyanGutenkunst/dadi/blob/33eed8cae4cd09df0b8cf49051c67f5da2aa3de9/dadi/Spectrum_mod.py#L281):

1. A single line containing the dimension of the fs array.
   - On the *same line*, the string `folded` or `unfolded` indicating the folding status.
   - Optional population labels in quotes (e.g., `"pop1"`).

2. A single line containing the array elements.

3. A single line containing the mask values corresponding to the array elements  
   (`1` = masked, `0` = unmasked).

#### How to generate these files

If genotype information is **not available**, the following function can be used to construct a `data_dict` as input for `dadi.Spectrum.from_data_dict`:

```python
from cyvcf2 import VCF

# Refer to `dadi.Misc.make_data_dict_vcf` in dadi
def read_vcf_return_dd(
        vcf: str,
        ac_tag: str = "AC",
        an_tag: str = "AN",
        popid: str = "EUR",
):

    dd = {}
    vcf = VCF(vcf)

    for variant in vcf:
        chrom, pos, ref, alt = variant.CHROM, variant.POS, variant.REF, variant.ALT[0]

        ac = variant.INFO.get(ac_tag)
        an = variant.INFO.get(an_tag)

        variant_key = f'{chrom}:{pos}:{ref}:{alt}'
        dd[variant_key] = {
            'segregating': (ref, alt),
            'context': f'-{ref}-',
            'outgroup_allele': '-',
            'outgroup_context': '---',
            'calls': {popid: (an-ac, ac)}
        }

    return dd
````

Example workflow:

```python
import dadi

# Generate data_dict using the nonsynonymous VCF of a specific gene group
nonsyn_dd = read_vcf_return_dd(nonsyn_vcf, ac_tag=AC_TAG, an_tag=AN_TAG, popid=POPID)

# Downsample and generate folded SFS
data_fs_folded_nonsyn = dadi.Spectrum.from_data_dict(
    nonsyn_dd, [POPID], [ns], polarized=False
)

# Save SFS
data_fs_folded_nonsyn.to_file(OUT_SFS)
```

For VCF files **with genotype information** (e.g., 1KG), `dadi.Misc.make_data_dict_vcf` and `dadi.Spectrum.from_data_dict` in `dadi` can be used directly to generate the data dictionary and SFS.

### `gnomadv4.nfe.folded.NS2000.*_1d_dfe_fits.result`

These files contain **estimated parameters of the DFE**, assuming a **gamma distribution of 2Ns**.

Columns:

| Column | Description                                 |
| ------ | ------------------------------------------- |
| 1      | Log-likelihood                              |
| 2      | Gamma shape parameter (α)                   |
| 3      | Gamma scale parameter (β)                   |
| 4      | θ<sub>ns</sub> used in `dadi.Inference.opt` |

#### How to generate these files

The full procedure for DFE inference is described in detail in the paper. Briefly:

1. Infer the demographic model using the **synonymous SFS of the whole CDS** (`scripts/demography_perturb_opt.py`).
2. Generate SFSs under different levels of selection (`scripts/DFE_cache.py`).
3. For each gene group:

   * Estimate θ<sub>s</sub> using the demographic model and synonymous SFS.
   * Scale θ<sub>s</sub> to θ<sub>ns</sub> using the site-count ratio calculated from **MuRaL-derived mutation rates**.
4. Assume a **gamma probability density function for 2Ns** and infer the parameters.

Following the `basic_workflow` in the `dadi` repository is recommended as a starting point.

## `1KG_sanders.trio_TRD_analysis.tsv`

This file contains **filtered variant records for trios** from the 1KG and Sanders exome datasets.

| Column                      | Description                                                                                                                                          |
| --------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------- |
| `child`, `father`, `mother` | Genotypes of the trio at this variant                                                                                                                |
| `site`                      | Variant information (`chromosome:position:ref:alt`)                                                                                                  |
| `consequence`               | Coding consequence                                                                                                                                   |
| `gene name`                 | Influencing gene(s)                                                                                                                                  |
| `child id`                  | Child ID of the trio                                                                                                                                 |
| `gnomad AF`                 | Allele frequency of this variant in gnomAD                                                                                                           |
| `dataset`                   | Dataset and gene group (`[dataset]_[gene group]`)                                                                                                    |
| `ref AF`                    | gnomAD AF of the reference allele, calculated as `1 − sum(AF of all SNVs at this site)`. Minor deviations may occur at sites overlapping with InDels |

### Data sources

* **1KG:** gnomAD HGDP + 1KG callset
  [https://gnomad.broadinstitute.org/data#v3-hgdp-1kg](https://gnomad.broadinstitute.org/data#v3-hgdp-1kg)

* **Sanders exome:**
  [https://www.ncbi.nlm.nih.gov/bioproject/PRJNA167318/](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA167318/)

* **gnomAD SNV allele frequencies:**
  [https://gnomad.broadinstitute.org/data#v4](https://gnomad.broadinstitute.org/data#v4)

## `Neville2025_*`

### `Neville2025_SupTable7_SpermExomeTargVars.hg38.vcf`

Variants (hg19 coordinates) were obtained from **“Supplementary Table 7 | Sperm Exome and Targeted NanoSeq Variants”** in the study by [Neville et al.](https://doi.org/10.1038/s41586-025-09448-3). Coordinates were lifted over to hg38 using **CrossMap**.

### `Neville2025_SupTable7_SpermExomeTargVars.hg38.refGene.exonic_variant_function`

ANNOVAR annotation results (GENCODE v44) for `Neville2025_SupTable7_SpermExomeTargVars.hg38.vcf`.

### `Neville2025_sperm_mutation_count.with_site.tsv`

Gene-level mutation count summary derived from  
`Neville2025_SupTable7_SpermExomeTargVars.hg38.refGene.exonic_variant_function`, with additional gene-group annotations and three columns calculated from [MuRaL](https://github.com/CaiLiLab/MuRaL)-predicted mutation rates:

| Column      | Description                                                      |
|-------------|------------------------------------------------------------------|
| `non_nsite` | Sum of mutation rates for all possible nonsense SNVs in the gene |
| `mis_nsite` | Sum of mutation rates for all possible missense SNVs in the gene |
| `ssite`     | Sum of mutation rates for all possible synonymous SNVs in the gene |