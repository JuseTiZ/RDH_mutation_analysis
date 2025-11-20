# Overview

If you have any questions, please feel free to contact me via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

  - [`Bagert_histone_mut.benchmark.tsv`](#bagert_histone_mutbenchmarktsv)
  - [`cRDH.MSA_MTR.obs_RGC.pred_MuRaL.csv`](#crdhmsa_mtrobs_rgcpred_muralcsv)
  - [`cRDH.PTM.intra_interaction.summary.pkl`](#crdhptmintra_interactionsummarypkl)
  - [`cRDH.secondary_structure_anno.json`](#crdhsecondary_structure_annojson)
  - [`gene.variant_record.gnomad_RGCME.csv`](#genevariant_recordgnomad_rgcmecsv)

## `Bagert_histone_mut.benchmark.tsv`

This file corresponds to **Supplementary Table 2**. It serves as a benchmark dataset linking missense substitutions in cRDHs to experimentally measured phenotypic effects.

### How this file was generated

A reproducible workflow follows the methods described in the original paper. Briefly:

1. **Enumerate all possible SNVs** across cRDHG coding sequences to produce a comprehensive VCF.
2. **Annotate functional consequence** using **ANNOVAR**.
3. For each missense SNV, **compute functional scores** (REVEL, AlphaMissense, MSA-MTR, MPC, ProSST, ESM-1v, etc.) using either VEP plugins or in-house Python parsing.
4. **Retrieve quantitative phenotypic effects** from
   *Bagert et al., Nat. Chem. Biol., 2021*
   DOI: [https://doi.org/10.1038/s41589-021-00738-1](https://doi.org/10.1038/s41589-021-00738-1)
5. For each amino acid substitution, **average metric values across all SNVs** that can generate that amino-acid change (unnecessary for some metrics).

### Resources for functional metrics

* **MSA-MTR**

  * Population-observed variants: **RGC-ME**
    [https://rgc-research.regeneron.com/me/home](https://rgc-research.regeneron.com/me/home)
  * Predicted mutation rate: **MuRaL**
    [https://github.com/CaiLiLab/MuRaL](https://github.com/CaiLiLab/MuRaL)

* **REVEL & AlphaMissense (via VEP)**
  Example VEP command:

  ```bash
  singularity exec vep.sif \
      vep --dir [dir] \
          --cache --cache_version 110 --offline --format vcf --vcf --force_overwrite --assembly GRCh38 \
          --input_file [input.vcf] \
          --output_file [output.vcf] \
          --compress_output bgzip \
          --plugin AlphaMissense,file=[path]/AlphaMissense_hg38.tsv.gz \
          --plugin REVEL,file=[path]/new_tabbed_revel_grch38.tsv.gz,no_match=1
  ```

* **phyloP100way conservation**
  [https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/](https://hgdownload.cse.ucsc.edu/goldenpath/hg38/phyloP100way/)

* **RGC-MTR**
  [https://doi.org/10.6084/m9.figshare.24587328](https://doi.org/10.6084/m9.figshare.24587328)

* **MPC**

  ```bash
  gsutil -m cp -r \
    "gs://gcp-public-data--gnomad/release/2.1.1/regional_missense_constraint/gnomad_v2.1.1_mpc_liftover_GRCh38.ht" \
    .
  ```

* **ProSST**
  [https://github.com/ai4protein/ProSST](https://github.com/ai4protein/ProSST)

* **ESM-1v functional prediction**
  Pipeline from: [https://github.com/ntranoslab/esm-variants](https://github.com/ntranoslab/esm-variants)

## `cRDH.MSA_MTR.obs_RGC.pred_MuRaL.csv`

This file reports **MSA-MTR** values for every position across cRDHs. It corresponds to **Supplementary Table 4**.

Values are derived from:

* **Observed variants** from RGC-ME
* **Predicted mutation rates** from MuRaL

A detailed description of the generation workflow can be found in the methods section of the manuscript.

## `cRDH.PTM.intra_interaction.summary.pkl`

This file contains **post-translational modification (PTM)** annotations and **intra-nucleosome interaction** profiles for each cRDH family.

### File structure

Loading it in Python returns a dictionary:

```python
{
  "H2A": [{...}, {...}],
  "H2B": [{...}, {...}],
  "H3": [{...}, {...}],
  "H4": [{...}, {...}],
}
```

Values include:

* PTM types at each residue
* Interaction types

### Color mapping used for PTMs & interactions

#### PTMs

```python
{
 'Phosphorylation': '#A74EB6',
 'Acetylation': '#2E7D2A',
 'Methylation': '#8CD436',
 'Ubiquitination': '#BB4E17',
 'Sumoylation': '#F29D33',
 'Glycosylation': '#3380C0',
 'None': 'azure'
}
```

#### Interaction categories

```python
{
  "None": "#FFF0ED",
  "vdW": "#EF9575",
  "H-bonds (HB)": "#88C5F8",
  "Salt bridges (SB)": "gold",
  "SB & HB": "#BE8FF4"
}
```

### Data sources

* **PTM information:**
  Retrieved from **PTMcosmos**
  [https://ptmcosmos.wustl.edu/](https://ptmcosmos.wustl.edu/)

* **Residue–residue contacts:**
  Computed using **PDBe Arpeggio**
  [https://github.com/PDBeurope/arpeggio](https://github.com/PDBeurope/arpeggio)

* **3D structural templates:**
  AlphaFold Server
  [https://alphafoldserver.com/](https://alphafoldserver.com/)

## `cRDH.secondary_structure_anno.json`

This JSON file provides **secondary structure annotations** for cRDH.

* Keys correspond to histone families (`H2A`, `H2B`, `H3`, `H4`).
* These annotations are used for visualization in structural plots.

### Reproducing these annotations

The secondary structure was parsed from the **7PFV nucleosome structure**. You can regenerate the annotation using `get_structure()` in `scripts/plot_structure.py`. The function identifies helix/sheet regions and returns their residue intervals.

The positions of the patches and tail domains were referenced from [Nacev et al.](https://doi.org/10.1038/s41586-019-1038-1) and [Armeev et al.](https://doi.org/10.1038/s41467-021-22636-9).

## `gene.variant_record.gnomad_RGCME.csv`

This table summarizes **gene-level variant count** across gnomAD and RGC-ME, stratified by functional consequence and AF bins.

### Columns include:

| Column                             | Description                                                              |
| ---------------------------------- | ------------------------------------------------------------------------ |
| `gene`                             | Gene name                                                                |
| `[dataset]_[consequence]_[AF bin]` | Count of variants in each combination |
| `non_nsite`, `mis_nsite`, `ssite`  | Predicted mutation-rate sums from MuRaL                                  |
| `is_cRDHGs`                        | Boolean indicating whether the gene belongs to the cRDH family          |

### How to generate

Data originate from:

* Public population datasets
  (gnomAD, RGC-ME)
* Functional annotation via **ANNOVAR**
* Mutation-rate prediction via **MuRaL**