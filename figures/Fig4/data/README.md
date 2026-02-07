# Overview

If you have any questions, please feel free to contact me via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

- [`ASD-related_gene.lst`](#asd-related_genelst)
- [`DD_denovo.RDH_distribution.tsv`](#dd_denovordh_distributiontsv)
- [`denovoWEST_DD-related_gene.lst`](#denovowest_dd-related_genelst)
- [`Fu2022_Nature_autism_denovo_snv-indel.ST7.per-gene_case-control-count.csv`](#fu2022_nature_autism_denovo_snv-indelst7per-gene_case-control-countcsv)
- [`gencodev44_DD_RGC_integration.tsv`](#gencodev44_dd_rgc_integrationtsv)
- [`ukb.RDH_burden.WBnR.association.continu.MSA-MTR_0.649_AF0.001.tsv`](#ukbrdh_burdenwbnrassociationcontinumsa-mtr_0649_af0001tsv)

## `ASD-related_gene.lst`

This file lists autism spectrum disorder (ASD)–associated genes identified using the TADA framework, as reported by [Fu et al.](https://doi.org/10.1038/s41588-022-01104-0).

## `DD_denovo.RDH_distribution.tsv`

This file summarizes the distribution of de novo mutations (DNMs) from developmental disorder (DD) patients in core RDH genes, together with relevant annotations from the DECIPHER database. This table corresponds to **Supplementary Table 5** in the paper.

**Data sources:**

- RGC-ME: <https://rgc-research.regeneron.com/me/home>
- DNMs in DD patients: <https://github.com/HurlesGroupSanger/DeNovoWEST>
- DECIPHER: <https://www.deciphergenomics.org>

## `denovoWEST_DD-related_gene.lst`

This file lists DD-associated genes identified using the DeNovoWEST framework, as reported by [Kaplanis et al.](https://doi.org/10.1038/s41586-020-2832-5).

## `Fu2022_Nature_autism_denovo_snv-indel.ST7.per-gene_case-control-count.csv`

This file was obtained from **Supplementary Table 7** of [Fu et al.](https://doi.org/10.1038/s41588-022-01104-0).

The definition of each column can be found in the original supplementary table.

## `gencodev44_DD_RGC_integration.tsv`

This table summarizes **gene-level variant counts** across gnomAD, RGC-ME, and DNMs from DD patients.

### Columns

| Column                     | Description                                           |
| -------------------------- | ----------------------------------------------------- |
| `gene_name`                | Gene name                                             |
| `[dataset]_[consequence]`  | Variant count for each dataset–consequence combination |

### How to generate

Data were integrated from the following sources:

- Public population datasets (gnomAD, RGC-ME)
- DNMs from DD patients (as described above)
- Functional annotation performed using **ANNOVAR**

## `ukb.RDH_burden.WBnR.association.continu.MSA-MTR_0.649_AF0.001.tsv`

This table summarizes the association between cognitive phenotypes and coding variant burden in RDH genes. It corresponds to **Supplementary Table 7** in the paper.