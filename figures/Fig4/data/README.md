# Overview

If you have any questions, please feel free to contact us via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

- [`ASD-related_gene.lst`](#asd-related_genelst)
- [`DD_denovo.RDH_distribution.tsv`](#dd_denovordh_distributiontsv)
- [`denovoWEST_DD-related_gene.lst`](#denovowest_dd-related_genelst)
- [`Fu2022_Nature_autism_denovo_snv-indel.ST7.per-gene_case-control-count.csv`](#fu2022_nature_autism_denovo_snv-indelst7per-gene_case-control-countcsv)
- [`gencodev44_DD_gnomADnfe_integration.tsv`](#gencodev44_dd_gnomadnfe_integrationtsv)
- [`human_cRDH_gene.gencodev44.lst`](#human_crdh_genegencodev44lst)
- [`RDH_all_possible_mis_syn_mutation_annotated.tsv`](#rdh_all_possible_mis_syn_mutation_annotatedtsv)
- [`ukb.RDH_burden.WBnR.association.continu.HistMTR_0.630_AF0.001.tsv`](#ukbrdh_burdenwbnrassociationcontinuhistmtr_0630_af0001tsv)


## `ASD-related_gene.lst`

This file lists autism spectrum disorder (ASD)–associated genes identified using the TADA framework, as reported by [Fu et al.](https://doi.org/10.1038/s41588-022-01104-0).

## `DD_denovo.RDH_distribution.tsv`

This file summarizes the distribution of de novo mutations (DNMs) from developmental disorder (DD) patients in core RDH genes, together with relevant annotations from the DECIPHER database. This table corresponds to **Supplementary Table 8** in the paper.

**Data sources:**

- RGC-ME: <https://rgc-research.regeneron.com/me/home>
- DNMs in DD patients: <https://github.com/HurlesGroupSanger/DeNovoWEST>
- DECIPHER: <https://www.deciphergenomics.org>

## `denovoWEST_DD-related_gene.lst`

This file lists DD-associated genes identified using the DeNovoWEST framework, as reported by [Kaplanis et al.](https://doi.org/10.1038/s41586-020-2832-5).

## `Fu2022_Nature_autism_denovo_snv-indel.ST7.per-gene_case-control-count.csv`

This file was obtained from ***Supplementary Table 7*** of [Fu et al.](https://doi.org/10.1038/s41588-022-01104-0).

The definition of each column can be found in the original supplementary table.

## `gencodev44_DD_gnomADnfe_integration.tsv`

This table summarizes **gene-level variant counts** based on gnomAD SNPs (NFE subset) and DNMs from DD patients.

### Columns

| Column                     | Description                                           |
| -------------------------- | ----------------------------------------------------- |
| `gene_name`                | Gene name                                             |
| `[dataset]_[consequence]`  | Variant count for each dataset–consequence combination |

### How to generate

Data were integrated from the following sources:

- gnomAD variant database: [https://gnomad.broadinstitute.org/data#v4](https://gnomad.broadinstitute.org/data#v4).
- DNMs from DD patients (as described above)
- Functional annotation performed using **ANNOVAR**

## `RDH_all_possible_mis_syn_mutation_annotated.tsv`

This table provides detailed annotations for all possible missense and synonymous SNVs in core RDH genes.

### Columns

| Column | Description |
|--------|-------------|
| `chr`, `pos`, `ref`, `alt` | Variant information |
| `gene name`, `consequence`, `amino-acid change`, `gene family` | Functional annotation of the variant |
| `consensus amino-acid`, `consensus idx` | Consensus position in the family alignment |
| `no-gap consensus idx` | Consensus index after removing alignment gaps |
| `mutation rate` | MuRaL-predicted mutation rate for this variant |
| `HistMTR` | HistMTR value for this amino-acid position |
| `DD DNM count` | Occurrence count in DNMs from DD patients |
| `is in gnomAD NFE` | Whether the variant is present in the gnomAD NFE subset |

## `ukb.RDH_burden.WBnR.association.continu.HistMTR_0.630_AF0.001.tsv`

This table summarizes the association between cognitive phenotypes and coding variant burden in core RDH genes. It corresponds to **Supplementary Table 9** in the paper.