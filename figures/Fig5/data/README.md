# Overview

If you have any questions, please feel free to contact us via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

- [`stemloop_hde_genomic3UTR.density_phylop_comp.tsv`](#stemloop_hde_genomic3utrdensity_phylop_comptsv)
- [`ukb.stemloop_burden.WBnR.aggregate_mutation_association.binary.tsv`](#ukbstemloop_burdenwbnraggregate_mutation_associationbinarytsv)

## `stemloop_hde_genomic3UTR.density_phylop_comp.tsv`

This table contains rare variant densities and PhyloP 100-way conservation scores for regulatory elements. For stem-loop structures and HDEs, values were calculated at the gene level for each RDH gene. For genomic 3′ UTRs, values were calculated for each continuous 3′ UTR fragment across the human genome.

### Data required to generate this file

1. A BigWig track of rare variant counts, as described in `figures/Fig1/data/README.md`.
2. PhyloP 100-way conservation scores, available at:  
   <https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phyloP100way/hg38.phyloP100way.bw>
3. A 3′ UTR BED file generated using [gencode_regions](https://github.com/saketkc/gencode_regions) and merged with [bedtools](https://github.com/arq5x/bedtools2).

## `ukb.stemloop_burden.WBnR.aggregate_mutation_association.binary.tsv`

This table summarizes the association between reproductive phenotypes and regulatory variant burden in RDH genes. It corresponds to **Supplementary Table 11** in the paper.