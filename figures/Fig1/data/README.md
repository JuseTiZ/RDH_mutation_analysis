# Overview

If you have any questions, please feel free to contact us via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.

  - [`Human.gencodev44.histone.tsv`](#humangencodev44histonetsv)
  - [`human.gencodeV44.pc_canonical.mut_info.tsv`](#humangencodev44pc_canonicalmut_infotsv)
  - [`human.ortholog_mouse_monkey.dnds_similarity.tsv`](#humanortholog_mouse_monkeydnds_similaritytsv)
  - [`multispecies_4dsite.1_mer_record.tsv`](#multispecies_4dsite1_mer_recordtsv)


## `Human.gencodev44.histone.tsv`

This file lists histone genes from the human GENCODE v44 annotation. Each record includes gene ID, gene name and histone type. Entries were compiled with reference to [Seal et al.](https://doi.org/10.1186/s13072-022-00467-2).

## `human.gencodeV44.pc_canonical.mut_info.tsv`

This table reports rare-variant densities for each human gene across various genomic regions (CDS, exon, transcript, etc.).

| Column                             | Description                                                              |
| ---------------------------------- | ------------------------------------------------------------------------ |
| `gene id`, `gene name`, `transcript id`, `gene type`                             | Gene identifiers and basic annotation                                 |
| `region` | Genomic region used to calculate rare-variant density |
| `bigwig name`  | Name of the bigWig track used to compute densities                       |
| `mean value`    | Calculated rare-variant density                                          |
| `non-nan length`, `nan percentage`    | Effective length after coverage/mappability filtering and the fraction masked as NaN |
| `chrom`, `start`    | Genomic coordinate of the annotated transcription start site           |

### How to generate

1. Create a bedGraph from gnomAD variant data (put the rare-variant count in the 4th column).
2. Apply coverage and mappability filters (for `4dsite_snv`, additionally filter 4-fold degenerate sites), then convert the result to a bigWig track.
3. Produce this table with `scripts/calcu_transcript_feature.py`.

### Data sources

1. gnomAD variant database: [https://gnomad.broadinstitute.org/data#v4](https://gnomad.broadinstitute.org/data#v4).
2. human GENCODE V44 annotation: [https://www.gencodegenes.org/human/release_44.html](https://www.gencodegenes.org/human/release_44.html).

## `human.ortholog_mouse_monkey.dnds_similarity.tsv`

This file provides inter-species metrics comparing human to mouse and human to monkey orthologs.

### How to generate

1. Identify orthologs using [OrthoFinder](https://github.com/davidemms/OrthoFinder) and standardized histone gene names.
2. Compute dN and dS with [kaks_calculator](https://ngdc.cncb.ac.cn/biocode/tools/BT000001). For sequence similarity we directly compare aligned peptide pairs:
```
def calculate_similarity(pep1, pep2):
    assert len(pep1) == len(pep2), "Peptide lengths do not match."
    matches = sum(1 for a, b in zip(pep1, pep2) if a == b)
    similarity = matches / len(pep1)
    return similarity
```


## `multispecies_4dsite.1_mer_record.tsv`

This file summarizes variant counts on 4-fold degenerate sites for RDH genes and other genes across multiple species.

| Column                             | Description                                                              |
| ---------------------------------- | ------------------------------------------------------------------------ |
| `gene group`     | Gene group (two possible values: `RDH` or `other_gene`)                 |
| `species` | Species name                                                            |
| `mutation`  | 1-mer mutation type (single base change)                                |
| `total num`    | Count of rare variants                                                  |
| `total site`    | Number of considered sites                                               |
| `total mask`    | Number of masked sites (used to compute corrected rare-variant density)  |

Sources for the data used to compute these statistics are referenced in the paper’s supplementary tables.
