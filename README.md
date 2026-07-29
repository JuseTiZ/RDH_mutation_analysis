# RDH_mutation_analysis

- [1. Introduction](#1-introduction)
- [2. Repository Structure](#2-repository-structure)
- [3. Environment](#3-environment)
- [4. Citation](#4-citation)
- [5. License](#5-license)
- [6. Contact](#6-contact)

## 1. Introduction

This repository contains the main datasets and scripts used in the analyses described in  **"Hypermutability of ultraconserved histone genes and its contribution to human disease"**.

## 2. Repository Structure

This repository is organized into two main directories:

- `figures/`  
  Contains the data and scripts used to generate the main figures presented in the paper.

- `scripts/`  
  Contains the primary scripts used to generate specific datasets analyzed in this study.

Detailed descriptions of the data files, scripts, and directory organization can be found in the `README.md` files within each corresponding folder.

## 3. Environment

Create and activate the tested Conda environment with:

```bash
conda env create -f environment.yml
conda activate rdhgs-mutation-analysis
```

The environment contains the dependencies used by all Python scripts and figure notebooks in this repository. `figures/Fig2/fig2_c_d.ipynb` additionally refers to external BigWig, GTF, and deepTools matrix files through machine-specific paths; replace those paths when running the notebook on another system.

The MuSiCal dependency used by `figures/Fig2/fig2_e_f.ipynb` is distributed under its own academic and non-commercial research license. The MIT license of this repository does not replace the licenses of third-party dependencies or datasets.

## 4. Citation

If you find this repository, data, or scripts helpful in your research, please cite:

```

xxx

```

## 5. License

This repository is available under the [MIT License](LICENSE). Third-party software and datasets retain their respective licenses and terms of use.

## 6. Contact

If you encounter any issues with this repository or have questions beyond the scope of the provided documentation, please feel free to contact us via **[e-mail](mailto:jiangzj6@mail2.sysu.edu.cn)** or open an **[issue on GitHub](https://github.com/JuseTiZ/RDHGs_mutation_analysis/issues)**.
