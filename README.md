# Fine-mapping with infinitesimal effects v1.2
This repo contains code for SuSiE-inf (v1.4) and FINEMAP-inf (v1.3)

## Getting started
Python3.9 is recommended

Clone this repo:
```
git clone https://github.com/FinucaneLab/fine-mapping-inf.git
cd fine-mapping-inf
```
Install dependencies:
```
setuptools
wheel
numpy
pandas
scipy
bgzip
```

We recommend that you install both SuSiE-inf and FINEMAP-inf.

To install SuSiE-inf:
```
cd susieinf
python setup.py bdist_wheel
pip install .
```
Similarly, to install FINEMAP-inf:
```
cd finemapinf
python setup.py bdist_wheel
pip install .
```
Once completed, run
```
python run_fine_mapping.py -h
```
to print a list of all command-line options.

## Tutorial
Please see [wiki](https://github.com/FinucaneLab/fine-mapping-inf/wiki) for tutorial.
**Note that both SuSiE-inf and FINEMAP-inf are developed for use in single-cohort fine-mapping with in-sample LD, for fine-mapping of meta-analyzed GWAS with reference panel LD, please see methods like [SLALOM](https://www.cell.com/cell-genomics/pdf/S2666-979X(22)00163-X.pdf) or [CARMA](https://www.nature.com/articles/s41588-023-01392-0).**

## Authors
Ran Cui
Zhou Fan

## Contact
For questions and comments, please contact Ran Cui at rancui@broadinstitute.org
