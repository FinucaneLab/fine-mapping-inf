# Fine-mapping with infinitesimal effects v1.0
This repo contains code for SuSiE-inf (v1.2) and FINEMAP-inf (v1.3)

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

## Author
Zhou Fan (Yale)

## Contact
For questions and comments, please contact Ran Cui at rancui@broadinstitute.org
