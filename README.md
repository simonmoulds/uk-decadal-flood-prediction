# Snakemake workflow: `uk-decadal-flood-prediction`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/simonmoulds/uk-decadal-flood-prediction/workflows/Tests/badge.svg?branch=main)](https://github.com/simonmoulds/uk-decadal-flood-prediction/actions?query=branch%3Amain+workflow%3ATests)

A Snakemake workflow for making decadal flood predictions using CMIP5/6 decadal hindcasts. The workflow reproduces the results shown in the research article: 

```
Moulds S, Slater LJ, Dunstone NJ, Smith DM (2023). Skillful decadal flood prediction. Geophysical Research Letters, 49, e2022GL100650. https://doi.org/10.1029/2022GL100650
```

## Usage

This workflow uses Snakemake. Users should install Snakemake using [Mamba/Conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba). It can then be run from the command line as follows: 

``` sh
snakemake --cores N --use-conda
```
where N is the number of cores you wish to use. Note that the workflow has been tested on MacOS (Monterey 12.6) and Linux systems only. On Windows it should work under the Windows Subsystem for Linux. 

The raw data needed to run the workflow is not available on this repository. It can be downloaded from the corresponding [Zenodo dataset](https://doi.org/10.5281/zenodo.6940449).

## Authors

This workflow was developed by:

- [Simon Moulds](https://github.com/simonmoulds)

## Acknowledgements

This work is funded by UK Research and Innovation grant MR/V022008/1 (SM and LJS). NJD and DMS were supported by the Met Office Hadley Centre Climate Programme funded by BEIS and Defra. The authors would like to acknowledge the use of the University of Oxford Advanced Research Computing (ARC) facility in carrying out this work (http://dx.doi.org/10.5281/zenodo.22558).

The workflow itself is based on the following research article:
```
Klaus B and Reisenauer S. An end to end workflow for differential gene expression using Affymetrix microarrays [version 2; peer review: 2 approved]. F1000Research 2018, 5:1384 (https://doi.org/10.12688/f1000research.8967.2)
```

## License

This workflow is licensed under the [MIT](LICENSE.md) license. 
