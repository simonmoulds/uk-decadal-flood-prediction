#!/usr/bin/env python3

from utils import ensemble_preprocessor

ensemble_preprocessor(
    snakemake.config['input_data_root'],
    snakemake.params['outputdir'],
    snakemake.config
)
