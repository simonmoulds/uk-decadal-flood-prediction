#!/usr/bin/env python3

from utils import observed_preprocessor

observed_preprocessor(
    snakemake.config['input_data_root'],
    snakemake.params['outputdir'],
    snakemake.config
)
