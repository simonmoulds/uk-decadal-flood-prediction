#!/usr/bin/env python3

from utils import ensemble_field_preprocessor

ensemble_field_preprocessor(
    snakemake.config['input_data_root'],
    snakemake.params['outputdir'],
    snakemake.input['stations'],
    snakemake.input['grid'],
    snakemake.config
)
