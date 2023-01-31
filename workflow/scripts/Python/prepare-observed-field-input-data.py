#!/usr/bin/env python3

from utils import observed_field_preprocessor

observed_field_preprocessor(
    snakemake.config['input_data_root'],
    snakemake.params['outputdir'],
    snakemake.input['stations'],
    snakemake.input['grid'],
    snakemake.config
)
