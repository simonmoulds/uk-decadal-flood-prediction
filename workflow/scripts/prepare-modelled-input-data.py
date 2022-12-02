#!/usr/bin/env python3

import click

from utils import _ensemble_preprocessor

@click.command()
@click.option('-i', '--inputdir', default='.', help='Input directory')
@click.option('-o', '--outputdir', default='.', help='Output directory')
@click.option('--config', default='config.yml', help='YAML configuration file')
def main(inputdir, outputdir, config):
    _ensemble_preprocessor(inputdir, outputdir, config)

if __name__ == '__main__':
    main()
