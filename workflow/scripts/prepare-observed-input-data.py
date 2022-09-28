#!/usr/bin/env python3

import click

from utils import _observed_preprocessor

@click.command()
@click.option('-i', '--input', 'config', default='config.yml', help='Input YAML configuration file')
@click.option('-o', '--outputdir', default='.', help='Output directory')
def main(config, outputdir):
    _observed_preprocessor(config, outputdir)

if __name__ == '__main__':
    main()
