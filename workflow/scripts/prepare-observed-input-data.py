#!/usr/bin/env python3

import click

from utils import _observed_preprocessor

@click.command()
@click.option('-i', '--inputdir', default='.', help='Input YAML configuration file')
@click.option('-o', '--outputdir', default='.', help='Output directory')
def main(inputdir, outputdir):
    _observed_preprocessor(inputdir, outputdir)

if __name__ == '__main__':
    main()
