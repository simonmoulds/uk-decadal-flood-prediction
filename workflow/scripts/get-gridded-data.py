# #!/usr/bin/env python3

# import pandas as pd
# import yaml
# import click

# from utils import _ensemble_preprocessor

# config = 'config/config.yml'
# with open(config, 'r') as f:
#     config = yaml.load(f, Loader=yaml.FullLoader)

# best_n = 20

# dirpath = config['ensemble_data']['cmip5']['grid_subdirectory']
# matched_forecast_error = pd.read_parquet('results/analysis/yr2to9_lag/matched_ensemble_error.parquet')
# matched_forecast_error = matched_forecast_error.set_index(['project', 'source_id', 'mip', 'member', 'init_year_matched'])
# matched_forecast_spec = matched_forecast_error.groupby('init_year')['error'].nsmallest(best_n).reset_index()
# yrs = matched_forecast_spec['init_year'].unique()
# for i in range(len(yrs)):
#     yr = yrs[i]
#     current_spec = matched_forecast_spec.loc[matched_forecast_spec['init_year'] == yr]
#     filepath = 'TODO'
#     # Read file

# @click.command()
# @click.option('-i', '--input', 'config', default='config.yml', help='Input YAML configuration file')
# @click.option('-o', '--outputdir', default='.', help='Output directory')
# def main(config, outputdir):
#     _ensemble_preprocessor(config, outputdir)

# if __name__ == '__main__':
#     main()
