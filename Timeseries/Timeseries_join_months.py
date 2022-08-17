import os
import sys

from han_func import join_monthly_files

year = sys.argv[1]
site = sys.argv[2]
sector = sys.argv[3]
climatology = True if sys.argv[4].lower()=='true' else False

# change to directory containing fpdm files
ts_path = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2')

clim_str = f'_climatology' if climatology else ''

print(f'Searching for files with pattern: {site}_{sector}_timeseries{clim_str}')

join_monthly_files(file_pattern = f'{site}_{sector}_timeseries{clim_str}',
                   year = year,
                   path = ts_path,
                   zip_months = False,
                   remove_months = True,
                   verbose = True)