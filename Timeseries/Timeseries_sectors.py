import os
import sys
import numpy as np
import xarray as xr

from acrg.convert import concentration
from acrg.name.name import footprints, flux, flux_for_HiTRes, timeseries_HiTRes

################# Input information #################
year = int(sys.argv[1])
month = None if str(sys.argv[2]).lower()=='none' else int(sys.argv[2])
site = str(sys.argv[3])
sector = str(sys.argv[4])

species = ['o2', 'co2']

print(f'Sector: {sector}')

HiTRes = True

start_month = 1 if month is None else month
end_month = 1 if month is None or month==12 else month+1
end_year = year+1 if month is None or month==12 else year

start = f'{year}-{str(start_month).zfill(2)}-01'
end = f'{end_year}-{str(end_month).zfill(2)}-01'

################# Get footprints #################
fp_ds = footprints(site, 'UKV',
                   start = start,
                   end = end, 
                   domain = 'EUROPE', 
                   species = 'co2',
                   HiTRes = True,
                   chunks = {'time': 50})

################# Get emissions #################
emissions_names = {spec: {'high_freq': f'{spec}-ukghg-ff-{sector}-1hr'} for spec in species}

print('\nemissions names:')
[print(f'{spec}: {emissions_name}') for spec, emissions_name in emissions_names.items()]

emissions = flux_for_HiTRes(domain = 'EUROPE', 
                            emissions_dict = emissions_names, 
                            start = start, 
                            end = end,
                            chunks = {'time': 50})

# make sure the lat & lon of the footprints and emissions are formatted the same
for spec, em_spec in emissions.items():
    for freq, em_freq in em_spec.items():
        for ll in ['lat', 'lon']:
            em_freq[ll] = fp_ds[ll]

################# Calculate timeseries #################
mf_ts = timeseries_HiTRes(flux_dict = emissions,
                          fp_HiTRes_ds = fp_ds, 
                          output_TS = True,
                          output_fpXflux = False)
for dv in mf_ts.data_vars:
    mf_ts[dv] = mf_ts[dv] / concentration('ppm')
    mf_ts[dv] = mf_ts[dv].assign_attrs({'units': '1e-6'})
if 'total' in mf_ts.data_vars:
    mf_ts = mf_ts.rename({'total': list(emissions.keys())[0]})

################# Save results #################
date_str = f'{year}' if month is None else \
           f'{year}{str(month).zfill(2)}'
filename = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2',
                        f'{site}_ff-ukghg-{sector}_timeseries_{date_str}.nc')

# save to netcdf
print(f'\n\n\nSaving to: {filename}')
mf_ts.to_netcdf(filename)
print('\nDone')