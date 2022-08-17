import os

from acrg.name import name

year = 2018
species = 'o2'
or_ff_model = 'gridfed-ukghg'

species_names = {'co2': 'co2-edgar-ukghg-ff-1hr',
                 f'ox_ratio': f'oxidativeratio-{or_ff_model}-1hr'}
species_names = {'co2': species_names['co2']} if species=='co2' else species_names
method = 'constant_flux' #'integrated_fp' #

# get the flux data
flux = {gas: name.flux(domain = 'EUROPE',
                       species = species_name,
                       start = f'{year}-01-01',
                       end = f'{year+1}-01-01',
                       chunks = {'time': 50})
        for gas, species_name in species_names.items()}

if species=='o2':
    # fix differences in floating points
    for ll in ['lat', 'lon']:
        flux['ox_ratio'][ll] = flux['co2'][ll]
# estimate o2 uptake
flux = flux['ox_ratio'] * flux['co2'] if species=='o2' else flux['co2'] 

# average the monthly fluxes and then reindex to hourly
if method == 'constant_flux':
    flux = flux.resample(time='1MS').mean(dim='time').reindex_like(flux, method='ffill')

# get the footprints
fp = name.footprints('WAO', met_model='UKV',
                     start = f'{year}-01-01',
                     end = f'{year+1}-01-01',
                     domain = 'EUROPE',
                     chunks = {'time': 50})
# multiply the footprints by the fluxes to estimate the timeseries
ts = (flux.flux * fp.fp).sum(["lat", "lon"])
ts = ts.compute()

# convert to a dataset
ts_ds = ts.to_dataset(name='flux')

# create a filename
or_ff_model = f'_{or_ff_model}' if species=='o2' else ''
file_name = f'{species}_WAO_ff{or_ff_model}_timeseries_constflux_{year}.nc' if method=='constant_flux' else \
            f'{species}_WAO_ff{or_ff_model}_timeseries_intfp_{year}.nc' if method=='integrated_fp' else None

# ave to netcdf
if file_name is not None:
    print(f'Saving to {file_name}') 
    ts_ds.to_netcdf(os.path.join('/user/work/vf20487/Timeseries/o2_co2', file_name))
else:
    print('Cannot find file name')

