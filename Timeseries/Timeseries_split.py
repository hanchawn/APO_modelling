'''
Run the CO2 forward model for one APO component

The inputs are imported from the bash script:

- year,
- month,
- site,
- source - e.g. ocean, bio, ff,
- biospheric model - e.g. orchidee, lpjguess, cardamom,
- ff model to use - e.g. edgar, ukghg
- oxidative ratio model to use - e.g. gridfed, ukghg,
- climatology - True or False

Outputs a netcdf file with the timeseries for all APO species associated the source specified

'''

import os
import sys
import numpy as np
import xarray as xr
from datetime import datetime
from dateutil.relativedelta import relativedelta

from acrg.convert import concentration
from acrg.name.name import footprints, flux, flux_for_HiTRes, timeseries_HiTRes

import apo_funcs

################# Input information #################
year = int(sys.argv[1])
month = None if str(sys.argv[2]).lower()=='none' else int(sys.argv[2])
site = str(sys.argv[3])
source = str(sys.argv[4])
bio_model = str(sys.argv[5])
ff_model = str(sys.argv[6])
oxidativeratio = None if str(sys.argv[7]).lower()=='none' else str(sys.argv[7])
climatology = True if str(sys.argv[8]).lower()=='true' else False

HiTRes = True if source!='ocean' else False
fp_species = 'co2' if HiTRes else None

# get the start and end date from the year and month above
time_diff = relativedelta(years=1, hours=-1) if month is None else relativedelta(months=1)
mth = 1 if month is None else month
start = datetime(year, mth, 1)
end = start + time_diff

print(f'Calculating timeseries between {start} and {end}')

################# Get footprints #################
height = '185magl' if site=='TAC' else None
fp_ds = footprints(site, 'UKV',
                   start = start,
                   end = end, 
                   height = height,
                   domain = 'EUROPE', 
                   species = fp_species,
                   HiTRes = HiTRes,
                   chunks = {'time': 50})

################# Get emissions #################
if 'ff' in source:
    ff_spec = ['co2', 'o2'] if oxidativeratio is None else ['co2']
    emissions_names = {spec: {'high_freq': f'{spec}-{ff_model}-ff-1hr'} for spec in ff_spec}
    if oxidativeratio:
        emissions_names['oxidativeratio'] = {'high_freq': f'oxidativeratio-{oxidativeratio}-1hr'}
elif source=='bio' and bio_model=='cardamom':
    emissions_names = {f'co2-{bio_type}': {'high_freq': f'co2-{bio_type}-cardamom-2hr'}
                       for bio_type in ['gpp', 'rtot']}
elif source=='bio' and bio_model=='orchidee':
    emissions_names = {f'co2-{bio_type}': {'high_freq': f'co2-{bio_type}-orchidee-3hr'}
                       for bio_type in ['gpp', 'rtot']}
else:
    # the species available for each ocean flux for each resolution
    species_sims = {'day': {'ecco': ['co2', 'o2'], 'jena': ['co2', 'o2'], 'nemo': ['co2', 'o2', 'n2']},
                    'mth': {'ecco': ['co2', 'o2'], 'nemo': ['co2', 'o2', 'n2']}}
    # determine whether a climatology is needed
    clim_str = {sim: '-climatology' if year>year or climatology else ''
                for sim, year in {'ecco': 2018, 'jena': 2020, 'nemo': 2015}.items()}

    # get all of the needed emissions names
    emissions_names = {res: {sim: {f'{spec}_{sim}_{res}': f'{spec}{clim_str[sim]}-{sim}-ocean-{res}'
                                   for spec in specs_sims}
                             for sim, specs_sims in specs_res.items()}
                       for res, specs_res in species_sims.items()}
    # flatten the dictionary
    emissions_names = apo_funcs.flatten_nested_dict(emissions_names, join_keys=False)
    emissions_names = {'_'.join(spec.split('_')[1:]): em_name for spec, em_name in emissions_names.items()}

print('\nemissions names:')
[print(f'{spec}: {emissions_name}') for spec, emissions_name in emissions_names.items()]

# cardamom fluxes go up to 2015 and orchidee up to 2018, if the year is after that use the last year of data
if any([all([source=='bio' and year>2015 and bio_model=='cardamom']),
        all([source=='bio' and year>2018 and bio_model=='orchidee'])]):
    year_flux = 2015 if bio_model=='cardamom' else 2018

    print(f'Finding {bio_model} {source} emissions data for {year_flux}')
    start_flux = datetime(year_flux, mth, 1)
    end_flux = start_flux + time_diff

    # time will be replaced to match footprint year
    change_time = True
    
else:
    start_flux = start
    end_flux = end
    change_time = False

# get the high time res fluxes
print(f'Getting flux data between {start_flux} and {end_flux}')
if source!='ocean':
    emissions = flux_for_HiTRes(domain = 'EUROPE', 
                                emissions_dict = emissions_names, 
                                start = start_flux, 
                                end = end_flux,
                                chunks = {'time': 50})

    if source=='ff' and oxidativeratio is not None:
        for ll in ['lat', 'lon']:
            emissions['oxidativeratio']['high_freq'][ll] = emissions['co2']['high_freq'][ll]
        emissions['o2'] = {'high_freq': emissions['oxidativeratio']['high_freq'] * emissions['co2']['high_freq']}
        emissions.pop('oxidativeratio', None)
        
else:
    print('Getting ocean flux data')
    emissions = {}
    for spec, em_spec in emissions_names.items():
        try:
            emissions[spec] = flux(domain = 'EUROPE', 
                                   species = em_spec, 
                                   start = start_flux, 
                                   end = end_flux)
        except OSError:
            print(f'Cannot find file for {em_spec}')

################# Adjust emissions #################
# if the year is >2015 and source is bio then replace the time to match the footprints
if change_time:
    for spec, em_spec in emissions.items():
        for freq, em_freq in em_spec.items():
            start_year = year if source=='ocean' or month!=1 else year-1
            new_start = datetime.strptime(em_freq.time.values[0].astype(str).split('T')[0], '%Y-%m-%d')
            new_start = f'{start_year}-{str(new_start.month).zfill(2)}-{str(new_start.day).zfill(2)}'
            
            # time step of the flux data
            time_step = int(em_freq.time.diff(dim="time").values.mean() / np.timedelta64(1, 'h'))
            print(f'Time step: {time_step}h')
            # new time array
            time = np.arange(new_start, end, time_step, dtype='datetime64[h]').astype('datetime64[ns]')

# make sure the lat & lon of the footprints and emissions are formatted the same
for spec, em_spec in emissions.items():
    for ll in ['lat', 'lon']:
        em_spec[ll] = fp_ds[ll]

################# Calculate timeseries #################
if source=='ocean':
    print('Data is low frequency: multiplying by footprint to estimate mol fraction')
    # reindex the flux to match the footprints
    # multiply by the footprint to give a map of the contirbution to mf a a site
    mf_map = {spec: em_spec.reindex_like(fp_ds.fp, method='ffill') * fp_ds.fp
              for spec, em_spec in emissions.items()}
    # sum the map to get the mf at a site
    mf_ts = {spec: mf.flux.sum(dim=['lat', 'lon']) / concentration('ppm') for spec, mf in mf_map.items()}
    for spec, mf in mf_map.items():
        if 'uncertainty' in mf.data_vars:
            mf_ts[f'{spec}_uncertainty'] = (mf.uncertainty**2).sum(dim=['lat', 'lon'])**0.5

    ts_tuple = {spec: (['time'], mf.compute().values,
                       {'units': '1e-6',
                        'flux file': f'{emissions_names[spec]}_EUROPE_{year}.nc'})
                for spec, mf in mf_ts.items()}
    mf_ts = xr.Dataset(data_vars = ts_tuple,
                       coords = {'time': fp_ds.time.values})

# calculate the mf contribution map using the HiTRes information
else:
    print('Data is high frequency: using timeseries_HiTRes to estimate mol fraction')
    emissions = {spec: {'high_freq': fl_spec['high_freq'].reindex_like(fp_ds.fp, method='ffill')}
                for spec, fl_spec in emissions.items()}
    mf_ts = timeseries_HiTRes(flux_dict = emissions,
                              fp_HiTRes_ds = fp_ds, 
                              output_TS = True,
                              output_fpXflux = False)
    for dv in mf_ts.data_vars:
        mf_ts[dv] = mf_ts[dv] / concentration('ppm')
        mf_ts[dv] = mf_ts[dv].assign_attrs({'units': '1e-6',
                                            'flux file': f'{emissions_names[spec]}_EUROPE_{year}.nc'})
    if 'total' in mf_ts.data_vars:
        mf_ts = mf_ts.rename({'total': list(emissions.keys())[0]})

################# Save results #################
ff_str = f'_{ff_model}' if ff_model!='edgar-ukghg' else ''
ox_ratio = ff_str if oxidativeratio is None else f'{ff_str}-{oxidativeratio}' if ff_str!='' else f'_{oxidativeratio}'
model_str = f'_{bio_model}' if source=='bio' else f'{ox_ratio}' if source=='ff' else ''
date_str = f'climatology_{year}' if climatology and month is None else \
           f'climatology_{year}{str(month).zfill(2)}' if climatology else \
           f'{year}' if month is None else \
           f'{year}{str(month).zfill(2)}'
filename = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2',
                        f'{site}_{source}{model_str}_timeseries_{date_str}.nc')

# assign attributes to the dataset
attrs = {'description': f'{source} mol fraction timeseries for {site}'}
if source=='ff':
    attrs['inventory'] = ff_model
    if oxidativeratio is not None:
        attrs['oxidative_ratio'] = oxidativeratio
    else:
        attrs['o2 description'] = 'calculated by applying oxidative ratios to each sector'
elif source=='bio':
    attrs['biospheric model'] = bio_model
else:
    attrs['ECCO year'] = 'climatology' if year>2018 or climatology else year
    attrs['NEMO year'] = 'climatology' if year>2015 or climatology else year
    attrs['Jena CarboScope year'] = 'climatology' if year>2015 or climatology else year
    print('')
if fp_species is not None:
    attrs['fp_species'] = fp_species
mf_ts = mf_ts.assign_attrs(attrs)

# save to netcdf
print(f'\n\n\nSaving to: {filename}')
mf_ts.to_netcdf(filename)
print('\nDone')