import os
import sys
import xarray as xr

from acrg.name import name

species = ['apo', 'co2', 'o2', 'delta_o2_n2']

# import variables from bash script
year = int(sys.argv[1])
month = None if str(sys.argv[2])=='none' else int(sys.argv[2])
site = str(sys.argv[3])
height = '185magl' if site=='TAC' else None

units = {'apo': 'per meg', 'co2': 'ppm', 'o2': 'ppm', 'delta_o2_n2': 'per meg'}


# define the start and end dates
start_month = 1 if month is None else month
end_month = 1 if month is None or month==12 else month+1
end_year = year+1 if month is None or month==12 else year

start = f'{year}-{str(start_month).zfill(2)}-01'
end = f'{end_year}-{str(end_month).zfill(2)}-01'

baseline = {}
for spec in species:
    print(f'\n{spec}')
    ########################### GET FOOTPRINTS ###########################
    print('Getting footprint')
    footprint = name.footprints(
                            site, met_model='UKV',
                            start = start,
                            end = end,
                            height = height,
                            domain = 'EUROPE',
                            species = 'co2',
                            chunks = {'time': 50})
    fpdm = {site: footprint}

    ################## Calculate boundary conditions ###################
    print('Loading boundary conditions')
    species_bc = 'co2-jenacarboscope' if spec.lower()=='co2' else spec
    bc = name.boundary_conditions(
                                domain='EUROPE',
                                species = species_bc,
                                start = start,
                                end = end,
                                chunks = {'time': 50})

    print('Calculating boundary conditions contribution')
    fpdm = name.add_bc(
                    fp_and_data = fpdm,
                    load_bc = True,
                    species = spec,
                    bc = bc)
    baseline_spec = fpdm[site].bc
        
    if spec=='apo':
        print('Converting to per meg')
        baseline_spec = baseline_spec * (1e6 / 209460)

        if all(abs(baseline_spec)>1e7):
            baseline_spec = baseline_spec * 1e-6
        
    elif all(abs(baseline_spec<1e-3)):
        print('Converting to ppm')
        baseline_spec = baseline_spec * 1e6

    baseline[spec] = baseline_spec

# calculate the delta(O2/N2) and oxygen timeseries
if all([spec in species for spec in ['co2', 'apo']]):
    baseline['delta_o2_n2'] = baseline['apo'] - 1.1 * (1e6 / 209460) * (baseline['co2'] - 350)
    baseline['o2'] = (baseline['delta_o2_n2'] * (209460 / 790190) * 1e-6 + (209460/790190)) * 790190

# convert to tuples for creating dataset
for spec, ts in baseline.items():
    baseline[spec] = (['time'], ts.values,
                      {'description': f'modelled boundary condition contribution to the mole fraction of {spec}',
                       'units': units[spec]})

# convert to dataset
fm_bc = xr.Dataset(data_vars = baseline,
                   coords = {'time': fpdm[site].time},
                   attrs = {'description': 'Baseline calculated from boundary conditions using the Jena Carboscope global simulation'})

# create filenames
fname = f'{site}_bc_timeseries_{year}.nc' if month is None else \
        f'{site}_bc_timeseries_{year}{str(month).zfill(2)}.nc'
filename = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2', fname)

# save to netcdf
print(f'\n\n\nSaving to: {filename}')
fm_bc.to_netcdf(filename)
print('\nDone')

