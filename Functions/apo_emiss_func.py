import os
import re
import glob
import datetime
import regionmask
import numpy as np
import pandas as pd
import xarray as xr
from calendar import monthrange

from acrg.name import name
from acrg.name.flux import write

import apo_funcs

def mask_land_ocean(lat, lon, mask_region):
    '''
    Create a mask to mask out either the land or the ocean

    Args
        lat, lon (xarray.DataArray)
            latitud and longitude grid to apply the mask to
        mask_region (str)
            region to mask out, either 'land' or 'ocean'
    
    Returns
        xarray.DataArray
            array with 0s for the masked regions and 1s for the unmasked regions
    '''
    masked = regionmask.defined_regions.natural_earth_v5_0_0.land_110 \
             if mask_region=='ocean' else \
             regionmask.defined_regions.natural_earth_v5_0_0.ocean_basins_50 \
             if mask_region=='land' else None
    if masked is None:
        print('Mask must either be land or ocean')
        return None
    
    mask = masked.mask(lon, lat)
    mask = mask.where(mask.isnull(), other=1)

    return mask

def climatology(data, start_year=None, end_year=None, res=None, flux_var=None, uncertainty=False, save=True, write_kwargs=None):
    '''
    Calculate a climatology from the data

    Args
        data (xarray.Dataset)
            the data from which to create a climatology
        start_year, end_year (float or int or str, optional)
            the first and last year to include data from for the climatology
        res (str, optional)
            the resolution of the data
        flux_var (str, optional)
            the name of the variable for which to create a climatology
        uncertainty (bool, optional)
            whether to calculate the variability across the years
            if True, calculates the standard deviation for each time step
        save (bool, optional)
            whether to save the climatology to netcdf
            if False the xarray Dataset is returned
        write_kwargs (dict, optional)
            keyword arguments for the function name.flux.write
    '''
    # get the first and last year to include data from in the climatology calculation
    start_year = pd.to_datetime(data.time[0].values).year if start_year is None else int(start_year)
    end_year = pd.to_datetime(data.time[-1].values).year if end_year is None else int(end_year)
    print(f'Calculating climatology using data between {start_year} and {end_year}')

    # calculate the time resolution
    if res is None:
        # get the data time step in hours
        time_step = int(data.time.diff(dim="time").values.mean()) / (1e9*60*60)
        print(time_step)
        # return a string of what the time step
        res = 'M' if round(time_step)==729 else 'D' if time_step==24 else 'h' if time_step==1 else None
        if res is None:
            print('Unknown time step')
            return None

    flux_var = 'flux' if flux_var is None else flux_var
    # remove leap days
    data = data.sel(time=~((data.time.dt.month == 2) & (data.time.dt.day == 29)))

    # slice the data into years
    data_years = {year: data[flux_var].sel(time=slice(f'{year}-01-01', f'{year}-12-31'))
                  for year in range(start_year, end_year+1)}

    # replace the datetime information on each year to be consistent
    time_climatology = np.arange('1900-01-01', '1901-01-01', dtype=f'datetime64[{res}]')
    for  data_year in data_years.values():
        data_year['time'] = time_climatology
    
    # add an extra dimension for the year and combine all arrays into one along this dimension
    data_years_da = {yr: n_yr.expand_dims({'year':[yr]}) for yr, n_yr in data_years.items()}
    data_years_da = xr.concat([n_yr for n_yr in data_years.values()], dim='year')

    # average the arrays across all years
    data_ave = data_years_da.mean(dim='year').to_dataset(name='flux')
    if uncertainty:
        data_ave['standard_deviation'] = data_years_da.std(dim='year')
    
    if save:
        # set some default inputs for write_kwargs if not given
        write_kwargs = dict() if write_kwargs is None else write_kwargs
        if 'flux_comments' not in write_kwargs.keys():
            write_kwargs['flux_comments'] = f'Climatology from the years {start_year} - {end_year}'
        if 'uncertainty_desc' not in write_kwargs.keys() and uncertainty:
            write_kwargs['uncertainty_desc'] = f'Climatology standard deviation from the years {start_year} - {end_year}'
        if 'regridder_used' not in write_kwargs.keys() and 'regridder_used' in data.attrs.keys():
            write_kwargs['regridder_used'] = data.attrs['regridder_used']
        if 'source' in write_kwargs and res!='M':
            write_kwargs['source'] = f"climatology-{write_kwargs['source']}" \
                                     if 'climatology' not in write_kwargs['source'] else write_kwargs['source']
        if 'prior_info_dict' not in write_kwargs.keys() and 'prior_file_1' in data.attrs.keys():
            write_kwargs['prior_info_dict'] = apo_funcs.get_prior_info(attributes = data.attrs)
            # write_kwargs['prior_info_dict'] = {data.attrs['prior_file_1']: [data.attrs['prior_file_1_version'],
            #                                                                 data.attrs['prior_file_1_raw_resolution'],
            #                                                                 data.attrs['prior_file_1_reference']]}
        if 'domain' not in write_kwargs.keys():
            write_kwargs['domain'] = 'EUROPE'
        
        data_std = data_ave['standard_deviation'] if uncertainty else None
        climatology = True if res=='M' else False
        try:
            write(lat = data_ave.lat.values,
                  lon = data_ave.lon.values,
                  time = data_ave.time.values,
                  flux = data_ave.flux.values,
                  uncertainty = data_std,
                  climatology = climatology,
                  **write_kwargs)
        except TypeError:
            print('Arguments missing from write_kwargs')
            return write_kwargs, data_ave
    else:
        return data_ave

def del_o2_n2(o2, n2, o2_base=None, o2_n2_base=None):
    '''
    Combine O2 and N2 fluxes into the delta(O2/N2) term

    Args
        o2 (numpy ndarray or xarray DataArray)
            oxygen flux, change in O2
        n2 (numpy ndarray or xarray DataArray)
            nitrogen flux, change in N2
            should have the same length as o2
        o2_base (float, int, numpy ndarray, or xarray DataArray)
            baseline oxygen value
            if an array, it should have the same length as o2
        o2_n2_base (float, int, numpy ndarray, or xarray DataArray)
            baseline delta(o2/n2) value
            if an array, it should have the same length as o2
    
    Returns
        numpy ndarray or xarray DataArray
            array of equal length and type as o2 and n2
    '''
    o2_n2_ref = 209460 / 790190
    o2_base = 0 if all([base is None for base in [o2_base, o2_n2_base]]) else \
              o2_base if o2_base is not None else \
              (o2_n2_base * o2_n2_ref * 1e-6 + o2_n2_ref) * 790190
    o2_n2_sample = (o2 + o2_base) / (n2 + 790190)

    delta_sample = (o2_n2_sample - o2_n2_ref) / o2_n2_ref * 1e6

    return delta_sample

## https://web.uvic.ca/~rhamme/Hamme_gas_solubility.pdf
solubility_coeffs = {'A_coeffs' : {'nitrogen': [6.42931, 2.92704, 4.32531, 4.69149]},
                     'B_coeffs' : {'nitrogen': [-7.44129e-3, -8.02566e-3, -1.46775e-2]}
                    }

# solubility_coeffs = {'A_coeffs' : {'nitrogen': [5.81, 3.2, 4.17, 5.1]},
#                      'B_coeffs' : {'nitrogen': [-7e-3, -7.7e-3, -1.14e-2]}
#                     }

def heat_capacity(specific_heat_capacity):
    '''
    Calculate the heat capacity for a substance given its specific heat
    
    Args
        specific_heat_capacity (int or float)
            specific heat of substance
    
    Returns
        int or float
            heat capacity, J kg-1 K-1
    
    '''
    heat_capacity_water = 4184
    return heat_capacity_water * specific_heat_capacity

def sea_heat_flux(heat_flux, solubility_deriv=None,
                  sea_heat_capacity=None, specific_heat_capacity=None, **kwargs):
    '''
    Estimate the sea-air flux of a gas from the sea-air heat flux
    Keeling 1993, equation 19
    
    Args
        heat_flux (numpy ndarray)
            net air-to-sea heat flux, J m-2 s-1
        solubility_deriv (float, int)
            temperature derivatie of solubility, mol m-3 K-1
            dC_eq/dT
            C_eq : mixed layer concentration in equilibrium with the atmosphere
            if not given then args for solubility_derivative() must be given
            e.g. sea_heat_flux(heat_flux, specific_heat_capacity=0.9,
                               temperature=10, salinity=35, species='nitrogen')
        heat_capacity
            heat capacity of seawater, J m-3 K-1
        
        denominator of solubility_derivative and heat_capacity can be different but should match
        
    '''
    if solubility_deriv is None:
        # if a value for the solubilty is not given, it can be calculated from the keyword args
        # these must be given as kwargs
        args = {arg : None if arg not in kwargs else kwargs[arg]
                for arg in ['A_coeffs', 'B_coeffs', 'species', 'out_units']}
        solubility_deriv = solubility_derivative(temperature = kwargs['temperature'],
                                                 salinity    = kwargs['salinity'],
                                                 species     = args['species'],
                                                 A_coeffs    = args['A_coeffs'],
                                                 B_coeffs    = args['B_coeffs'],
                                                 out_units   = args['out_units'])
        
    sea_heat_capacity = heat_capacity(specific_heat_capacity) if sea_heat_capacity is None else sea_heat_capacity
    
    flux = solubility_deriv * heat_flux / sea_heat_capacity
    
    return flux

def solubility_oxygen(temperature, salinity):
    a = 14.161 - 0.3943 * temperature + 7.714e-3 * temperature**2 - 6.46e-5 * temperature**3
    b = 0.0841 - 2.56e-3 * temperature + 3.74e-5 * temperature**3
    
    return a - salinity * b

def solubility(temperature, salinity, A_coeffs=None, B_coeffs=None, species=None,
               out_units='mol_L'):
    '''
    Estimate gas solubility in seawater, dependent on temperature
    C_g is the gas concentration in micro mol/L
    https://web.uvic.ca/~rhamme/Hamme_gas_solubility.pdf
    
    Lower temperature --> higher solubility
    
    Args
        temperature (int or float or numpy.ndarray)
            water temperature, degrees celsius
        salinity (int or float)
            salinity of water, PSS
            0 for pure water
        A_coeffs, B_coeffs (list or numpy.ndarray)
            Coefficients for solubility equation
            e.g. (Hamme & Emerson, 2004)
            N2: A_coeffs = [6.42931, 2.92704, 4.32531, 4.69149]
                B_coeffs = [-7.44129e-3, -8.02566e-3, -1.46775e-2]
        species (str)
            gas species
            Used to get the A_coeffs and B_coeffs, if they are in the dict above
        out_units (str)
            units for C_g, dafults to mol/L
    
    Returns
        float or numpy.ndarray
            Gas solubility for given temperature and salinity
            If temperature is an array, an array will returned
    '''
    out_units = 'mol_L' if out_units is None else out_units
    
    T_s     = np.log((298.15 - temperature) / (273.15 + temperature))
    
    A_coeffs = A_coeffs if A_coeffs is not None else solubility_coeffs['A_coeffs'][species]
    B_coeffs = B_coeffs if B_coeffs is not None else solubility_coeffs['B_coeffs'][species]
    
    log_C_1 = A_coeffs[0] + A_coeffs[1]*T_s + A_coeffs[2]*T_s**2 + A_coeffs[3]*T_s**3
    log_C_2 = B_coeffs[0] + B_coeffs[1]*T_s + B_coeffs[2]*T_s**2
    
    log_C   = log_C_1 + salinity * log_C_2
    
    C = np.exp(log_C) if out_units=='micro_mol_L' else \
        np.exp(log_C) * 1e-6
    
    return C

def solubility_derivative(temperature, salinity, A_coeffs=None, B_coeffs=None, species=None,
                          out_units='mol_L'):
    '''
    Derivative of the equation for gas solubility in seawater, above
    
    Args
        temperature (int or float or numpy.ndarray)
            water temperature, Kelvin
        salinity (int or float)
            salinity of water, PSS
            0 for pure water
        A_coeffs, B_coeffs (list or numpy.ndarray)
            Coefficients for solubility equation
            e.g. (Hamme & Emerson, 2004)
            N2: A_coeffs = [6.42931, 2.92704, 4.32531, 4.69149]
                B_coeffs = [-7.44129e-3, -8.02566e-3, -1.46775e-2]
        species (str)
            gas species
            Used to get the A_coeffs and B_coeffs, if they are in the dict above
        out_units (str)
            units for C_g, dafults to mol/L
            note : L = kg
            
    Returns
        float or numpy.ndarray
            Derivative of the gas solubility for given temperature and salinity
            If temperature is an array, an array will returned
    '''
    out_units = 'mol_L' if out_units is None else out_units

    A_coeffs = A_coeffs if A_coeffs is not None else solubility_coeffs['A_coeffs'][species]
    B_coeffs = B_coeffs if B_coeffs is not None else solubility_coeffs['B_coeffs'][species]
    
    d_T_s_d_T = (571.3) / ((temperature-571.3) * temperature)
    T_s       = np.log((571.3 - temperature) / temperature)
    
    # by the chain rule: dy/dx = dy/du * du/dx
    d_log_C_1_d_t = (A_coeffs[1] + 2*A_coeffs[2]*T_s + 3*A_coeffs[3]*T_s**2) * d_T_s_d_T
    d_log_C_2_d_t = (B_coeffs[1] + 2*B_coeffs[2]*T_s) * d_T_s_d_T * salinity
    
    # differentiating exponents: y = e^x, dy/dx = e^x
    # chain rule: dy/dt = dy/dx * dx/dt = e^x * dx/dt
    C_g = solubility(temperature, salinity, A_coeffs, B_coeffs, species, out_units='micro_mol_L')
    d_C_d_t = (d_log_C_1_d_t + d_log_C_2_d_t) * C_g
    
    d_C_d_t = d_C_d_t if out_units=='micro_mol_L' else \
              d_C_d_t * 1e-6
    
    return d_C_d_t

def apo_base(oxygen_sample, nitrogen_sample, oxygen_ref, nitrogen_ref,
             mf_oxygen, mf_co2, oxidative_ratio_bio):
    rat_o2_n2 = {'sample' : oxygen_sample / nitrogen_sample,
                 'ref'    : oxygen_ref / nitrogen_ref}
    
    del_o2_n2 = ((rat_o2_n2['sample'] / rat_o2_n2['ref']) - 1) * 1e6
    
    del_apo = del_o2_n2 + (oxidative_ratio_bio / mf_oxygen) * (mf_co2 - 363.29)
    
    return del_apo

def apo_species_split(ocean_o2=None, ff_co2=None, ff_o2=None,
                      ocean_co2=None, ocean_n2=None,
                      oxidative_ratio_ff=None, oxidative_ratio_bio=1.1,
                      mf_o2=0.2094, mf_n2=0.78084, convert=False):
    '''
    Args:
        ocean_o2 (numpy.ndarray or xarray.DataArray)
            field of the ocean oxygen flux
        ff_co2 (numpy.ndarray or xarray.DataArray)
            field of the fossil fuel co2 flux
        ocean_co2 (numpy.ndarray or xarray.DataArray)
            field of the ocean co2n flux
        ocean_n2 (numpy.ndarray or xarray.DataArray)
            field of the ocean nitrogen flux
        mf_oxygen (float or int)
            mol fraction of oxygen in dry air
            default = 0.2904
        mf_n2 (float or int)
            mol fraction of nitrogen in dry air
            default = 0.78084
    '''    
    apo_split = {}
    oxidative_ratio_bio = -oxidative_ratio_bio if oxidative_ratio_bio>0 else oxidative_ratio_bio
    
    if ff_co2 is not None and oxidative_ratio_ff is not None:
        oxidative_ratio_ff = -oxidative_ratio_ff if (oxidative_ratio_ff.fillna(1)>0).all() else \
                             oxidative_ratio_ff
        apo_split['co2_ff'] = (oxidative_ratio_ff - oxidative_ratio_bio) * ff_co2 / mf_o2
    elif ff_co2 is not None and ff_o2 is not None:
        ff_o2 = -ff_o2 if (ff_o2.fillna(1)>0).all() else ff_o2
        apo_split['co2_ff'] = (ff_o2 - oxidative_ratio_bio * ff_co2) / mf_o2

    if ocean_co2 is not None:
        apo_split['co2_ocean'] = oxidative_ratio_bio * ocean_co2 / mf_o2

    if ocean_o2 is not None:
        apo_split['o2_ocean'] = ocean_o2 / mf_o2

    if ocean_n2 is not None:
        apo_split['n2_ocean'] = -ocean_n2 / mf_n2
    
    if convert:
        apo_split = {source: apo_source * 1e6 for source, apo_source in apo_split.items()}

    return apo_split

def apo(ocean_o2, ff_co2, ocean_co2, ocean_n2, ff_o2=None,
        oxidative_ratio_ff=None, oxidative_ratio_bio=1.1,
        mf_o2=0.2094, mf_n2=0.78084, land_ocean_split=False,
        convert=False):
    '''
    Args:
        ocean_o2 (numpy.ndarray)
            field of the ocean oxygen flux
        ff_co2 (numpy.ndarray)
            field of the fossil fuel co2 flux
        ocean_co2 (numpy.ndarray)
            field of the ocean co2n flux
        ocean_n2 (numpy.ndarray)
            field of the ocean nitrogen flux
        mf_oxygen (float or int)
            mol fraction of oxygen in dry air
            default = 0.2904
        mf_n2 (float or int)
            mol fraction of nitrogen in dry air
            default = 0.78084
    '''
    apo_split = apo_species_split(ocean_o2 = ocean_o2,
                                  ff_co2 = ff_co2,
                                  ocean_co2 = ocean_co2,
                                  ocean_n2 = ocean_n2,
                                  ff_o2 = ff_o2,
                                  oxidative_ratio_ff = oxidative_ratio_ff,
                                  oxidative_ratio_bio = oxidative_ratio_bio,
                                  mf_o2 = mf_o2,
                                  mf_n2 = mf_n2,
                                  convert = convert)

    if apo_split['o2_ocean'].time.diff(dim="time").values.mean() != apo_split['co2_ocean'].time.diff(dim="time").values.mean():
        apo_ocean_o2 = apo_funcs.combine_diff_resolution(apo_split['o2_ocean'], apo_split['co2_ocean'],
                                                        method = 'add', verbose = 'False')
    else:
        apo_ocean_o2 = apo_split['o2_ocean'] + apo_split['co2_ocean']

    if apo_ocean_o2.time.diff(dim="time").values.mean() != apo_split['n2_ocean'].time.diff(dim="time").values.mean():
        apo_ocean = apo_funcs.combine_diff_resolution(apo_ocean_o2, apo_split['n2_ocean'],
                                                     method = 'add', verbose = 'False')
    else:
        apo_ocean = apo_ocean_o2 + apo_split['n2_ocean']
    
    if land_ocean_split:
        return {'ocean': apo_ocean, 'land': apo_split['co2_ff']}
    
    else:
        if apo_ocean.time.diff(dim="time").values.mean() != apo_split['co2_ff'].time.diff(dim="time").values.mean():
            apo = apo_funcs.combine_diff_resolution(apo_ocean, apo_split['co2_ff'],
                                                   method = 'add', verbose = 'False')
        else:
            apo = apo_ocean + apo_split['co2_ff']
        return apo

def daily_to_hourly(data, timescale=1, year=None, lat=None, lon=None, time=None):
    '''
    Split monthly data into smaller chunks
    
    Args:
        data (numpy.ndarray)
            monthly data to be split
        timescale (int or float)
            timescale to be used for split data in hours, defaults to 1 hour
        year (int or float)
            year of data
        lat, lon, time (numpy.ndarray)
            coordinates for output DataArray
            if all of these are given the split data will be returned as an xarray.DataArray
            if any of these are not given the split data will be returned as an numpy.ndarray
    
    Returns:
        numpy.ndarray or xarray.DataArray
            array of data split into new timescales
            if a DataArray is returned then it will have the lat, lon, and time coords given
            data shape:
            xarray.DataArray : lat, lon, time
            nump.ndarray : time, lat, lon
    '''
    if 'xr' not in dir(): import xarray as xr
    if 'np' not in dir(): import numpy as np
    # define the month lengths
    days = 365 if year is None or year%4!=0 else 366
    
    # calculate the number of splits required
    num_split  = int(24 / timescale)
    
    flux_split = [data[:,:,dd].tolist() for dd in range(days)]
    flux_split = ([np.array([flux_split[dd]] * num_split) for dd in range(days)])
    flux_split = np.concatenate(flux_split, axis=0)
    
    if all([coord is not None for coord in [lat, lon]]):
        if time is None and year is not None:
            time = np.arange(np.datetime64(f'{year}-01-01T00:00'),
                             np.datetime64(f'{year+1}-01-01T00:00'),
                             timescale, dtype='datetime64[h]')
            
        elif time is None and year is None:
            print('A time array or year must be given to create an xarray.DataArray')
            print('Returning split data as a numpy.ndarray')
            return flux_split
        
        flux_da  = xr.DataArray(data   = flux_split,
                                coords = dict(time = time, lat = lat, lon = lon),
                                dims   = ['time', 'lat', 'lon'])
        flux_da  = flux_da.transpose('lat','lon','time')
        
        return flux_da
    
    else:
        return flux_split
    
def monthly_to_hourly(data, timescale=1, year=None, lat=None, lon=None, time=None):
    '''
    Split monthly data into smaller chunks
    
    Args:
        data (numpy.ndarray)
            monthly data to be split
        timescale (int or float)
            timescale to be used for split data in hours, defaults to 1 hour
        year (int or float)
            year of data
        lat, lon, time (numpy.ndarray)
            coordinates for output DataArray
            if all of these are given the split data will be returned as an xarray.DataArray
            if any of these are not given the split data will be returned as an numpy.ndarray
    
    Returns:
        numpy.ndarray or xarray.DataArray
            array of data split into new timescales
            if a DataArray is returned then it will have the lat, lon, and time coords given
            data shape:
            xarray.DataArray : lat, lon, time
            nump.ndarray : time, lat, lon
    '''
    if 'xr' not in dir(): import xarray as xr
    if 'np' not in dir(): import numpy as np
    # define the month lengths
    feb    = 28 if year is None or year%4!=0 else 29
    months = [31, feb, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    
    # calculate the number of splits required
    num_split  = int(24 / timescale)
    
    flux_split = [data[:,:,mm].tolist() for mm, month in enumerate(months)]
    flux_split = ([np.array([flux_split[mm]]*(month*num_split)) for mm, month in enumerate(months)])
    flux_split = np.concatenate(flux_split, axis=0)
    
    if all([coord is not None for coord in [lat, lon]]):
        if time is None and year is not None:
            time = np.arange(np.datetime64(f'{year}-01-01T00:00'),
                             np.datetime64(f'{year+1}-01-01T00:00'),
                             timescale, dtype='datetime64[h]')
            
        elif time is None and year is None:
            print('A time array or year must be given to create an xarray.DataArray')
            print('Returning split data as a numpy.ndarray')
            return flux_split
        
        flux_da  = xr.DataArray(data   = flux_split,
                                coords = dict(time = time, lat = lat, lon = lon),
                                dims   = ['time', 'lat', 'lon'])
        flux_da  = flux_da.transpose('lat','lon','time')
        
        return flux_da
    
    else:
        return flux_split
    
# @jit(nopython=True)
def TNO_emissions_map(emissions_source, map_shape, longitude_index, latitude_index,
                      adjust=0, check_fails=False, verbose=True, updates_perc=1,
                      cat_index=None, split=False, cat_index_req=None):
    '''
    Convert a 1d array of emissions with a source dimension to a 2d array with 
    latitude and longitude dimensions
    
    Args
        emissions_source (numpy.ndarray or list)
            emissions for each source
        map_shape (tuple)
            shape of output map, should have latitude and longitude axes
            should be (lat_len, lon_len)
        longitude_index, latitude_index (numpy.ndarray or list)
            index of the longitude and latitude of each source
            should have length equal to emissions_source
        adjust (int)
            adjustment to account for 0 or 1 indexing of longitude_index and latitude_index
            i.e. if longitude_index and latitude_index are 1 indexed adjust = -1
        verbose (bool)
            whether to print updates on progress
        updates_perc (int)
            if verbose==True, this is the percentage completion at which updates will be printed
        cat_index (np.ndarray or list)
            indices of emissions categories
        cat_index_req (int)
            index of emissions category required - used to create source specific maps,
            if None then a total emissions map is produced
    
    Returns
        dict
            keys will be the source indices if split==True or 'total' if not
            values will be emission maps
    '''
    cats          = list(set(cat_index)) if cat_index is not None and cat_index_req is None else \
                    [cat_index_req] if cat_index is not None else None
    
    # create a 2d map with lat & lon coords to fill with emissions values
    emissions_map = {'tot': np.zeros(map_shape)} if not split or cat_index_req is None else \
                    {cat_index_req: np.zeros(map_shape)} if cat_index_req is not None else \
                    {cat : np.zeros(map_shape) for cat in cats}
    perc          = int(updates_perc * 0.01 * len(emissions_source))
    
    # loop through emissions sources to add to map
    for ee, emission in enumerate(emissions_source):
        # get the index of the lat & lon for the source
        ind_lat = latitude_index  + adjust
        ind_lon = longitude_index + adjust
        ind_cat = cat_index[ee] if cat_index is not None else 'tot'
        
        fails = []
        try:
            # add emissions to map
            if cat_index_req is None or ind_cat==cat_index_req:
                emissions_map[ind_cat][ind_lat][ind_lon] = emissions_map[ind_cat][ind_lat][ind_lon] + float(emission)
                
        except:
            fails.append(ee)
            
        if verbose and ee%perc==0: print(updates_perc, 'per cent completed')
        else: pass
            
    if verbose: print(len(fails), 'fails')
    else: pass
    
    return emissions_map


class TNO_emissions:
    def __init__(self, path, year, check_files=False):
        '''
        Args
            path (str)
                path to directory containing inventory data files
            year (str or int)
                year of data required
            check_files (bool, optional)
                if True, the filenames are added to the object
        '''
        self.year = int(year)
        self.path = path
        
        filenames = glob.glob(os.path.join(path, f'*{year}*.nc'))
        if check_files: self.filenames = filenames
        
        if len(filenames)>0:
            self.emissions_ds = name.open_ds(filenames[0])
        else:
            print(f'Cannot find TNO file for {year} in directory : {path}')
        
    
    def get_time_profiles(self, path_timeprofiles=None, timeprofile_names=None):
        '''
        Args
            path_timeprofiles (str, optional)
                path to directory containing time profile files
                only needed if different to self.path
            timeprofile_names (dict, optional)
                timeprofile_dict keys which describe the scaling for each dataframe
                keys should be 'month', 'day', 'hour'
                defaults to {'month' : 'months-in-year', 'day' : 'days-in-week', 'hour' : 'hours-in-day'}
        Returns
            pandas.Dataframe
                timeprofile dataframes for disaggregating into monthly, daily,
                and hourly timescales with keys matching time_names values
        '''
        path_timeprofiles = path_timeprofiles if path_timeprofiles is not None else self.path
        self.timeprofile_names = {'month' : 'months-in-year', 'day' : 'days-in-week', 'hour' : 'hours-in-day'} if \
                                  timeprofile_names is None else timeprofile_names
        
        time_profile_files = {name : glob.glob(os.path.join(path_timeprofiles, 'timeprofiles', f'*{name}*.csv'))[0]
                              for name in self.timeprofile_names.keys()}
        df_timeprofile     = {name : pd.read_csv(time_profile_files[name], sep=';', header=6, index_col=2, engine='python')
                              for name in self.timeprofile_names.keys()}
        df_timeprofile     = {name : df_timeprofile[name].drop(labels=df_timeprofile[name].columns.values[:2], axis=1)
                              for name in self.timeprofile_names.keys()}

        # reformat time period names to take spaces out
        time_periods       = {name : df_timeprofile[name].columns.values for name in self.timeprofile_names.keys()}
        new_headers        = {name : {time : re.sub(' ', '', time) for time in time_periods[name]}
                              for name in self.timeprofile_names.keys()}
        [df_timeprofile[name].rename(columns=new_headers[name], inplace=True) for name in self.timeprofile_names.keys()];
        
        self.timeprofiles = df_timeprofile
        self.sectors      = {str(cat): ind+1 for ind, cat in
                             enumerate(self.timeprofiles['month'].index)}
    
    
    def disaggregate_TNO(self, sector, species='co2_ff', resolution='hourly', data_type='map', flux_var=None,
                         save_lower_res=False, path_timeprofiles=None, timeprofile_names=None):
        '''
        Disaggregate annual TNO emissions estimates into higher time resolution

        Args
            sector (str)
                sector name e.g. 'A_PublicPower'
                should correspond to row names in self.timeprofiles dataframes
            species (str, optional)
                species to disaggregate, defaults to 'co2_ff'
            resolution (str, optional)
                time resolution to disggregate into
                one of monthly, daily, or hourly (case insensitive)
                defaults to hourly
            data_type (str, optional)
                form of data to be disaggregated, either 'map' or 'array'
                'map'   : emissions are in a map of latitude & longitude
                'array' : emissions are in a 1D array for each source
                defaults to 'map'
            flux_var (str, optional)
                name of the flux variable in the dataset
                if None, will use 'flux'
            save_lower_res (bool, optional)
                if True then lower time resolution disaggregations will be saved
                as well as the required resolution
                defaults to False
            path_timeprofiles (str, optional)
                path to directory containing time profile files
                only needed if different to self.path
            timeprofile_names (dict, optional)
                timeprofile_dict keys which describe the scaling for each dataframe
                keys should be 'month', 'day', 'hour'
                defaults to {'month' : 'months-in-year', 'day' : 'days-in-week', 'hour' : 'hours-in-day'}
        Returns
            dict
                disaggregated emissions into month names, day names, and hour numbers
                (depending on required resolution)
        '''
        
        self.get_time_profiles(path_timeprofiles=path_timeprofiles, timeprofile_names=timeprofile_names)
        
        self.resolution  = resolution
        
        ### get the emissions for the correct sector and species
        if data_type.lower()=='array':
            sector_index     = self.sectors[sector]

            emissions_annual = getattr(self.emissions_ds, species).values
            emissions_annual = emissions_annual[self.emissions_ds.emission_category_index.values==sector_index]
            self.source_number = self.emissions_ds.source[self.emissions_ds.emission_category_index.values==sector_index]
        elif data_type.lower()=='map':
            flux_var         = 'flux' if flux_var is None else flux_var
            emissions_annual = getattr(self.emissions_ds, flux_var)
        else:
            print('Unknown data type, shuld be either map or array')
        
        sector_split = sector.split(' ') if ' ' in sector else sector.split('_') if '_' in sector else [sector]
        sector       = [sec for sec in self.sectors.keys() if all([string in sec for string in sector_split])][0]
        ### get the timeprofiles for the correct sector
        timeprofile      = {time : self.timeprofiles[time].loc[sector] for time in self.timeprofiles.keys()}
        
        ### disaggregate the annual emissions into monthly
        emissions_month  = {month : emissions_annual * 
                                    timeprofile['month'][timeprofile['month']].index[mm]/12.
                            for mm, month in enumerate(timeprofile['month'].index)}
        if resolution.lower()=='monthly' or save_lower_res:
            self.emissions_month = emissions_month
        
        ### disaggregate the monthly emissions into daily
        if resolution.lower() in ['hourly', 'daily']:
            emissions_day = {month :
                             {day: emissions_month[month] * 
                                   timeprofile['day'][timeprofile['day'].index[dd]]/7.
                              for dd, day in enumerate(timeprofile['day'].index)}
                             for month in timeprofile['month'].index}
            if resolution.lower()=='daily' or save_lower_res:
                self.emissions_day = emissions_day
        
        ### disaggregate the daily emissions into hourly
        if resolution.lower()=='hourly':
            emissions_hr  = {month :
                             {day: 
                              {hour: emissions_day[month][day] * 
                                     timeprofile['hour'][timeprofile['hour'].index[hh]]/24.
                               for hh, hour in enumerate(timeprofile['hour'].index)}
                              for day in timeprofile['day'].index}
                             for month in timeprofile['month'].index}
            
            self.emissions_hour = emissions_hr
    
    def time(self, resolution=None):
        resolution = self.resolution if resolution is None else resolution
        time_step  = 'M' if resolution.lower()=='monthly' else \
                     'D' if resolution.lower()=='daily' else 'h'
        self.time_array = np.arange(f'{self.year}-01-01', f'{self.year+1}-01-01', dtype=f'datetime64[{time_step}]')
    
    def disaggregate_TNO_time(self, sector, species='co2_ff', resolution='hourly', data_type='map',
                              flux_var=None, save_lower_res=False, path_timeprofiles=None,
                              timeprofile_names=None, save=False, filename=None):
        '''
        Disaggregate annual TNO emissions estimates into higher time resolution

        Args
            sector (str)
                sector name e.g. 'A_PublicPower'
                should correspond to row names in self.timeprofiles dataframes
            species (str, optional)
                species to disaggregate, defaults to 'co2_ff'
            resolution (str, optional)
                time resolution to disggregate into
                one of monthly, daily, or hourly (case insensitive)
                defaults to hourly
            data_type (str, optional)
                form of data to be disaggregated, either 'map' or 'array'
                'map'   : emissions are in a map of latitude & longitude
                'array' : emissions are in a 1D array for each source
                defaults to 'map'
            flux_var (str, optional)
                name of the flux variable in the dataset
                if None, will use 'flux'
            save_lower_res (bool, optional)
                if True then lower time resolution disaggregations will be saved
                as well as the required resolution
                defaults to False
            path_timeprofiles (str, optional)
                path to directory containing time profile files
                only needed if different to self.path
            timeprofile_names (dict, optional)
                timeprofile_dict keys which describe the scaling for each dataframe
                keys should be 'month', 'day', 'hour'
                defaults to {'month' : 'months-in-year', 'day' : 'days-in-week', 'hour' : 'hours-in-day'}
        Returns
            xarray.Dataset
                disaggregated emissions with an associated time coordinate
        '''
        if data_type.lower() not in ['map', 'array']:
            print('Data type unknown, should be either map or array')
            return None
        
        self.disaggregate_TNO(sector            = sector,
                              species           = species,
                              resolution        = resolution,
                              save_lower_res    = save_lower_res,
                              path_timeprofiles = path_timeprofiles,
                              timeprofile_names = timeprofile_names,
                              data_type         = data_type,
                              flux_var          = flux_var)
        self.time()
        
        if resolution.lower()=='monthly':
            # extract emissions from dict into a list
            emissions_array = self.emissions_month.values()
        
        else:
            ### the emissions dict is split into month, day of the week, and hour
            ### we need the day name associated with each date to determine the emissions throughout each month
            
            # list of day names and months
            days   = list(self.emissions_day['jan'].keys())*6 if resolution.lower()=='daily' else \
                     list(self.emissions_hour['jan'].keys())*6
            months = self.emissions_day.keys() if resolution.lower()=='daily' else self.emissions_hour.keys()
            
            # day of the week for which each month starts
            month_start_day = {month : datetime.datetime.strptime(f'01 {mm+1} {self.year}', '%d %m %Y').weekday()
                               for mm, month in enumerate(months)}
            # number of days in each month
            month_num_days  = {month : monthrange(self.year, mm+1)[1] for mm, month in enumerate(months)}
            # the names of the days of the week throughout each month
            # used to extract emissions from the dict
            month_days      = {month : days[month_start_day[month] : month_num_days[month] + month_start_day[month]]
                               for mm, month in enumerate(months)}
            
            if resolution.lower()=='daily':
                # convert emissions from a dict into a list and stack into a 2d array
                emissions_daily = [[emissions_days[day] for day in month_days[month]]
                                  for month, emissions_days in self.emissions_day.items()]
                emissions_array = np.row_stack(emissions) if data_type.lower()=='array' else \
                                  np.stack(emissions)
                
            else:
                # stack hourly emissions for each day into a 2d array
                emissions_hrly  = {month:
                                   {day: np.row_stack([emissions for emissions in emissions_day.values()])
                                         if data_type.lower()=='array' else
                                         np.stack([emissions for emissions in emissions_day.values()])
                                    for day, emissions_day in emissions_month.items()}
                                   for month, emissions_month in self.emissions_hour.items()}
                # convert emissions from a dict into a list and stack into a 2d array
                emissions_hrly  = [[emissions[day] for day in month_days[month]]
                                   for month, emissions in emissions_hrly.items()]
                emissions_array = np.row_stack([np.row_stack(emis) for emis in emissions_hrly]) \
                                  if data_type.lower()=='array' else \
                                  np.concatenate(np.concatenate(emissions_hrly))
        
        self.emissions_disag = emissions_array
        if data_type.lower()=='array':
            self.emissions_disag    = xr.Dataset({'flux': (['time', 'source'], emissions_array)},
                                         coords = {'time'   : self.time_array,
                                                   'source' : self.source_number},
                                                   attrs = self.emissions_ds.attrs)
        
        else:
            lat = self.emissions_ds.lat if 'lat' in self.emissions_ds else self.emissions_ds.latitude
            lon = self.emissions_ds.lon if 'lon' in self.emissions_ds else self.emissions_ds.longitude
            self.emissions_disag    = xr.Dataset({'flux': (['time','lat','lon'],
                                                        emissions_array)},
                                              coords = {'lat'  : lat,
                                                        'lon'  : lon,
                                                        'time' : self.time_array },
                                                        attrs = self.emissions_ds.attrs)
        
        if save:
            filename = os.path.join(outdir, f'{species}{sector}_EUROPE_{resolution}_{self.year}.nc')\
                       if filename is None else filename
        
            print(f'Saving map to {filename}')

            self.emissions_disag.to_netcdf(filename)
            
    
    def disaggregate_map_all_sectors(self, resolution='hourly', path_timeprofiles=None,
                                     timeprofile_names=None, filenames=None):
        '''
        Disaggregate TNO emissions maps for all sectors into a higher resolution and
        save to netcdf
        
        Args
            resolution (str, optional)
                time resolution to disggregate into
                one of monthly, daily, or hourly (case insensitive)
                defaults to hourly
            path_timeprofiles (str, optional)
                path to directory containing time profile files
                only needed if different to self.path
            timeprofile_names (dict, optional)
                timeprofile_dict keys which describe the scaling for each dataframe
                keys should be 'month', 'day', 'hour'
                defaults to {'month' : 'months-in-year', 'day' : 'days-in-week', 'hour' : 'hours-in-day'}
        '''
        sectors   = list(TNO_emiss_map.emissions_ds.data_vars)
        filenames = [None]*len(sectors) if filenames is None else filenames
        
        for ss, sector in enumerate(sectors):
            self.disaggregate_TNO_time(sector            = sector,
                                       resolution        = resolution,
                                       data_type         = 'map',
                                       save_lower_res    = False,
                                       path_timeprofiles = path_timeprofiles,
                                       timeprofile_names = timeprofile_names,
                                       save              = True,
                                       filename          = filenames[ss])
    
    def create_map(self, species='co2_ff', split=False, sector=None, adjust=-1, verbose=True, updates_perc=1):
        '''
        Convert a 1d array of emissions with a source dimension to a 2d array with 
        latitude and longitude dimensions

        Args
            species (str, optional)
                species to create a map for, defaults to 'co2_ff'
            split (bool, optional)
                if True the emissions are split into separate sources
            sector (str, optional)
                required sector name e.g. 'A_PublicPower'
                should correspond to row names in self.timeprofiles dataframes
                if None all sectors are returned
            adjust (int, optional)
                adjustment to account for 0 or 1 indexing of longitude_index and latitude_index
                i.e. if longitude_index and latitude_index are 1 indexed adjust = -1
            verbose (bool, optional)
                whether to print updates on progress
            updates_perc (int, optional)
                if verbose==True, this is the percentage completion at which updates will be printed

        Returns
            dict
                keys will be the source indices if split==True or 'total' if not
                values will be emission maps
        '''
        if sector is not None:
            print(f'Creating emissions map for {sector}')
        elif split:
            print('Creating separate emissions maps for each sector')
        else:
            print('Creating total emissions map')
        sectors_ind   = {str(sec).split("'")[1] : ss
                         for ss, sec in enumerate(self.emissions_ds.emis_cat_name.values)}
        index_req     = sectors_ind[sector] if sector else None
        
        sectors_index = list(set(self.emissions_ds.emission_category_index.values)) \
                        if sector is None else [index_req]
        map_shape     = (len(self.emissions_ds.latitude.values),
                         len(self.emissions_ds.longitude.values))
        
        # create a 2d map with lat & lon coords to fill with emissions values
        emissions_map = {'tot'     : np.zeros(map_shape)} if not split and sector is None else \
                        {index_req : np.zeros(map_shape)} if sector is not None else \
                        {sec       : np.zeros(map_shape) for sec in sectors_index}
        
        emissions_all = getattr(self.emissions_ds, species).values
        perc          = int(updates_perc * 0.01 * len(emissions_all))
        
        fails = []
        # loop through emissions sources to add to map
        for ee, emission in enumerate(emissions_all):
            # get the index of the lat & lon for the source
            ind_lat = int(self.emissions_ds.latitude_index[ee].values)  + adjust
            ind_lon = int(self.emissions_ds.longitude_index[ee].values) + adjust
            ind_sec = int(self.emissions_ds.emission_category_index.values[ee]) + adjust \
                      if split or sector is not None else 'tot'

            try:
                # add emissions to map
                if sector is None or ind_sec==index_req:
                    emissions_map[ind_sec][ind_lat][ind_lon] = emissions_map[ind_sec][ind_lat][ind_lon] + float(emission)

            except:
                fails.append(ee)

            if verbose and ee%perc==0:
                perc_done = int((ee+1) / len(emissions_all) * 100)
                print(perc_done, 'per cent completed')
            else: pass

        if verbose: print(len(fails), 'fails')
        else: pass

        self.emissions_map = emissions_map

    def create_map_file(self, species='co2_ff', split=False, sector=None, adjust=-1, verbose=True, updates_perc=1,
                        filename=None, outdir=None):
        self.create_map(species      = species,
                        split        = split,
                        sector       = sector,
                        adjust       = adjust,
                        verbose      = verbose,
                        updates_perc = updates_perc)
        
        sector_names = {sec : str(self.emissions_ds.emis_cat_name[sec-1].values).split("'")[1]
                              if sec!='tot' else sec
                        for sec in self.emissions_map.keys()}
        
        attributes   = getattr(self.emissions_ds, species).attrs
        map_dict     = {sector_names[sec]: (['lat','lon'], em_map, attributes)
                        for sec, em_map in self.emissions_map.items()}
        
        map_ds       = xr.Dataset(data_vars = map_dict,
                                  coords    = {'lat'  : self.emissions_ds.latitude.values,
                                               'lon'  : self.emissions_ds.longitude.values},
                                  attrs     = self.emissions_ds.attrs)
        
        outdir   = self.path if outdir is None else outdir
        
        sector   = '_total' if sector is None and not split else '' if sector is None else f'_{sector}'
        filename = filename if filename is not None else \
                   os.path.join(outdir, f'{species}{sector}_EUROPE_{self.year}.nc')
        
        print(f'Saving map to {filename}')
        
        map_ds.to_netcdf(filename)