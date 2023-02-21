import os, glob
import zipfile
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

from acrg.name import name
from acrg.convert import concentration

def make_date_string(year, months=None):
    '''
    Create a string to use for file naming from the input dates

    Inputs
    ------
    year: (int or float)
    months: (list or np.ndarray)

    Outputs
    -------
    string
        string combining the dates given
    '''
    year = int(year)
    months = [months] if type(months) in [int, float] else months
    date_str = f'{year}' if months==[mm for mm in range(1, 13)] or months is None else \
               f'{year}{str(months[0]).zfill(2)}-{str(months[-1]).zfill(2)}' if len(months)>2 else \
               f'{year}{str(months[0]).zfill(2)}-{year}{str(months[-1]).zfill(2)}' if len(months)==2 else \
               f'{year}{str(months[0]).zfill(2)}'
    
    return date_str

def get_timeseries(sites, year, sources, path=None, climatology=False, drop_monthly=True, verbose=True):
    '''
    '''
    path = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2') if path is None \
           else path
    
    ts_data_all = {}
    for site in sites:
        if verbose:
            print(f'Finding timeseries for {site}')
        # paths to data
        ts_path = os.path.join('/user', 'work', 'vf20487', 'Timeseries', 'o2_co2')

        # find the files
        clim_str = '_climatology_' if climatology else '_'
        ts_files = {source: glob.glob(os.path.join(ts_path, f'*{site}_{source}_timeseries{clim_str}{year}.nc'))
                    for source in sources}
        [print(f'{source}: {file_source}') for source, file_source in ts_files.items()]

        # open the files
        ts_data_dict = {source: name.open_ds(file_source[0]) for source, file_source in ts_files.items()
                        if len(file_source)>0}
        # create a new dataset with the time coords
        time = list(list(ts_data_dict.values())[0].values())[0].time.values
        ts_data = xr.Dataset(coords={'time': time})
        # loop through the timeseries, reformat the variable name, and add to the dataset
        for source, ts_data_source in ts_data_dict.items():
            for data_var in ts_data_source.data_vars:
                dv_name = f'{data_var.split("_")[0]}_ocean_{"_".join(data_var.split("_")[1:])}' if source=='ocean' \
                          else f'{data_var}_{source}'
                var_name = f'co2_{dv_name}' if not any([spec in dv_name for spec in ['co2', 'n2', 'o2', 'apo']]) else dv_name
                ts_data[var_name] = ts_data_source[data_var]

        # make sure the o2 ff sign is correct
        if 'o2_ff' in ts_data.data_vars:
                ts_data['o2_ff'] = -ts_data.o2_ff if (ts_data.o2_ff.fillna(1)>0).all() else ts_data.o2_ff
        # # make sur the units are correct
        for vv in ts_data.data_vars:
            ts_data[vv] = ts_data[vv] / concentration('ppm') if all(abs(ts_data[vv])<1e-4) else ts_data[vv]
            if '-' in vv:
                if verbose: print(f"renaming {vv}: {vv.replace('-', '_')}")
                ts_data = ts_data.rename({vv: vv.replace('-', '_')})

        # drop uncertainties, differences
        vars_drop = [dv for dv in ts_data.data_vars if any([any([ss in dv for ss in ['unc', 'diff']])])]
        vars_drop2 = [dv for dv in ts_data.data_vars if dv not in vars_drop and all([dv.split('_')[0]!=spec for spec in ['apo', 'co2', 'n2', 'o2']])]
        vars_drop3 = [dv for dv in ts_data.data_vars if 'mth' in dv] if drop_monthly else []
        vars_drop = vars_drop+vars_drop2+vars_drop3
        if verbose and len(vars_drop)>0:
            print(f'\nDropping variables: \n{vars_drop}')
        ts_data = ts_data.drop_vars(vars_drop)

        if drop_monthly:
            vars_rename = {dv: ''.join(dv.split('_day')) for dv in ts_data.data_vars if 'day' in dv}
            # if 'ocean' in sources and 'co2_ocean_nemo_mth' in ts_data.data_vars:
            #     vars_rename['co2_ocean_nemo_mth'] = 'co2_ocean_nemo'
            ts_data_all[site] = ts_data.rename(vars_rename)
        else:
            ts_data_all[site] = ts_data

        if verbose:
            print(f'\nData variables: \n{[dv for dv in ts_data_all[site].data_vars]}\n---------------------------\n')

    return ts_data_all

def get_prior_info(attributes, num_priors=None, keep_all=False):
    '''
    Get all of the attributes from different fields

    Inputs
    ------
    attributes: dict
        nested dictionary where the outer key is the sourc or sector, and the inner dict
        contains the attributes
    num_priors: int or float, optional
        the maximum number of priors used for any given field
        searches for f'prior_file_{num}' in the attrbutes dictionary
    
    Outputs
    -------
    dict
        the prior_info_dict needed for flux.write
        contains version, reoslution, and reference for each input fields
    '''
    attributes = {'all': attributes} if type(attributes[list(attributes.keys())[0]])!=dict else attributes

    num_priors = 10 if num_priors is None else num_priors
    prior_attrs = ['version', 'raw_resolution', 'reference']
    # get the prior info from each field attributes dict
    prior_info_dict = {source: {attr_source[f'prior_file_{num}']: [attr_source[f'prior_file_{num}_{attr}']
                                                                   for attr in prior_attrs]
                                for num in range(0, num_priors) if f'prior_file_{num}' in attr_source}
                       for source, attr_source in attributes.items()}
    # flatten dictionary, also removes repated entries
    join_keys = True if keep_all else False
    separator = '-' if keep_all else '+'
    prior_info_dict = flatten_nested_dict(prior_info_dict, separator=separator, join_keys=join_keys)

    return prior_info_dict

def flatten_nested_dict(nested_dict, separator=None, join_keys=True):
    '''
    Convert a nested dictionary into a flat dictionary

    i.e.
    {'a': {'one': 1}, 'b': {'two': 2}} --> {'a_one': 1, 'b_two': 2}

    Inputs
    ------
    nested_dict: dict
        nested dictionary to be flattened
    separator: str
        separator to use when joining the inner and outer keys
        defaults to '_'
    join_keys: bool
        if True then the keys from the inner and outer dictionary are joined together
        if False then the keys from the inner dictionary are used
    
    Outputs
    -------
    dict
        flat dictionary where the keys are the keys of the inner and outer dictionaries joined
    '''
    separator = '_' if separator is None else separator
    df = pd.json_normalize(nested_dict, sep=separator)
    flat_dict = df.to_dict(orient="records")[0]
    
    if not join_keys:
        flat_dict = {separator.join(key.split(separator)[1:]): value for key, value in flat_dict.items()}

    return flat_dict

def join_monthly_files(file_pattern:str, year:int, path:str=None, zip_months:bool=False,
                       remove_months:bool=False, verbose:bool=True, join:str='_'):
    '''
    Merge monthly files into a file for the whole year

    Inputs
    ------
    file_pattern: str
        start of filename, i.e. filename = f'{file_pattern}_{year}{month}.nc
        where months are ['01', '02', '03', '04', '05', '06',
                          '07', '08', '09', '10', '11', '12']
    year: int or str
        year associated with files, should be part of the filename
    path: str or PosixPath
        path to files
    zip_months: bool
        if True, monthly files will be zipped and unzipped files removed
        default False
    remove_months: bool
        if True, monthly files will be removed
        default False
    verbose: bool
        default True
    
    Outputs
    -------
    Saves a netcdf file with the merged data from the monthly files
    '''
    # change working dir
    if path is not None:
        cwd = os.getcwd()
        os.chdir(path)
    
    filename = os.path.join(path, f'{file_pattern}{join}{year}.nc')
    if os.path.isfile(filename):
        os.remove(filename)

    if verbose:
        print(f'Searching {os.getcwd()}')
        print(f'  for files with pattern: {file_pattern}{join}{year}*.nc')
    
    # find the file for each month
    files = glob.glob(os.path.join(path, f'{file_pattern}{join}{year}*.nc')) if path is not None else \
            glob.glob(f'{file_pattern}{join}{year}*.nc')
    
    if verbose:
        print(f'Monthly files found:')
        [print(f' --- {ff}') for ff in sorted(files)]

    # check that there are files for all months
    if len(files)<12:
        raise OSError('Could not find files for all months')
    
    # extract monthly data and merge all months
    data = [xr.open_dataset(ff) for ff in files]
    # data = [open_ds(ff) for ff in files]
    data_all = xr.combine_by_coords(data)

    # save the merged data
    if verbose: print(f'Saving to: {filename}')
    data_all.to_netcdf(filename)
    
    # zip the monthly files
    if zip_months:
        if verbose: print('Zipping monthly files')
        with zipfile.ZipFile(f'{file_pattern}{join}{year}.zip', 'w') as zipF:
            for ff in files:
                zipF.write(ff, compress_type=zipfile.ZIP_DEFLATED)
    
    # remove the monthly files
    if zip_months or remove_months:
        if verbose: print('Removing monthly files')
        [os.remove(ff) for ff in files]
    
    # go back to initial working dir
    if path is not None:
        os.chdir(cwd)

    if verbose: print('Done')

    return None
            

def show_map(data, ax=None, coast_color='black', lat=None, lon=None, crop_uk=False,
             im_kwargs={}, colorbar=True, colorbar_kwargs=None, colorbar_label_kwargs=None,
             colorbar_ticks=None, colorbar_tick_kwargs=None, fig_kwargs=None, ax_kwargs=None,
             show=True, reverse_cmap=False, return_im=False):
    '''
    Inputs
    ------
    data
        data to show, should be 2d with lat and lon dimensions
    ax: matplotlib axis, optional
        axis to show image in
    lat, lon: slice, optional
        lat and lon slice to crop to
        e.g. lat = slice(45, 63)
             lon = slice(-15, 8)
    crop_uk: bool, optional
        crop to UK if True
    coast_color: string, optional
        color of coastlines, default = 'black'
    im_kwargs: dict, optional
        arguments for formatting the image
        can include any argument for plt.contourf
        default = {'cmap': 'inferno'}
    colorbar: bool
        if True add a colorbar, defaults to horizontal with 0.05 pad
    colorbar_kwargs: dict, optional
        arguments for the colorbar
        can include any argument for plt.colorbar
        default = {'orientation': 'horizontal', 'pad': 0.05}
    colorbar_label_kwargs: dict, optional
        arguments for the colorbar label
        can include any argument for colorbar.set_label
        e.g. {'label': 'label', 'fontsize': 14}
    
    fig_kwargs: dict, optional
        arguments for creating the figure
        default = {'figsize': (10, 10)}
    ax_kwargs: dict, optional
        arguments for creating the axis
        default = {'projection': ccrs.PlateCarree()}
    
    '''
    # define lat and lon for UK if crop_uk=True
    lat = slice(45, 63) if crop_uk else lat
    lon = slice(-15, 8) if crop_uk else lon

    data = [data] if type(data)!=list else data
    # if not given create an axis to show map in
    if ax is None:
        fig_kwargs = {} if fig_kwargs is None else fig_kwargs
        ax_kwargs = {} if ax_kwargs is None else ax_kwargs
        ax_kwargs['projection'] = ax_kwargs['projection'] if 'projection' in ax_kwargs else ccrs.PlateCarree()
        fig_kwargs['figsize'] = fig_kwargs['figsize']if 'figsize' in fig_kwargs else (10*len(data), 10)
        fig_kwargs['nrows'] = fig_kwargs['nrows'] if 'nrows' in fig_kwargs else 1
        fig_kwargs['ncols'] = fig_kwargs['ncols'] if 'ncols' in fig_kwargs else len(data)
        fig, axes = plt.subplots(subplot_kw=ax_kwargs, **fig_kwargs)
        axes = np.array([axes]) if type(axes)!=np.ndarray else axes
    else:
        axes = np.array([ax]) if type(ax)!=np.ndarray else ax
    
    # crop the data if required
    data = [dat.loc[dict(lat=lat, lon=lon)] for dat in data] if all([lat, lon]) else data
    # show the data
    im_kwargs = [im_kwargs] if type(im_kwargs)==dict else im_kwargs
    for im_kwarg in im_kwargs:
        im_kwarg['cmap'] = im_kwarg['cmap'] if 'cmap' in im_kwarg else 'inferno'
        im_kwarg['cmap'] = plt.cm.get_cmap(im_kwarg['cmap']).reversed() if reverse_cmap else im_kwarg['cmap']
        im_kwarg['levels'] = im_kwarg['levels'] if 'levels' in im_kwarg else 60
    im_kwargs = im_kwargs*len(data) if len(im_kwargs)==1 else im_kwargs
    ims = [axes.flatten()[dd].contourf(dat.lon,
                                       dat.lat,
                                       dat,
                                       transform = ccrs.PlateCarree(),
                                       **im_kwargs[dd])
           for dd, dat in enumerate(data)]

    # add coastlines
    coast_color = [coast_color]*len(data) if type(coast_color)==str else coast_color*len(data) if len(coast_color)== 1 else coast_color
    coast = [ax.coastlines(resolution='50m', color=coast_color[aa], linewidth=1)
             for aa, ax in enumerate(axes.flatten())]

    # add colorbar
    if colorbar:
        colorbar_kwargs = {} if colorbar_kwargs is None else colorbar_kwargs
        if 'orientation' not in colorbar_kwargs:
            colorbar_kwargs['orientation'] = 'horizontal'
        if 'pad' not in colorbar_kwargs:
            colorbar_kwargs['pad'] = 0.05
        if 'label' in colorbar_kwargs:
            colorbar_label_kwargs = {} if colorbar_label_kwargs is None else colorbar_label_kwargs
            colorbar_label_kwargs['label'] = colorbar_kwargs['label']
            colorbar_kwargs.pop('label')
        colorbar_ticks = colorbar_ticks if type(colorbar_ticks) in [list, np.ndarray] or colorbar_ticks==None else [colorbar_ticks]
        cbars = [plt.colorbar(ims[aa], ax=ax, **colorbar_kwargs) for aa, ax in enumerate(axes.flatten())] if colorbar_ticks is None else \
                [plt.colorbar(ims[aa], ax=ax, ticks=colorbar_ticks[aa], **colorbar_kwargs) for aa, ax in enumerate(axes.flatten())]
        if colorbar_label_kwargs is not None:
            for cc, cbar in enumerate(cbars):
                # create a new dict for the colorbar label kwargs for this axis
                label_kwargs = {}
                for key, arg in colorbar_label_kwargs.items():
                    # convert the arg to alist so that it's length can be checked
                    arg = [arg] if type(arg)!=list else arg
                    if len(arg) not in [1, len(data)]:
                        raise ValueError(f'The length of colorbar_label_kwargs["{key}"] should be 1 or len(data)')
                    label_kwargs[key] = arg[0] if len(arg)==1 else colorbar_label_kwargs[key][cc]
                cbar.set_label(**label_kwargs)
        if colorbar_tick_kwargs is not None:
            if 'fontsize' in colorbar_tick_kwargs.keys():
                colorbar_tick_kwargs['labelsize'] = colorbar_tick_kwargs['fontsize']
                colorbar_tick_kwargs.pop('fontsize', None)
            [cbar.ax.tick_params(**colorbar_tick_kwargs) for cbar in cbars]
    
    if show:
        plt.show()
    
    if ax is None:
        if return_im:
            return fig, axes, ims
        else:
            return fig, axes
    else:
        if return_im:
            return axes, ims
        else:
            return axes


def convert_latlon_metres(lat1, lon1, lat2, lon2, units=None, verbose=True):
    '''
    Find the distance between two points in latitude and longitude
    
    Inputs
    ------
    lat1, lon1, lat2, lon2: int or float
        latitude and longitude of the two points
    units: str, optional
        units to return the distance in
        should be 'm', 'metres', 'meters', 'km', 'kilometres', 'kilometers'
    '''
    units = 'km' if units is None else units
    if verbose:
        print(f'Converting to units of {units}')
    
    if units not in ['m', 'metres', 'meters', 'km', 'kilometres', 'kilometers']:
        print('Please select units from: [m, metres, meters, km, kilometres, kilometers]')
        return None
    
    radius_earth = 6378.137    # radius of the earth, km
    radius_earth = radius_earth * 1e3 if units in ['m', 'metres', 'meters'] else radius_earth
         
    # angle distance between latitudes and longitudes, radians
    dLat = np.deg2rad(lat2 - lat1)
    dLon = np.deg2rad(lon2 - lon1)
 
    # convert to radians
    lat1 = np.deg2rad(lat1)
    lat2 = np.deg2rad(lat2)
 
    # apply haversine formulae
    a = ((np.sin(dLat / 2))**2 + (np.sin(dLon / 2))**2 * np.cos(lat1) * np.cos(lat2))
    
    c = 2 * np.arcsinh(np.sqrt(a))
    return radius_earth * c

def latlon_distance(lat, lon, units=None, verbose=True):
    '''
    Get the distance between 2 points in km
    
    Inputs
    ------
    lat, lon: list
        latitude and longitude of the 2 points, degrees
    units: str, optional
        units to return the distance in, defaults to km
        should be 'm', 'metres', 'meters', 'km', 'kilometres', 'kilometers'
    verbose: bool, optional
        default: True
    '''
    units = 'km' if units is None else units
    if verbose:
        print(f'Converting to units of {units}')
    
    if units not in ['m', 'metres', 'meters', 'km', 'kilometres', 'kilometers']:
        print('Please select units from: [m, metres, meters, km, kilometres, kilometers]')
        return None
    
    latlon = {key: [np.deg2rad(ll) for ll in latlon_deg]
              for key, latlon_deg in {'lat': lat, 'lon': lon}.items()}
    d_latlon = {key: abs(ll[1] - ll[0]) for key, ll in latlon.items()}

    radius_earth = 6378.137    # radius of the earth, km
    radius_earth = radius_earth * 1e3 if units in ['m', 'metres', 'meters'] else radius_earth

    # apply haversine formulae
    a = (np.sin(d_latlon['lat'] / 2))**2
    b = (np.sin(d_latlon['lon'] / 2))**2 * np.cos(latlon['lat'][0]) * np.cos(latlon['lat'][1])
    c = 2 * np.arcsinh(np.sqrt(a+b))
    return radius_earth * c

def convert_degrees_metres(lat, d_lat, d_lon, units=None, verbose=True):
    '''
    Convert the size of a pixel from degrees to km
    
    Inputs
    ------
    lat: list or int or float
        latitude of the centre of the pixel, degrees
    d_lat, d_lon: int or float
        size of the pixel, degrees
    units: str, optional
        units to return the distance in, defaults to km
        should be 'm', 'metres', 'meters', 'km', 'kilometres', 'kilometers'
    verbose: bool, optional
        default: True
    '''
    lon = [0, d_lon]
    lat = [lat - d_lat/2, lat + d_lat/2]
    dist = latlon_distance(lat, lon, units=units, verbose=verbose)
    return dist

def combine_diff_resolution(data_1, data_2, method='add', data_1_res=None, data_2_res=None, verbose=True):
    '''
    Combine datasets which have different time resolutions
    Avoids having to resample data to the same resolution in order to add,
    subtract, multiply, or divide them
    
    Args
        data_1, data_2 (xarray.DataArray or xarray.Dataset)
            Data to combine
            Must have resolution of 1 hour, 1 day, 1 month, or 1 year
        method (str)
            Method by which to combine the datasets
            Can be:
                'add':      data_1 + data_2 (default)
                'subtract': data_1 - data_2
                'multiply': data_1 * data_2
                'divide':   data_1 / data_2
    
    Returns
        xarray.DataArray or xarray.Dataset
            has the same shape as resolved_data
    '''
    if 'np' not in dir(): import numpy as np
    data = {dd: dat for dd, dat in enumerate([data_1, data_2])}
    # calculate the time step for each dataset
    time_step = {0: data_1_res, 1: data_2_res}
    time_step = {dd: (dat.time[1] - dat.time[0]).values 
                     if time_step[dd] is None else
                     np.timedelta64(int(time_step[dd][0]), time_step[dd][1]).astype('timedelta64[ns]')
                 for dd, dat in data.items()}
    
    if time_step[0]==time_step[1]:
        # if the time steps are equal then apply the method as normal
        if method.lower()=='multiply':
            return data_1 * data_2
        elif method.lower()=='add':
            return data_1 + data_2
        elif method.lower()=='subtract':
            return data_1 - data_2
        elif method.lower()=='divide':
            return data_1 / data_2
        else:
            print(f'Method not recognised: {method}')
            return(None)
    
    # work out which dataset is resolved and which is integrated
    data_val = {'resolved'  : [dd for dd, step in time_step.items() if step==np.nanmin(np.array(list(time_step.values())))][0],
                'integrated': [dd for dd, step in time_step.items() if step==np.nanmax(np.array(list(time_step.values())))][0]}
    data = {res: data[dd] for res, dd in data_val.items()}
    
    # work out the time scale for each dataset
    time_scale = {res: 'Month' if time_step[dd].astype('timedelta64[M]').astype(int)>0 else 
                       'Dayofyear' if time_step[dd].astype('timedelta64[D]').astype(int)>0 else 'hour'
                  for res, dd in data_val.items()}
    time_step  = {res: time_step[dd].astype(f'timedelta64[{time_scale[res][0]}]').astype(int)
                  for res, dd in data_val.items()}
    
    if verbose:
        print(f'Method to combine datasets: {method}',
              f'\nIntegrated data has time step of {time_step["integrated"]} {time_scale["integrated"].lower()}')
    
    # function to apply to the resolved data
    def combine_method(resolved, method):
        # select the correct time slice from the integrated data
        integrated = data['integrated'].sel(time=resolved.time[0])
        
        if method.lower()=='multiply':
            return resolved * integrated
        elif method.lower()=='add':
            return resolved + integrated
        elif method.lower()=='subtract_high':
            return integrated - resolved
        elif method.lower()=='subtract_low':
            return resolved - integrated
        elif method.lower()=='divide_low':
            # divide by the low res data
            return resolved / integrated
        elif method.lower()=='divide_high':
            return integrated / resolved
        else:
            print(f'Method not recognised: {method}')
            return(None)
    
    # group the data by the same time scale as the integrated da
    grouped_data = data['resolved'].groupby(f'time.{time_scale["integrated"].lower()}')
    
    # the required method depends on which order the data was given in
    method = method if method in ['multiply', 'add'] else \
             '_'.join([method, 'high']) if (data['integrated']==data_1).all() else \
             '_'.join([method, 'low']) if (data['integrated']==data_2).all() else \
             method
    
    # apply the combine method to each slice of the grouped data
    output = grouped_data.map(combine_method, method=method)
    
    return output

def multiply_2d_3d(array_2d, array_3d, output='array',
                   lon_var='lon', lat_var='lat', time_var='time'):
    '''
    Multiply a 2D array by each slice of a 3D array
    Multiplies along the 3rd dimension of array_3d
    
    Args:
        array_2d, array_3d (numpy.ndarray or xarray.DataArray)
            arrays with 2 and 3 dimensions respectively
        output (str)
            output type, either 'array', 'dataarray', or 'xarray'
            defaults to 'array'
            (case insensitive)
    
    Returns:
        numpy.ndarray
            array with the same shape as array_3d
    '''
    if 'np' not in dir(): import numpy as np
    
    if len(array_2d.shape)>2:
        print(f'array_2d must have 2 dimensions')
        return None
    
    # multiply the 2d array by each slice of the 3d array
    multiplied_array = np.einsum('ij,ijk->ijk', array_2d, array_3d)
    
    # if output is dataarray then convert the array to an xarray DataArray
    if output.lower() in ['dataarray', 'xarray']:
        import xarray as xr
        multiplied_ds = xr.DataArray(data = multiplied_array,
                                     dims=["lat", "lon", "time"],
                                     coords=dict(lon  = getattr(array_3d, lon_var),
                                                 lat  = getattr(array_3d, lat_var),
                                                 time = getattr(array_3d, time_var)))
        return multiplied_ds
    else:
        return multiplied_array


def round_sig_fig(x, n):
    if 'np' not in dir(): import numpy as np
    xr = (np.floor(np.log10(np.abs(x)))).astype(int)
    xr=10.**xr*np.around(x/10.**xr,n-1)   
    return xr 


def centre_timestamp(dataset):
    if 'pd' not in dir() : import pandas as pd
    if 'np' not in dir() : import numpy as np
    if len(dataset.time) == 12:
        time_pd = pd.date_range(start = "%s-01-01" %int(dataset.time[0]['time.year'].values),
                                periods = len(dataset.time), freq = 'MS')
        time    = [np.datetime64(t) for t in time_pd]
    else:
        diff    = (dataset.time[1].values - dataset.time[0].values)/2
        time    = dataset.time.values - pd.Timedelta(diff)
    return time


def hour2month(dataset, timeframe='month'):
    if 'pd' not in dir() : import pandas as pd
    if 'dt' not in dir() : import datetime as dt
    month_data = dataset.groupby('time.'+timeframe)
    monthly    = month_data.mean('time')

    times    = dataset.time.values
    timelist = pd.date_range(start=pd.to_datetime(times[0]), end = pd.to_datetime(times[-1]), freq='1M')
    new_time = [dt.datetime(i.year, i.month, 1) for i in timelist]
    monthly  = monthly.update({'month' : new_time})
    monthly  = monthly.rename({'month' : 'time'})
    monthly  = monthly.transpose('latitude','longitude','time')
    
    return monthly


def data_split_hours(data, year, scaling_factors, time_space=2):
    if 'np' not in dir() : import numpy as np
    date_range = [366 if year%4==0 else 365][0]   # check if year is a leap year
    split_data = np.zeros([data.shape[0],
                           data.shape[1],
                           12*date_range])
    num = 0
    for day in range(date_range):
        for hour in np.arange(0, 24, time_space):
            # update the split_data array with the scaled data for each hour
            split_data[:,:,num] = data * scaling_factors[hour].values
            num+=1
    
    return split_data


def hour_split(data, years, scaling_factors, time_space=2, max_year=2018,
               species_obs_type=None, verbose=False):
    '''
    Split data into intervals
    
    Parameters
    ----------
    data: dict
    years: list
    scaling: factors
    time_space: int or float
        time between intervals, defaults to 2
    max_year: int or float
        the most recent year for which there is data
    '''
    if 'np' not in dir() : import numpy as np
    split_data = {}
    
    if species_obs_type and verbose:
        print ('Working on', species_obs_type)
    for year in years:
        flag = 0
        if verbose: print (year)
        if year in list(data.keys()):
            split_data[year] = data_split_hours(data[year], year, scaling_factors)
#             date_range = [366 if year%4==0 else 365][0]   # check if year is a leap year
#             # make array of zeros for the scaled data for every hour in the year
#             split_data[year] = np.zeros([data[year].shape[0],
#                                          data[year].shape[1],
#                                          12*date_range])
#             num = 0
#             for day in range(date_range):
#                 for hour in np.arange(0, 24, time_space):
#                     # update the split_data array with the scaled data for each hour
#                     split_data[year][:,:,num] = data[year] * scaling_factors[hour].values
#                     num+=1
        elif year>max_year:
            if year%4==0:
                leap_yr_day_hrs = [int(24./float(time_space) * 59), 
                                   int(24./float(time_space) * 60)]
                ds_flux_1      = split_data[max_year][:,:,0:leap_yr_day_hrs[1]]
                ds_flux_2      = split_data[max_year][:,:,leap_yr_day_hrs[1]:split_data[max_year].shape[-1]]
                ds_flux_insert = split_data[max_year][:,:,leap_yr_day_hrs[0]:leap_yr_day_hrs[1]]
                ds_flux_new    = np.append(ds_flux_1, ds_flux_insert, axis=2)
                ds_flux_new    = np.append(ds_flux_new, ds_flux_2, axis=2)
                split_data[year] = ds_flux_new
            else:
                split_data[year] = split_data[max_year]
    return split_data


def regrid_mask(data, latitude, longitude, latitude_out, longitude_out):
    '''
    Regrid data and mask pixels from EDGAR grid where there is no information from the NAEI grid
    
    Parameters
    ----------
    data : dict
        data for each year
    latitude, longitude : dict
        lat and long associated with the data for each year
    latitude_out, longitdue_out : array
        required output lat and long
    
    Result
    ------
    flux_regrid : dict
        regridded data for each year
    '''
    if 'np' not in dir() : import numpy as np
    from acrg_grid.regrid import regrid2d
    years       = data.keys()
    flux_regrid = {year: regrid2d(data[year],
                                  latitude[year], longitude[year],
                                  latitude_out, longitude_out)
                   for year in years}
    
    # Setting the fill value for pixels from EDGAR grid where there is no information from the NAEI grid
    [np.ma.set_fill_value(flux_regrid[year], 0) for year in years]
    
    # Grab the regridded data, where values that were masked are set to 0.
    flux_regrid = {year: flux_regrid[year][0].data for year in years}
    return flux_regrid


def MHD_baseline(start_date = "2019-01-01", end_date = "2019-06-30", show_plot=False,
                 obs_directory=None):
    '''
    Read in baseline data. The baseline data are daily baseline values of the gas calculate by Alistair Manning.
    These files are currently on Snowy: /dagage2/agage/metoffice/baselines_MHD/daily/
    Code then interpolates based on these values to the time step of requested
    The start and end dates are important. And it assumes an hourly time step.     
    '''
    if 'os' not in dir() : import os
    if 'pd'  not in dir() : import pandas as pd
    if 'plt' not in dir() and show_plot : import matplotlib.pyplot as plt
    if 'interpolate' not in dir() : from scipy import interpolate
    
    baselines             = pd.read_csv(os.path.join(obs_directory, 'MH_G_co2_day.txt'),
                                        delim_whitespace=True)
    baselines['datetime'] = pd.to_datetime(baselines['Yr']*10000+baselines['Mt']*100+baselines['Dy'],
                                           format='%Y%m%d')
    baselines.set_index('datetime')
    timeseries     = baselines['datetime']
    timeseries     = (ts - np.datetime64(start_date)) / np.timedelta64(1, 'h')
    timeseries_all = pd.date_range(start_date, end_date, freq='H')
    timeseries_all = (ts_all - np.datetime64(start_date)) / np.timedelta64(1, 'h')
    
    y              = baselines['Conc']
    
    # Extends the baseline dataset to equal the last value in the dataset until the end of the requested period.
    # Otherwise extrapolation can go in unexpected directions.
    if timeseries.iloc[-1] < ts_all[-1]:
        y_add          = pd.Series([y.iloc[-1]])
        timeseries_add = pd.Series([ts_all[-1]])
        y              = y.append(y_add, ignore_index=True)
        timeseries     = timeseries.append(timeseries_add, ignore_index=True)
        
    f = interpolate.PchipInterpolator(timeseries, y, extrapolate=True)
    # Get the times from the dataset which will have the same times as the modelled ch4 enhancements
    MHD_co2_baseline = f(timeseries_all)
    if show_plot:
        fig = plt.figure()
        fig.plot(timeseries, y)
        fig.plot(timeseries_all, MHD_co2_baseline)
        plt.show()
    return (f)