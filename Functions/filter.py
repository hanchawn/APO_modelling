import json
import pandas as pd

from acrg.config.paths import Paths

with open(Paths.acrg / "data/site_info.json") as f:
    site_info=json.load(f,object_pairs_hook=dict)

def local_solar_time(time, sitelon):
    """
    Returns hour of day as a function of local solar time
    relative to the Greenwich Meridian. 
    """
    # convert lon to [-180,180], so time offset is negative west of 0 degrees
    if sitelon > 180:
        sitelon = sitelon - 360.
    time = time + pd.Timedelta(minutes=float(24*60*sitelon/360.))
    hours = time.to_pandas().index.hour
    return hours

def filter_time(dataset, site, hour_start, hour_end, keep_missing=False):
    """
    Filter the data for time during the given range

    Args
        dataset: xarray.DataArray or xarray.Dataset
            dataset to be filtered
        sitelon: int or float
            longitude at the site
        hour_start, hour_end: int
            start and end hour between which to keep the data
            time of the 24 hr clock
            e.g. for daytime: hour_start=11, hour_end=15
        keep_missing: bool
            whether to keep the filtered times as nans
    
    Returns
        xarray.DataArray or xarray.Dataset
    """
    sitelon = site_info[site][list(site_info[site].keys())[0]]['longitude']

    hours = local_solar_time(time=dataset.time, sitelon=sitelon)
    ti = [i for i, h in enumerate(hours) if h>=hour_start and h<=hour_end]
        
    if keep_missing:
        dataset_temp = dataset[dict(time = ti)]   
        dataset_out = dataset_temp.reindex_like(dataset)
        return dataset_out
    else:
        return dataset[dict(time = ti)]