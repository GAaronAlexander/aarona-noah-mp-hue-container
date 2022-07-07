import s3fs
import xarray as xr
import xesmf
import numpy as np
import pandas as pd
import os

"""
Import your specific initial conditions: 
Inputs: 
Geogrid file - created from WPS
Beginning date - YYYYMMDDHH
Save Location - a path for where to save data

You must have set up an era5 account to be able to use this script, as 
we are downloading a few data-files to create this, and it requires you to have a 
copernicous account with appropriate set-up. 
"""


def get_region_era5(era5,geo):
    
    """
    Takes an era5 file, and a geogrid file, and then calculates the start and end points 
    of the era5 file based on the geogrid file. Adds 3 to each edge of the ERA5 data for 
    a safety protocol
    """
    
    ## determine the nearest coordinates of the geogrid file, and then add 3
    lonsv, latsv = np.meshgrid((era5.lon.values+180)%360-180,era5.lat.values)
    ## this is ERA5 specific!
    
    ## get the id of the lons from geogrid
    LON_start = np.abs(lonsv - geo.lon[0,0].values)
    LON_end = np.abs(lonsv - geo.lon[-1,-1].values)

    idx_lonstart = np.where(LON_start == LON_start.min())
    idx_lonend = np.where(LON_end == LON_end.min())
    input_lonstart = idx_lonstart[1][0]-3
    input_lonend = idx_lonend[1][0]+3


    ## get the id of the lats from geogrid
    LAT_start = np.abs(latsv -geo.lat[0,0].values)
    LAT_end = np.abs(latsv - geo.lat[-1,-1].values)

    idx_latstart = np.where(LAT_start == LAT_start.min())
    idx_latend = np.where(LAT_end == LAT_end.min())
    input_latstart = idx_latstart[0][0]+3
    input_latend = idx_latend[0][0]-3
    
    return(input_lonstart,input_lonend,input_latend,input_latstart) # order due to the way era5 defines coordinates

def transform_era5_to_dataarray_2d(dsubset):
    
    """
    takes a zarr dataarray, and creates a dataset that has the correct lats and lons to 
    feed into HRLDAS
    """
    
    lonsv, latsv = np.meshgrid((dsubset.lon.values+180)%360-180,dsubset.lat.values) 
    # meshgrid is based on the notation of ERA5 data! !!SHOULD BE CHANGED IF USING ANOTHER analysis forcing

    data_array_input = xr.Dataset(
        {
        name_variables[0]:(["y","x"],dsubset[0,1:-1,1:-1].values),
        }
        ,coords={
            "lat": (["y","x"],latsv[1:-1,1:-1]),
            "lon": (["y","x"], lonsv[1:-1,1:-1]),
            "lat_b": (["y_b","x_b"], 0.25*(latsv[:-1,:-1]+latsv[1:,1:]+latsv[1:,:-1]+latsv[:-1,1:])),
            "lon_b": (["y_b","x_b"], 0.25*(lonsv[:-1,:-1]+lonsv[1:,1:]+lonsv[1:,:-1]+lonsv[:-1,1:])),
            
            
        }
    )   
    return(data_array_input)

def transform_era5_to_dataarray_3d(dsubset3d,name):
    
    """
    takes a zarr dataarray, and creates a dataset that has the correct lats and lons to 
    feed into HRLDAS
    """
    
    lonsv, latsv = np.meshgrid((dsubset3d.lon.values+180)%360-180,dsubset3d.lat.values) 
    # meshgrid is based on the notation of ERA5 data! !!SHOULD BE CHANGED IF USING ANOTHER analysis forcing

    data_array_input = xr.Dataset(
        {
        name:(["time","y","x"],dsubset3d[:,1:-1,1:-1].values),
        }
        ,coords={
            "lat": (["y","x"],latsv[1:-1,1:-1]),
            "lon": (["y","x"], lonsv[1:-1,1:-1]),
            "lat_b": (["y_b","x_b"], 0.25*(latsv[:-1,:-1]+latsv[1:,1:]+latsv[1:,:-1]+latsv[:-1,1:])),
            "lon_b": (["y_b","x_b"], 0.25*(lonsv[:-1,:-1]+lonsv[1:,1:]+lonsv[1:,:-1]+lonsv[:-1,1:])),
            
            
        }
    )   
    return(data_array_input)


def get_regridder(grid_input,grid_out,option="bilinear"):
    """
    this function takes the two grids, and gets a regridder that we can use to get all gridded data!
    depends on xemsf.
    Option should be the method you want, default is "bilinear"
    """
    regridder_weights=xesmf.Regridder(grid_input, grid_out, option)
    return(regridder_weights)
    
