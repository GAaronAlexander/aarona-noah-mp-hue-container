import s3fs
import xarray as xr
import xesmf
import numpy as np
import pandas as pd
import cdsapi
import os

#### IMPORTANT 
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
        '2mTemp':(["y","x"],dsubset[0,1:-1,1:-1].values),
        }
        ,coords={
            "lat": (["y","x"],latsv[1:-1,1:-1]),
            "lon": (["y","x"], lonsv[1:-1,1:-1]),
            "lat_b": (["y_b","x_b"], 0.25*(latsv[:-1,:-1]+latsv[1:,1:]+latsv[1:,:-1]+latsv[:-1,1:])),
            "lon_b": (["y_b","x_b"], 0.25*(lonsv[:-1,:-1]+lonsv[1:,1:]+lonsv[1:,:-1]+lonsv[:-1,1:])),
            
            
        }
    )   
    return(data_array_input)

def transform_era5_to_dataarray_2d_netcdf(dsubset,name):
    
    """
    takes a zarr dataarray, and creates a dataset that has the correct lats and lons to 
    feed into HRLDAS
    """
    
    lonsv, latsv = np.meshgrid((dsubset.longitude.values+180)%360-180,dsubset.latitude.values) 
    # meshgrid is based on the notation of ERA5 data! !!SHOULD BE CHANGED IF USING ANOTHER analysis forcing

    data_array_input = xr.Dataset(
        {
        name:(["y","x"],dsubset[0,1:-1,1:-1].values),
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
    
    lonsv, latsv = np.meshgrid((dsubset3d.longitude.values+180)%360-180,dsubset3d.latitude.values) 
    # meshgrid is based on the notation of ERA5 data! !!SHOULD BE CHANGED IF USING ANOTHER analysis forcing

    data_array_input = xr.Dataset(
        {
        name:(["soil_layers_stag","time","y","x"],dsubset3d[:,:,1:-1,1:-1].values),
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
    

def get_era5_initial_conditions(start_date):
    
    """
    Takes a start-date, queries to get a new netcdf file and loads it
    and then returns a dataset that has been loaded into an xarray dataset
    DATE must be in the format YYYYMMDDHH
    """
    
    c = cdsapi.Client()

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'year': f'{start_date[0:4]}',
            'month': f'{start_date[4:6]}',
            'day': f'{start_date[6:8]}',
            'time': f'{start_date[8:]}:00',
            'variable': [
                'skin_reservoir_content', 'skin_temperature', 'snow_depth',
                'soil_temperature_level_1', 'soil_temperature_level_2', 'soil_temperature_level_3',
                'soil_temperature_level_4', 'volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2',
                'volumetric_soil_water_layer_3', 'volumetric_soil_water_layer_4',
            ],
        },
        'download.nc')

    dset = xr.open_dataset('./download.nc')

    return(dset)

def get_initial_LAI(geogrid,date_start,name):
    """
    generates a time series of LAI based on the day of year of our start date,
    LAI is linearly interpolated, as it is done in HRLDAS main driver
    """
    date = pd.to_datetime('2010040100',format='%Y%m%d%H')

    LAI = geogrid[name].values 
    LAI = np.append(l,l[:,0,:,:].reshape(1,1,l.shape[2],l.shape[3]),axis=1)

    months_index = pd.date_range(
                "1999-01-01",
                periods=13,
                freq=pd.DateOffset(months=1),
            )
    LAI_xr= xr.DataArray(data=l,dims={'Time':l.shape[0],'Months':months_index,'north_south':114,'east_west':99},coords={'Months':months_index})
    LAI_xr_daily = LAI_xr.resample(Months='1D').interpolate("linear")

    final_return = LAI_xr_daily.isel(Months=date.dayofyear)
    return(final_return)

def soil_coarse_stack(dset_era5,name):
    '''
    This generates a coarse domain (era5 grid) soil layer dataset. 
    It assumes that you are using the standard soil layer distribution within 
    NOAH-MP, which is 0-10cm, 10-40cm, 40-100cm, and 100-200cm.
    ERA5 stagger is Layer 1: 0 - 7cm, Layer 2: 7 - 28cm, Layer 3: 28 - 100cm, Layer 4: 100 - 289cm
    
    we just push in the 3rd and 4th layers, since there is no overlapping
    '''
    
    first_layer_era5 = dset_era5[name+'1']
    second_layer_era5 = dset_era5[name+'2']
    third_layer_era5 = dset_era5[name+'3']
    fourth_layer_era5 = dset_era5[name+'4']
    
    first_layer_noah = first_layer_era5*7/10 + second_layer_era5* 3/10 
    second_layer_noah = second_layer_era5*18/30 + third_layer_era5*12/30
    third_layer_noah = third_layer_era5
    fourth_layer_noah = fourth_layer_era5
    
    soil_data_return = xr.concat([first_layer_noah ,second_layer_noah ,third_layer_noah ,fourth_layer_noah],dim="soil_layers_stag") 
    return(soil_data_return)


############################### Begin typical code
def run(start_date,save_location,geo_file):
    
    #create output filename
    output_filename  = f'HRLDAS_setup_{start_date}_d1'
    
    ## 
    data_variables = [f'Times',f'XLAT_M',f'XLONG_M',f'SOILTEMP',f'HGT_M',f'SEAICE',f'MAPFAC_MX',f'MAPFAC_MY',
                      f'GREENFRAC',f'GREENFRAC',f'LAI12M',f'LANDMASK',f'LU_INDEX',f'SCT_DOM',f'sd',f'src',
                      f'skt',f'DZS',f'ZS',f'stl',f'swvl'] ## these are a combo of the names of variables for geogrid 
                      ## or era5 data downloaded
    
    HRLDAS_name_variables = [f'Times',f'XLAT',f'XLONG',f'TMN',f'HGT',f'SEAICE',f'MAPFAC_MX',f'MAPFAC_MY',
                            f'SHDMAX',f'SHDMIN',f'LAI',f'XLAND',f'IVGTYP',f'ISLTYP',f'SNOW',f'CANWAT',
                            f'TSK',f'DZS',f'ZS',f'TSLB',f'SMOIS'] ## these are what HRLDAS is looking for (no change)
    
    #These are the units that are saved in the initial conditions
    units = {
                        'XLAT' : 'degrees_north',
                        'XLONG' : 'degrees_east',
                        'TMN' : 'K',
                        'HGT' :'m',
                        'SEAICE' : '-',
                        'MAPFAC_MX' : '-',
                        'MAPFAC_MY' : '-',
                        'SHDMAX' : '%',
                        'SHDMIN' : '%',
                        'LAI' : 'm^2/m^2',
                        'XLAND' : '-',
                        'IVGTYP' : '-',
                        'ISLTYP' : '-',
                        'SNOW' : 'kg/m^2', #mm since a kg of water over a square meter is a depth of a mm
                        'CANWAT' : 'kg/m^2', #mm since a kg of water over a square meter is a depth of a mm
                        'TSK' : 'K', # surface radiative temperature
                        'DZS' : 'm', # depth of the individual layers
                        'ZS' : 'm', # depth of the the entire layer 
                        'TSLB' : 'K', # soil layer temperature
                        'SMOIS' : 'm^3/m^3', # soil layer volumetric soil mositure
                    }
    
    geogrid = xr.open_dataset(geo_file) # load in geofile

    ##grab the geogrid file to create the regridder (only once)
    ds_out = xr.Dataset(
        {
            "lat": (["y","x"],geogrid.XLAT_M[0].values),
            "lon": (["y","x"], geogrid.XLONG_M[0].values),
            "lat_b": (["y_b","x_b"],geogrid.XLAT_C[0].values),
            "lon_b": (["y_b","x_b"],geogrid.XLONG_C[0].values)
        }
    )
    
    ## load in a SINGLE time zarr file that is independent of the time loop below. This will allow us to calculate the 
    ## boundaries of the 
    fs = s3fs.S3FileSystem(anon=True) ### This is needed to be able to access wihtouth erroring
    fmap = s3fs.S3Map(f's3://era5-pds/zarr/2010/07/data/air_temperature_at_2_metres.zarr', s3=fs) ## path to era5 data
    dset_t = xr.open_zarr(fmap, consolidated=True) #grab a single zarr file 
    
    
    ##grab the subset of data
    subset_lon_start,subset_lon_end,subset_lat_start,subset_lat_end = get_region_era5(dset_t,ds_out)



    # Example on how to use the output
    # of era5 region dset_subset_td = dset_td[data_varaibles[2]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]
    _temp =  dset_t['air_temperature_at_2_metres'][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]
    data_array_era5 = transform_era5_to_dataarray_2d(_temp)
                    
    ## now we get the regridder weights (only need this once)
    regridder_era5_to_geogrid = get_regridder(data_array_era5,ds_out)
    
    
    ## We now download the data of interest that we need from the era5 servers
    ## these are all used for 3D mapping!
    data_era5_INIT = get_era5_initial_conditions(start_date)
    
    
    ## now we build the dataset of initial conditions based on the 
    ## the pre-processing steps that are in HRLDAS
    data_set_to_save={}
    
    # load in the geogrid files first, these do not have to have anything reset
    data_set_to_save[HRLDAS_name_variables[0]] = geogrid[data_variables[0]] #Xlat
    data_set_to_save[HRLDAS_name_variables[1]] = geogrid[data_variables[1]] #Xlat
    data_set_to_save[HRLDAS_name_variables[2]] = geogrid[data_variables[2]] #Xlong
    data_set_to_save[HRLDAS_name_variables[3]] = geogrid[data_variables[3]] #TMN
    data_set_to_save[HRLDAS_name_variables[4]] = geogrid[data_variables[4]] #HGT
    
    ## determine if there is sea-ice (likely this would be good to load in another dataset in future
    ## instead of creating a boolean mask by yourself like I am doing here)
    isice_NLCD = xr.where(geogrid.LU_INDEX==22,1,0)
    isice_USGS = xr.where(geogrid.LU_INDEX==15,1,0)
    isice_all = xr.where(isice_NLCD + isice_USGS >=1,1,0)
    
    iswater_USGS = xr.where(geogrid.LU_INDEX==17,1,0)
    iswater_NLCD = xr.where(geogrid.LU_INDEX==21,1,0)
    iswater_all = xr.where(iswater_USGS + iswater_NLCD >=1,1,0)
    
    isice_seaonly = xr.where(iswater_all + isice_all ==2,1,0)
    
    data_set_to_save[HRLDAS_name_variables[5]] = isice_seaonly #
    data_set_to_save[HRLDAS_name_variables[5]].attrs = {'units':units[HRLDAS_name_variables[5]]}
    
    data_set_to_save[HRLDAS_name_variables[6]] = geogrid[data_variables[6]] #mapping coefficientsx
    data_set_to_save[HRLDAS_name_variables[7]] = geogrid[data_variables[7]] #mapping coefficientsy
    data_set_to_save[HRLDAS_name_variables[8]] = geogrid[data_variables[8]].max(axis=1) #Maximum Green Fraction
    data_set_to_save[HRLDAS_name_variables[9]] = geogrid[data_variables[9]].min(axis=1) #Minimum Green Fraction
    
    #LAI needs to be interpolated to the day of year
    var_LAI = get_initial_LAI(geogrid,start_date,data_variables[10])
    
    data_set_to_save[HRLDAS_name_variables[10]] = var_LAI
    data_set_to_save[HRLDAS_name_variables[10]].attrs = {'units':units[HRLDAS_name_variables[10]]}
    
    data_set_to_save[HRLDAS_name_variables[11]] = geogrid[data_variables[11]] # land mask
    data_set_to_save[HRLDAS_name_variables[12]] = geogrid[data_variables[12]] # dominant land use catagory
    data_set_to_save[HRLDAS_name_variables[13]] = geogrid[data_variables[13]] # dominant soil type catagory
    
    ## snow variables (needs to be mm so multiply by 1000)
    ERA5_snow_=  transform_era5_to_dataarray_2d_netcdf(data_era5_INIT[data_variables[14]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end],HRLDAS_name_variables[14]) * 1000
    data_set_to_save[HRLDAS_name_variables[14]] = regridder_era5_to_geogrid(ERA5_snow_)[HRLDAS_name_variables[14]]
    data_set_to_save[HRLDAS_name_variables[14]] = data_set_to_save[HRLDAS_name_variables[14]].reset_coords(drop=True)
    data_set_to_save[HRLDAS_name_variables[14]] = data_set_to_save[HRLDAS_name_variables[14]].swap_dims({'y':'south_north','x':'west_east'})
    data_set_to_save[HRLDAS_name_variables[14]].attrs = {'units':units[HRLDAS_name_variables[14]]}
    data_set_to_save[HRLDAS_name_variables[14]] = data_set_to_save[HRLDAS_name_variables[14]].expand_dims('Time')
    
    
    ## canwat (needs to be mm so multiply by 1000)
    ERA5_canwat_=  transform_era5_to_dataarray_2d_netcdf(data_era5_INIT[data_variables[15]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end],HRLDAS_name_variables[15]) * 1000
    data_set_to_save[HRLDAS_name_variables[15]] = regridder_era5_to_geogrid(ERA5_canwat_)[HRLDAS_name_variables[15]]
    data_set_to_save[HRLDAS_name_variables[15]] = data_set_to_save[HRLDAS_name_variables[15]].reset_coords(drop=True)
    data_set_to_save[HRLDAS_name_variables[15]] = data_set_to_save[HRLDAS_name_variables[15]].swap_dims({'y':'south_north','x':'west_east'})
    data_set_to_save[HRLDAS_name_variables[15]].attrs = {'units':units[HRLDAS_name_variables[15]]}
    data_set_to_save[HRLDAS_name_variables[15]] = data_set_to_save[HRLDAS_name_variables[15]].expand_dims('Time')
    
    ## tsk 
    ERA5_tsk_=  transform_era5_to_dataarray_2d_netcdf(data_era5_INIT[data_variables[16]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end],HRLDAS_name_variables[16]) 
    data_set_to_save[HRLDAS_name_variables[16]] = regridder_era5_to_geogrid(ERA5_tsk_)[HRLDAS_name_variables[16]]
    data_set_to_save[HRLDAS_name_variables[16]] = data_set_to_save[HRLDAS_name_variables[16]].reset_coords(drop=True)
    data_set_to_save[HRLDAS_name_variables[16]] = data_set_to_save[HRLDAS_name_variables[16]].swap_dims({'y':'south_north','x':'west_east'})
    data_set_to_save[HRLDAS_name_variables[16]] = data_set_to_save[HRLDAS_name_variables[16]].expand_dims('Time')
    data_set_to_save[HRLDAS_name_variables[16]].attrs = {'units':units[HRLDAS_name_variables[16]]}
    
    #### The following would need to be adjsuted for finer resolution treatment of soil properties
    data_set_to_save[HRLDAS_name_variables[17]] = xr.DataArray([[0.1,0.3,0.6,1]],dims=['Time','soil_layers_stag'])
    data_set_to_save[HRLDAS_name_variables[17]].attrs = {'units':units[HRLDAS_name_variables[17]]}
    
    data_set_to_save[HRLDAS_name_variables[18]] = xr.DataArray([[0.05,0.25,0.7,1.5]],dims=['Time','soil_layers_stag'])
    data_set_to_save[HRLDAS_name_variables[18]].attrs = {'units':units[HRLDAS_name_variables[18]]}
    
    ## soil temperature  
    soil_temperature_era5 = soil_coarse_stack(data_era5_INIT,data_variables[19])
    ERA5_stemp_ = transform_era5_to_dataarray_3d(soil_temperature_era5[:,:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end],data_variables[19])
    data_set_to_save[HRLDAS_name_variables[19]] = regridder_era5_to_geogrid(ERA5_stemp_)[data_variables[19]]
    data_set_to_save[HRLDAS_name_variables[19]] = data_set_to_save[HRLDAS_name_variables[19]].reset_coords(drop=True)
    data_set_to_save[HRLDAS_name_variables[19]] = data_set_to_save[HRLDAS_name_variables[19]].swap_dims({'time':'Time','y':'south_north','x':'west_east'})
    data_set_to_save[HRLDAS_name_variables[19]].attrs = {'units':units[HRLDAS_name_variables[19]]}
    
     ##  soil moisture 
    soil_mois_era5 = soil_coarse_stack(data_era5_INIT,data_variables[20])
    ERA5_smois_ = transform_era5_to_dataarray_3d(soil_mois_era5[:,:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end],data_variables[20])
    data_set_to_save[HRLDAS_name_variables[20]] = regridder_era5_to_geogrid(ERA5_smois_)[data_variables[20]]
    data_set_to_save[HRLDAS_name_variables[20]] = data_set_to_save[HRLDAS_name_variables[20]].reset_coords(drop=True)
    data_set_to_save[HRLDAS_name_variables[20]] = data_set_to_save[HRLDAS_name_variables[20]].swap_dims({'time':'Time','y':'south_north','x':'west_east'})
    data_set_to_save[HRLDAS_name_variables[20]].attrs = {'units':units[HRLDAS_name_variables[20]]}
    
    
    
    attrs = {'TITLE':'Output from Python Re-grid',
                     'WEST-EAST_GRID_DIMENSION':geogrid.attrs['WEST-EAST_GRID_DIMENSION'],
                     'SOUTH-NORHT_GRID_DIMENSION':geogrid.attrs['SOUTH-NORTH_GRID_DIMENSION'],
                     'DX':geogrid.attrs['DX'],
                     'DY':geogrid.attrs['DY'],
                     'TRUELAT1':geogrid.attrs['TRUELAT1'],
                     'TRUELAT2':geogrid.attrs['TRUELAT2'],
                     'LA1':geogrid.attrs['corner_lats'][0],
                     'LO1':geogrid.attrs['corner_lons'][0],
                     'STAND_LON':geogrid.attrs['STAND_LON'],
                     'MAP_PROJ':geogrid.attrs['MAP_PROJ'],
                     'MMINLU':geogrid.attrs['MMINLU'],
                    }
    print(data_set_to_save)
    data_set_final = xr.Dataset(data_vars=data_set_to_save,attrs=attrs)
    
    if save_location.startswith('s3://'):
           # First write the file locally in the current directory
           data_set_final.to_netcdf(file_name)

           # Now, use aws cli to upload
           # Make sure there isn't an extra slash
           # at the end
           save_location = save_location.strip('/')
           output_path = f'{save_location}/{file_name}'
           os.system(f'aws s3 cp {file_name} {output_path}')

    else:
        data_set_final.to_netcdf(save_location+output_filename)
## end of run block
## basic argument parser given Luke M. 
if __name__ == '__main__':
    # Set up an argument parser object
    from argparse import ArgumentParser
    parser = ArgumentParser()
    
    # For each argument we expect, we need to tell ArgumentParser
    # to expect it

    # Here, add the 'start_date' argument
    # Note that when doing long command line flags (e.g., --start-date),
    # argparse will strip out the leading double-hyphen and convert all other hyphens
    # to underscores in the variable name.  Also, add a help string to define what this argument is.
    parser.add_argument('--start-date', type=str, help="Start date of simulation in format YYYY-MM-DD.  Must be beginning of a month")


    # Argument for save_location
    parser.add_argument('--save-location', type=str, default='./', help="Path to directory location to save outputs")
    
    # Argument for geogrid file
    parser.add_argument('--geogrid-file', type=str, help="Path and name of geogrid file (GEO_EM) from WPS")
    
    # Now, tell argparse to parse command line arguments
    args = parser.parse_args()
    
    # Pass this to our run function
    run(start_date=args.start_date, save_location=args.save_location,geo_file=args.geogrid_file)

    
    