import s3fs
import xarray as xr
import xesmf
import numpy as np
import pandas as pd
import os

"""
This code assumes that you are going to start at the BEGINNING OF A MONTH to generate the inputs. 
"""

############## options to set examples
#start_date_want = '2010-06-01' #set
#end_date_want = '2010-11-01' #set
#freq_want = '1H' #era5 is 1H
#save_location = './' #save location for data
#geo_file = './geogrid-files/geo_em.d01.milwaukee.nc' #location and name of file
#################


######################## Define all important functions

def cal_specific_humid(Tdew,P):
    """
    calculate the specific humidity in kg/kg from ERA5 inputs. 
    assumes that Tdew is in Kelving and P is in pascals (the raw data)
    
    """
    
    
    e = 6.112*np.exp((17.67*(Tdew-273.15))/((Tdew-273.15)+243.5))
    q = ((0.622*e)/((P/100) - 0.378*e)) #kg per kg
    
    return(q)

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

def transform_era5_to_dataarray_2d(dsubset,name):
    
    """
    takes a zarr dataarray, and creates a dataset that has the correct lats and lons to 
    feed into HRLDAS
    """
    
    lonsv, latsv = np.meshgrid((dsubset.lon.values+180)%360-180,dsubset.lat.values) 
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
    
######################################################################
# begin actual code (run function)

def run(start_date, end_date, freq_want, save_location, geo_file):
    
    
    
    ## create a list of inputs to loop over (data is ERA5, NAME is the names HRLDAS is looking for)
    data_variables = [f'air_temperature_at_2_metres',f'surface_air_pressure',f'dew_point_temperature_at_2_metres',
                     f'eastward_wind_at_10_metres',f'northward_wind_at_10_metres',
                     f'integral_wrt_time_of_surface_direct_downwelling_shortwave_flux_in_air_1hour_Accumulation',
                     f'precipitation_amount_1hour_Accumulation',f'strd'] ## these are what ERA5 is defining (zarr format)
    
    name_variables = [f'T2D',f'PSFC',f'Q2D',f'U2D',f'V2D',f'SWDOWN',f'RAINRATE','LWDOWN'] ## these are what HRLDAS is looking for (no change)
    
    variable_name_longwave = 'surface_thermal_radiation_downwards' ## this needs to be used to access the ERA5 data that was pre-downloaded

    # fs = s3fs.S3FileSystem() ### This is needed to be able to access wihtouth erroring
    # fmap = fs.open(geo_file, mode='rb') ## path to era5 data
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
    fmap = s3fs.S3Map(f's3://era5-pds/zarr/2010/07/data/{data_variables[0]}.zarr', s3=fs) ## path to era5 data
    dset_t = xr.open_zarr(fmap, consolidated=True) #grab a single zarr file 
    
    
    ##grab the subset of data
    subset_lon_start,subset_lon_end,subset_lat_start,subset_lat_end = get_region_era5(dset_t,ds_out)



    # Example on how to use the output
    # of era5 region dset_subset_td = dset_td[data_varaibles[2]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]
    _temp =  dset_t[data_variables[0]][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]
    data_array_era5 = transform_era5_to_dataarray_2d(_temp,name_variables[0])

    ## now we get the regridder weights (only need this once)
    regridder_era5_to_geogrid = get_regridder(data_array_era5,ds_out)
    regridder_era5_to_geogrid_conserve = get_regridder(data_array_era5,ds_out,"conservative")

    ## works up to here!
    dates_to_loop = pd.date_range(start=start_date,end=end_date, freq=freq_want)
    date_year = str(dates_to_loop[0].year)
    date_month = str(dates_to_loop[0].month).zfill(2)

    variables_to_save = {} # create an empty dictionary 

    for ii, date in enumerate(dates_to_loop): # we loop over each of these values

        # naming convention is of the format
        # YYYYMMDDHH.LDASIN_DOMAIN1

        output_file_name = str(date.year).zfill(4) + str(date.month).zfill(2) + str(date.day).zfill(2) + str(date.hour).zfill(2) + '.LDASIN_DOMAIN1'

        ## increment the longwave data

        date_month_new = str(date.month).zfill(2)

        ## determine if we need to load in a new data set
        if ii == 0:
            load_new = True
        elif date_month_new == date_month:
            load_new = False
        else:
            date_month = date_month_new
            load_new = True
            print(date_month)

        if (date.month % 3  == 1) or (ii == 0):
            load_LWDOWN = True
        else:
            load_LWDOWN = False


        ## if we should load the data  LOAD IT
        if load_new:
            for era5_name, output_name in zip(data_variables, name_variables):

                # do we need to load in longwave data? It is collected every three months
                if (output_name == 'LWDOWN') and (load_LWDOWN): #LWDOWN branch
                    
                    # which ending do we need?
                    if (date.month>= 1) and (date.month<=3):
                        end_string = f'_01-03.nc'
                    elif(date.month>= 4) and (date.month<=6):
                        end_string = f'_04-06.nc'
                    elif(date.month>= 7) and (date.month<=9):
                         end_string = f'_07-09.nc'
                    else:
                        end_string = f'_10-12.nc'

                    # create the file that we need to load in
                    file_name_LW = f'{variable_name_longwave}_'+date_year+end_string

                    # actually load it in (ERA5 SPECIFIC)
                    os.system(f'aws s3 cp s3://jupiter-reference-data/era5/rlds/{file_name_LW} .')
                    dLW = xr.open_dataset(file_name_LW)[era5_name]
                    dLW = dLW.rename({'longitude':'lon','latitude':'lat','time':'time'})
                    
                    # We divide by 3600 to get the W/m^2
                    dset_loop_sm = dLW[:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]/3600
                    dset_data_array = transform_era5_to_dataarray_3d(dset_loop_sm,output_name)

                    variables_to_save[output_name] = regridder_era5_to_geogrid(dset_data_array[output_name])#ensures this is a dataarray
                    
                    ## meta data and fixing the attributes
                    variables_to_save[output_name].attrs = {'units':'W/m^2'}
                    variables_to_save[output_name] = variables_to_save[output_name].reset_coords(drop=True)
                    variables_to_save[output_name] = variables_to_save[output_name].swap_dims({'time':'Time','y':'south_north','x':'west_east'})
                    

                else:
                    if (output_name == 'LWDOWN'):
                        continue
    
                    ## load in the data (ERA5 SPECIFIC)
                    fs = s3fs.S3FileSystem(anon=True) ### This is needed to be able to access wihtouth erroring
                    fmap = s3fs.S3Map(f's3://era5-pds/zarr/{date_year}/{date_month}/data/{era5_name}.zarr', s3=fs)
                    dset_loop = xr.open_zarr(fmap, consolidated=True)


                    # if rain or if solar radiation, we divide by 3600. if anything else, we divide by 1 (keep it the same)
                    if (output_name == 'SWDOWN'):
                        divisor = 3600
                    elif (output_name == 'RAINRATE'):
                        divisor = 3600/1000 #rainfall is in m and this is a divisor, so flip the 1000/3600 that we would expect
                    else:
                        divisor = 1

                    # subset the data
                    dset_loop_sm = dset_loop[era5_name][:,subset_lat_start:subset_lat_end,subset_lon_start:subset_lon_end]/divisor
                    dset_data_array = transform_era5_to_dataarray_3d(dset_loop_sm,output_name)
                   

                    ## regrid
                    ## do we need to conserve the data? 
                    if (output_name == 'RAINRATE'):
                        variables_to_save[output_name] = regridder_era5_to_geogrid_conserve(dset_data_array[output_name])#makes this a dataarray
                    else:
                        variables_to_save[output_name] = regridder_era5_to_geogrid(dset_data_array[output_name])#makes this a dataarray
                    
                    ## meta data and fixing the attributes
                    units = {
                        'T2D' : 'K',
                        'PSFC' : 'Pa',
                        'Q2D' : 'Kg/kg',
                        'U2D' :'m/s',
                        'V2D' : 'm/s',
                        'SWDOWN' : 'W/m^2',
                        'RAINRATE' : 'mm/s',
                    }
                    variables_to_save[output_name].attrs = {'units' : units[output_name]}

                        
                    variables_to_save[output_name] = variables_to_save[output_name].reset_coords(drop=True)
                    variables_to_save[output_name] = variables_to_save[output_name].swap_dims({'time':'Time','y':'south_north','x':'west_east'})
                        
                    

            # get the specific humidity, and not dew point temp
            variables_to_save['Q2D'] = cal_specific_humid(variables_to_save['Q2D'],variables_to_save['PSFC'])
            ## end of the for loop to load in data if needed


        #how do we index the ZARR data?
        if ii == 0:
            index_not_LW = ii
        elif load_new: #we need to restart!
            index_not_LW = 0
        else:
            index_not_LW += 1


        if ii == 0: # we should start out at the correct time 
            if date.month == 1 or date.month == 4 or date.month == 7 or date.month == 10:
                index_LW = 0
            elif date.month == 2 or date.month == 5 or date.month == 8 or date.month == 11:
                index_LW = 24*date.days_in_month - 1 #minus 1 because of the zero index
            else:
                x = date-pd.DateOffset(months=1)
                index_LW = 24*(date.days_in_month + x.days_in_month)- 1 #minus 1 because of the zero index
        elif load_LWDOWN:
            index_LW = 0 # reset with the new netcdf
        else:
            index_LW += 1


        data_set_save = xr.Dataset(
            {
                name_variables[0]:variables_to_save[name_variables[0]].isel(Time=[index_not_LW]),
                name_variables[1]:variables_to_save[name_variables[1]].isel(Time=[index_not_LW]),
                name_variables[2]:variables_to_save[name_variables[2]].isel(Time=[index_not_LW]),
                name_variables[3]:variables_to_save[name_variables[3]].isel(Time=[index_not_LW]),
                name_variables[4]:variables_to_save[name_variables[4]].isel(Time=[index_not_LW]),
                name_variables[5]:variables_to_save[name_variables[5]].isel(Time=[index_not_LW]),
                name_variables[6]:variables_to_save[name_variables[6]].isel(Time=[index_not_LW]),
                name_variables[7]:variables_to_save[name_variables[7]].isel(Time=[index_LW]),
            },
            attrs = {'TITLE':'Output from Python Re-grid',
                     'WEST-EAST_GRID_DIMENSION':geogrid.attrs['WEST-EAST_GRID_DIMENSION'],
                     'SOUTH-NORTH_GRID_DIMENSION':geogrid.attrs['SOUTH-NORTH_GRID_DIMENSION'],
                     'DX':geogrid.attrs['DX'],
                     'DY':geogrid.attrs['DY'],
                     'TRUELAT1':geogrid.attrs['TRUELAT1'],
                     'TRUELAT2':geogrid.attrs['TRUELAT2'],
                     'LA1':geogrid.attrs['corner_lats'][0],
                     'LO1':geogrid.attrs['corner_lons'][0],
                     'STAND_LON':geogrid.attrs['STAND_LON'],
                     'MAP_PROJ':geogrid.attrs['MAP_PROJ'],
                     'GRID_ID':geogrid.attrs['grid_id'],
                     'ISWATER':geogrid.attrs['ISWATER'],
                     'ISURBAN':geogrid.attrs['ISURBAN'],
                     'ISICE':geogrid.attrs['ISICE'],
                     'MMINLU':geogrid.attrs['MMINLU'],
            }


        )

        if save_location.startswith('s3://'):
           # First write the file locally in the current directory
           data_set_save.to_netcdf(output_file_name)

           # Now, use aws cli to upload
           # Make sure there isn't an extra slash
           # at the end
           save_location = save_location.strip('/')
           output_path = f'{save_location}/{output_file_name}'
           os.system(f'aws s3 cp {output_file_name} {output_path}')

        else:
            data_set_save.to_netcdf(f'{save_location}{output_file_name}')
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

    # Argument for end date
    parser.add_argument('--end-date', type=str, help="End date of simulation in format YYYY-MM-DD.  Must be beginning of a month")

    # Argument for frequency
    # Here we add a default of '1H', which will be assigned
    # if the user does not provide any argument to this flag
    # (could add defaults for the others too if you wanted)
    parser.add_argument('--freq', type=str, default='1H', help="Frequency of boundary condition. ERA5 defaults to 1H")

    # Argument for save_location
    parser.add_argument('--save-location', type=str, default='./', help="Path to directory location to save outputs")
    
    # Argument for geogrid file
    parser.add_argument('--geogrid-file', type=str, help="Path and name of geogrid file (GEO_EM) from WPS")
    
    # Now, tell argparse to parse command line arguments
    args = parser.parse_args()
    
    # Pass this to our run function
    run(start_date=args.start_date, end_date=args.end_date, freq_want=args.freq, save_location=args.save_location,geo_file=args.geogrid_file)
