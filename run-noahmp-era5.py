import os 
from namelist_utils import read_namelist, write_namelist


def run(ICBC_location,save_location,start_date,start_hour,number_runday,LSM_timestep,output_timestep,geogrid_file):
    
    """
    # assumes that you have already set up the namelist physics options in the namelist.hrldas.draft file
    # can be easily added as a keyword if you need to submit a suite of differetn run options. 
    """
    start_date = f'{start_date[0:4]}{start_date[5:7]}{start_date[8:10]}{start_hour}'
    
    #create output filename
    setup_filename  = f'HRLDAS_setup_{start_date}_d1'
    
    ## copy data to run directory
    os.system(f'aws s3 sync {ICBC_location} /home/jupiter/model/noahmp/hrldas/run/')
#     os.system(f'rm -r /home/jupiter/model/noahmp/hrldas/run/namelist.hrldas')
    
#     ## we first need to adjust the namelist paths
#     v = read_namelist('/home/jupiter/model/noahmp/hrldas/run/namelist.hrldas.draft')
    
#     v['NOAHLSM_OFFLINE']['HRLDAS_SETUP_FILE'] = f'./{setup_filename}'
#     v['NOAHLSM_OFFLINE']['START_YEAR'] = int(start_date[0:4])
#     v['NOAHLSM_OFFLINE']['START_MONTH'] = int(start_date[5:7])
#     v['NOAHLSM_OFFLINE']['START_DAY'] = int(start_date[8:10])
#     v['NOAHLSM_OFFLINE']['START_HOUR'] = int(start_hour)
#     v['NOAHLSM_OFFLINE']['KDAY'] = int(number_runday)
#     v['NOAHLSM_OFFLINE']['NOAH_TIMESTEP'] = int(LSM_timestep)
#     v['NOAHLSM_OFFLINE']['OUTPUT_TIMESTEP'] = int(output_timestep)
#     v['NOAHLSM_OFFLINE']['geogrid_file_name_for_mosaic'] = geogrid_file
    
#     write_namelist(v,"/home/jupiter/model/noahmp/hrldas/run/namelist.hrldas")
    
    #call data
    os.system('mpirun --allow-run-as-root /home/jupiter/model/noahmp/hrldas/run/hrldas.exe > output.log')
    
    ### copy data to save
    os.system(f'aws s3 sync --exclude "*" --include "*LDASOUT_DOMAIN1" . {save_location}')
    
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
    parser.add_argument('--ICBC-location', type=str, help="Path to s3 bucket where ICBC was created")

    # Argument for end date
    parser.add_argument('--save-location', type=str, help="Path to s3 bucket to save output")
    
     # Argument for end date
    parser.add_argument('--start-date', type=str, help="Start of simulation in format YYYY-MM-DD.")
    
    # Argument for end date
    parser.add_argument('--start-hour', type=str, help="End Hour of simulation. Format of HH.  ")
    
    # Argument for end date
    parser.add_argument('--number-runday', type=str, help="Number of days you would like to run, integer")
    
    # Argument for end date
    parser.add_argument('--LSM-timestep', type=str, default='1800',help="Timestep in seconds that the LSM should run. Must divide into 3600")
    
    # Argument for end date
    parser.add_argument('--output-timestep', type=str, default='1800', help="Output timestep in seconds. Must divide into 3600 and be larger then LSM-timestep")
    
    # Argument for geogrid file
    parser.add_argument('--geogrid-file', type=str, help="Path and name of geogrid file (GEO_EM) from WPS")
    
    args = parser.parse_args()
   # Pass this to our run function
    run(ICBC_location=args.ICBC_location,save_location=args.save_location,start_date=args.start_date,
        start_hour=args.start_hour,number_runday=args.number_runday,LSM_timestep=args.LSM_timestep,
        output_timestep=args.output_timestep,geogrid_file=args.geogrid_file)

