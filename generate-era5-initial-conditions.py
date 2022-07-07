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
