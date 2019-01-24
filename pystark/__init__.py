import os
import scipy.io.netcdf as netcdf

# define paths
root = os.path.abspath(os.path.join(os.path.realpath(__file__), '..', '..'))
data_dir_path = os.path.join(root, 'tabulated_data')
stehle_netcdf_file_path = os.path.join(data_dir_path, 'stehle', 'stehle_tables.nc')
stehle_raw_data_path = os.path.join(data_dir_path, 'stehle', 'raw')
rosato_path = os.path.join(data_dir_path, 'rosato/')
rosato_database_path = os.path.join(rosato_path, 'database')

# load tabulated_data from 'stehle_tables.nc'
nc = netcdf.netcdf_file(stehle_netcdf_file_path, 'r')

from . balmer_lineshape import *

try:
    from . rosato_wrapper import *
except:
    print('--pystark: Rosato model could not be imported')

from . tools import *
from . demo import demo





