import os
import scipy.io.netcdf as netcdf

root = os.path.dirname(os.path.realpath(__file__))
stehle_netcdf_file_path = os.path.join(root, 'stehle_data', 'stehle_tables.nc')
stehle_raw_data_path = os.path.join(root, 'stehle_data', 'raw')

# load data from 'stehle_tables.nc'
nc = netcdf.netcdf_file(stehle_netcdf_file_path, 'r')
