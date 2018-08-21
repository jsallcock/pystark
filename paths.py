import os
import scipy.io.netcdf as netcdf

root = os.path.dirname(os.path.realpath(__file__))

data_dir_path = os.path.join(root, 'data')

stehle_netcdf_file_path = os.path.join(data_dir_path, 'stehle', 'stehle_tables.nc')
stehle_raw_data_path = os.path.join(data_dir_path, 'stehle', 'raw')

rosato_database_path = os.path.join(data_dir_path, 'rosato', 'database')

# print(stehle_netcdf_file_path)
# load data from 'stehle_tables.nc'
nc = netcdf.netcdf_file(stehle_netcdf_file_path, 'r')
