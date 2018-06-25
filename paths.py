import os

print('hello')

root = os.path.dirname(os.path.realpath(__file__))
netcdf_data_path = os.path.join(root, 'stehle_tables.nc')
raw_data_path = os.path.join(root, 'raw_data')

print(root)
