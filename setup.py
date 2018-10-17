import setuptools
import numpy as np
import subprocess
import os.path
from Cython.Build import cythonize

# use f2py to create the python wrapper for Rosato's fortran routine.

# subprocess.check_output('cd ' + str(rosato_database_path), shell=True)
# subprocess.run('mkdir test', shell=True)
try:
    root = os.path.dirname(os.path.realpath(__file__))
    rosato_database_path = os.path.join(root, 'tabulated_data', 'rosato')
    os.chdir(rosato_database_path)
    subprocess.run('f2py -c LS_DATA_read_f2py.f90 -m LS_DATA_read_f2py', shell=True)
    subprocess.check_output(['f2py'])
except Exception as e:
    print('pystark install error:')
    print(e)


setuptools.setup(
    name='pystark',
    version='0.1',
    author='Joseph Allcock',
    description='Calculate Stark-Zeeman-Doppler broadened spectral lineshapes for the hydrogen Balmer series',
    url='https://git.ccfe.ac.uk/jallcock/pystark',
    packages=setuptools.find_packages(),
)


