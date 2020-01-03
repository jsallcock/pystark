import setuptools
import subprocess
import os.path


# try using f2py to create the python wrapper for Rosato's fortran routine.
root = os.path.dirname(os.path.realpath(__file__))
rosato_database_path = os.path.join(root, 'tabulated_data', 'rosato')
os.chdir(rosato_database_path)
# f2pys = ['f2py3.5', 'f2py3.6', 'f2py']
f2pys = ['f2py', 'f2py3.6']
for f2py in f2pys:
    try:
        subprocess.run(f2py + ' -c -m rosato_f90_funcs LS_DATA_read_f2py.f90', shell=True)
        subprocess.check_output(['f2py'])
        break
    except Exception as e:
        print('pystark install error:')
        print(e)

os.chdir(root)
### it is important to change location back to the root directory, or the package is
# looked for in the wrong place when importing.

setuptools.setup(
    name='pystark',
    version='0.1',
    author='Joseph Allcock, James Harrison',
    description='Calculate Stark-Zeeman-Doppler broadened spectral lineshapes for the hydrogen Balmer series',
    url='https://git.ccfe.ac.uk/jallcock/pystark',
    packages=setuptools.find_packages(),
)


