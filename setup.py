import setuptools
import subprocess
import os.path


# try using f2py to create the python wrapper for Rosato's fortran routine.
try:
    root = os.path.dirname(os.path.realpath(__file__))
    rosato_database_path = os.path.join(root, 'tabulated_data', 'rosato')
    os.chdir(rosato_database_path)
    subprocess.run('f2py -c LS_DATA_read_f2py.f90 -m LS_DATA_read_f2py', shell=True)
    # subprocess.check_output(['f2py'])

    os.chdir(root)
    ### For some reason, it is important to change location back to the root directory, or the package is
    # looked for in the wrong place when importing.

except Exception as e:
    print('pystark install error:')
    print(e)




setuptools.setup(
    name='pystark',
    version='0.1',
    author='Joseph Allcock, James Harrison',
    description='Calculate Stark-Zeeman-Doppler broadened spectral lineshapes for the hydrogen Balmer series',
    url='https://git.ccfe.ac.uk/jallcock/pystark',
    packages=setuptools.find_packages(),
)


