# pystark

Calculate Stark-Zeeman-Doppler broadened spectral lineshapes for the hydrogen Balmer series. Four lineshape models are available.

Authorship:

- James Harrison
- Joseph Allcock (joseph.allcock@durham.ac.uk for contact)

### Prerequisites

- numpy
- scipy
- matplotlib

### Setup

In the top-level cloned pystark directory run:

```
pip install -e .
```

from the terminal. Let me know if there are any problems.



### Functions

- **pystark.demo():** Run this to see what is and isn't working. Demonstrates how to use each of the lineshape models, produces an example plot. 

- **pystark.simple_profile():** Using either Griem's scaling and assuming a Lorentzian lineshape, or else using [B. Lomanowski's
parameterised model](http://iopscience.iop.org/article/10.1088/0029-5515/55/12/123028/meta "Bart's paper") of the tabulated Stehle data.

- **pystark.stehle_profile():** Use the tabulated Stehle data (Stark-Doppler). [paper](https://lerma.obspm.fr/~stehle/Articles/1999AAS140Stehle.pdf)
- **pystark.rosato_profile():** Manually convolves the Rosato tabulated profiles (Stark-Zeeman) with the Doppler profile. [paper](https://www.sciencedirect.com/science/article/pii/S0022407316305325). This one is largely a wrapper for some fortran subroutines which do the interpolation. Compilation is (in theory) handled by setup.py using 'f2py'. It should give you a clear error message if something goes wrong.



