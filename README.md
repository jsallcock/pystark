# pystark

Calculate Stark / Stark-Zeeman / Stark-Zeeman-Doppler lineshapes for the hydrogen Balmer series.

Authorship:

- James Harrison
- Joseph Allcock


# Functions

- **simple_profile():** Using either Griem's scaling and assuming a Lorentzian lineshape, or else using [B. Lomanowski's
parameterised model](http://iopscience.iop.org/article/10.1088/0029-5515/55/12/123028/meta "Bart's paper") of the tabulated Stehle data.

- **stehle_profile():** Use the tabulated Stehle data (Stark-Doppler). [paper](https://lerma.obspm.fr/~stehle/Articles/1999AAS140Stehle.pdf)
- **rosato_profile():** Use the tabulated Rosato data (Stark-Zeeman). Will add Doppler broadening to this. [paper](https://www.sciencedirect.com/science/article/pii/S0022407316305325)


# Prerequisites

- numpy
- scipy
- matplotlib



