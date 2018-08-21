# pystark

Calculate Stark / Stark-Zeeman / Stark-Zeeman-Doppler lineshapes for the hydrogen Balmer series.

Authorship:

- James Harrison
- Joseph Allcock


# Functions

- **simple_profile():** Using either Griem's scaling and assuming a Lorentzian lineshape, or else using [B. Lomanowski's
parameterised model](http://iopscience.iop.org/article/10.1088/0029-5515/55/12/123028/meta "Bart's paper") of the tabulated Stehle data.

- **stehle_profile():** Use the tabulated Stehle data (Stark-Doppler)
- **rosato_profile():** Use the tabulated Rosato data (Stark-Zeeman). Will add Doppler broadening to this.


# Prerequisites

- numpy: 1.13.1
- scipy: 0.19.1
- matplotlib: 2.0.2



