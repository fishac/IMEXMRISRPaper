**Implicit-Explicit Multirate Infinitesimal Stage-Restart Methods**

This is the public repository for the paper _Implicit-Explicit Multirate Infinitesimal Stage-Restart Methods_ by Alex Fish, Daniel Reynolds, and Steven Roberts.

# Requirements

This code base was written in C++11, Python 3, and Matlab.
We used the `g++` compiler for C++ compilation.

C++ was used for the implementation of the numerical methods and controllers, and requires the usage of the [Armadillo](http://arma.sourceforge.net/) library.
Python was used for postprocessing and generating plots, and require the [Matplotlib](https://matplotlib.org/) and [Pandas](https://pandas.pydata.org/) libraries.
Matlab was used for generating "true" solutions at fixed points, solutions accurate to strict tolerances, when analytical solutions were not available.
Bash was used for simplifying command-line instructions.

# Common Functionality

Nonlinear solvers and various other common functionality is found in `src/common`.

# Multirate Methods

Code for taking multirate method steps and adaptively evolving solutions is found in `src/methods/MRI`. 
Individual methods are found in `src/methods/MRI/instances`.

# Singlerate (DIRK/ERK) Methods

Code for taking singlerate (DIRK/ERK) method steps and adaptively evolving solutions is found in `src/methods/DIRK`.
Individual methods are found in `src/methods/DIRK/instances`.

# Controllers

Code for controller behavior is found in `src/controllers`.

# Problems

Problems with various splittings are found in `src/problems`.

# Resources

True solution data, generated from Matlab's `ode15s` is stored in `src/resources`.

# Output

Output from the drivers is stored to disc in `src/output`.

# Drivers

Main functions are found in `src/drivers`.

# Postprocessing

Plotting and other data postprocessing is found in `src/postprocessing/output`.
