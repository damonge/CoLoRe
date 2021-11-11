# CoLoRe

CoLoRe is a parallelized code to generate fast 3D realization of a wide variety of cosmological observations.

## 1 Methods

The methods used by CoLoRe are described in [Ramirez-Perez et al. 2021](https://arxiv.org/abs/2111.05069).

When in doubt, bear in mind that by default CoLoRe uses the following units:
 - Lenghts: Mpc/h
 - Angles: degrees


## 2 Compilation and usage

To compile CoLoRe, open the Makefile and edit it according to your
system. The default options (except for the paths to the external
libraries) should work for most systems.

### 2.1 Compilation options
The following compilation flags, tunable in the Makefile, can be modified
to control the behaviour of CoLoRe:

- CoLoRe may be very memory-demanding. To minimize the memory overhead use single precision floating point (`USE_PRECISION = yes`).
- OpenMP parallelization is enabled by setting the option `USE_OMP` to "yes".
- MPI parallelization is enabled by setting the option `USE_MPI` to "yes".
- The default bias model in CoLoRe is exponential. To select a different bias model add `-D_BIAS_MODEL_2` (exp-truncated) or `-D_BIAS_MODEL_3` (truncated).

### 2.2 Dependencies
CoLoRe uses 4 external packages:
 - GSL. The GNU Scientific Library (tested for versions 3.*)
 - FFTW. The Fastest Fourier Transform of the West (versions 3.*)
 - CFITSIO. FITS format library. This package is optional.
 - HDF5. HDF5 format library. This package is optional.
The paths to the corresponding headers and libraries should be correctly
set in the Makefile.

### 2.3 Running the code
Once the Makefile has been editted, typing 'make' should generate
the executable 'CoLoRe'. To run CoLoRe just type

> mpirun -np <number-of-nodes> ./CoLoRe <param_file>

where <param_file> is the path to the parameter file described in
section 3.


## 3 Parameter file and examples.

The behaviour of CoLoRe is mainly controlled by the input param file. The
param file is basically a set of name-value pairs. Any blank lines, and
anything beyond a #-symbol will be ignored.   We provide a [sample param
file](param_example.cfg) that includes all the input parameters
needed by CoLoRe. The comments included in this file explain the meaning
and functionality of these parameters.

We also provide an ipython [notebook](example_CoLoRe.ipynb) that demonstrates
how to generate all the different probes implemented in CoLoRe. The notebook
also exemplifies how to interpret the different outputs.


## 4 License

CoLoRe is distributed under the GPL license (see COPYING in the root
directory). We kindly ask you to cite the companion paper
[Ramirez-Perez et al. 2021](TBD) when using the code.


## 5 Contact

Regarding bugs, suggestions, questions or petitions, feel free to contact
the authors:
    David Alonso: david.alonso@physics.ox.ac.uk
    Cesar Ramirez: cramirez@ifae.es
